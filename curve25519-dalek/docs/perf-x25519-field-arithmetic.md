# X25519 field arithmetic: ADX/BMI2 on x86_64 and the wasm32 backend

This document records a measurement pass over the field arithmetic on the
X25519 path, on the two targets a downstream consumer of this fork uses:
**x86_64 with ADX/BMI2** and **wasm32**. It covers what was measured, what was
changed, and — at least as importantly — the two fronts that were measured and
*not* changed, with the numbers that closed them.

Everything here was produced inside this repository. The external profile that
motivated the work is quoted only as motivation; none of its numbers are
reproduced or relied on.

---

## 1. Motivation and what the code actually does

The work started from a profile of a WhatsApp protocol client using this crate
through `wacore-libsignal` (ping-pong messaging, `perf` on real Zen 4 hardware
with `avx2`/`sha_ni`/`adx`/`bmi2`, two X25519 operations per message):

| symbol | self |
| --- | ---: |
| `serial::u64::field::FieldElement51::pow2k` | 30.8% |
| `<serial::u64::field::FieldElement51 as Mul>::mul` | 25.3% |
| `RustCryptoProvider::x25519_agreement` | 7.7% |
| `sha2::sha256::x86_sha::compress` | 3.3% |

That is ~56% of the client in serial field arithmetic, in a binary that had the
AVX2 backend compiled in. The claim to be checked was that this crate's AVX2
backend does not reach X25519 at all.

### Confirmed: no vectorized path reaches X25519

* `src/field.rs` resolves `FieldElement` through a `cfg_if!` whose branches are
  `serial::u64::FieldElement51`, `fiat_u64::FieldElement51`,
  `serial::u32::FieldElement2625` and `fiat_u32::FieldElement2625`. **No branch
  resolves to a vectorized type.**
* `src/montgomery.rs` uses `crate::field::FieldElement`, so it inherits that.
* `grep -rn montgomery src/backend/vector/` returns nothing.
* The vector backend is wired into `backend.rs` only for
  `pippenger`, `straus`, `precomputed_straus` and `variable_base` — i.e.
  Edwards/Ristretto variable-base and multiscalar multiplication.
* `MontgomeryPoint::mul_base` goes through `EdwardsPoint::mul_base`, which uses
  the **serial** `EdwardsBasepointTable`, not the vector backend.

Measured confirmation of the last two points: forcing
`curve25519_dalek_backend="serial"` on a machine where the default is `"simd"`
changes neither X25519 operation (see §4.1, `serial` row) —
`mul_clamped` 59.97 µs vs 59.84 µs and `mul_base_clamped` 18.81 µs vs 18.77 µs,
both inside run-to-run noise.

And the other half of the claim, measured rather than assumed — the vector
backend *is* live in this build and does deliver where it applies.
`EdwardsPoint::vartime_double_scalar_mul_basepoint`, the
signature-verification shape and one of the operations the vector backend does
cover:

| backend | `vartime_double_scalar_mul_basepoint` |
| --- | ---: |
| `"simd"` (AVX2, the default here) | **38 491 ns** |
| `"serial"` | 47 949 ns (+24.6%) |

**Conclusion: on x86_64, the AVX2 backend is working — it is worth ~20% on
ed25519/ristretto variable-base and multiscalar work — and it does nothing at
all for either X25519 operation.** A binary with AVX2 compiled in and 531 AVX2
instructions present, spending its time in `serial::u64`, is the expected
behaviour rather than a misconfiguration: runtime detection succeeds, the
vectorized code runs, and the Montgomery ladder simply never calls it.

### Field operations per `mul_clamped`

Counted from `src/montgomery.rs`. One `differential_add_and_double`
(`montgomery.rs:428`) performs:

* squarings: `t4`, `t5`, `t11`, `t12` → **4**
* multiplications: `t7`, `t8`, `t13` (`APLUS2_OVER_FOUR * t6`), `t14`, `t16`,
  `t17` → **6**
* additions/subtractions: `t0`, `t1`, `t2`, `t3`, `t6`, `t9`, `t10`, `t15` → 8

(`t13` multiplies by the constant 121666. It was a full field multiplication
when this was counted; §5 replaces it with a specialized one costing about a
fifth as much. The tables in this section describe the code as it was at the
start of the work.)

`mul_clamped` → `Mul<&Scalar> for &MontgomeryPoint` → `mul_bits_be(bits_le().rev().skip(1))`
= **255** steps. The tail is `ProjectivePoint::as_affine`, which is one `invert`
plus one `mul`; `invert` is `pow22501` (249 squarings via `pow2k`, 10 muls —
`self * &t1` and nine `&t_i * &t_j`) plus `pow2k(5)` and one mul.

| | multiplications | squarings |
| --- | ---: | ---: |
| ladder, 255 × `differential_add_and_double` | 1530 | 1020 |
| final `invert` + `as_affine` | 12 | 254 |
| **total per `mul_clamped`** | **1542** | **1274** |

These are exposed as constants in `src/bench_internals.rs` so the benchmarks can
use them.

**Cross-check against measurement.** Using the pre-change numbers from §4.1
(default flags): 1530 × 26.11 ns + 1020 × 22.51 ns = 62.9 µs for the ladder, plus
254 × 14.37 ns + 12 × 26.11 ns = 3.9 µs for the inversion, giving a predicted
66.9 µs against a measured 59.8 µs. The model is 12% high, which is the expected
sign and magnitude: within one ladder step the two squarings `t4`/`t5` are
independent of each other, as are the two multiplications `t7`/`t8`, so the real
step is faster than a strictly serial sum of latencies. The op count is
therefore corroborated, and it is what turns "25% of the client in `mul`" into
"1542 `mul` per X25519".

---

## 2. Machine, toolchain, and method

**Host** (a shared 4-vCPU KVM guest, *not* the Zen 4 of the motivating profile):

```text
Architecture:   x86_64
Vendor ID:      GenuineIntel
Model name:     Intel(R) Xeon(R) Processor @ 2.80GHz
CPU family:     6   Model: 85   Stepping: 7      (Cascade Lake / Skylake-SP)
CPU(s):         4   (KVM guest, hypervisor: KVM)
L1d: 128 KiB (4 inst.)  L2: 4 MiB (4 inst.)  L3: 33 MiB (1 inst.)
Relevant flags: avx2 bmi1 bmi2 adx fma aes
                avx512f avx512dq avx512cd avx512bw avx512vl avx512_vnni
                (no avx512ifma, so the AVX-512 IFMA backend is not selected;
                 no sha_ni)
```

**Toolchain:** `rustc 1.94.1 (e408947bf 2026-03-25)`, `cargo 1.94.1`,
`node v22.22.2`, target `wasm32-unknown-unknown` (rust-std installed via
rustup).

### Noise control

The host is a shared virtual machine and Criterion's bootstrapped **mean** was
not usable on it: the very first baseline run produced 95% intervals like
`[67.4 µs, 75.8 µs]` (±12%), and one synthetic bench came back `[232 µs, 386 µs]`
— a 60% spread. Steal time was intermittently non-zero.

So the primary instrument is a purpose-built driver,
`benches/harness/`, which reports the **minimum** over repetitions rather than
the mean. Interference can only ever make a repetition slower, so the minimum is
the closest available estimate of the noise-free cost; the median and the
min→max spread are printed next to it so the reader can judge how quiet the host
was. Every run below is `--reps 15` on a calibrated iteration count targeting
~20 ms per repetition, pinned with `taskset -c 2`.

Where a number carries weight it is also confirmed with Criterion, which is a
completely independent harness (§4.3).

All kernels are **dependent chains** — the output of each operation feeds the
next. The Montgomery ladder is latency bound, so a chained kernel is the honest
model; it also removes the need for `black_box`, which does not exist in the
same form on wasm32.

### What was built

| file | purpose |
| --- | --- |
| `benches/x25519_field.rs` | Criterion benches: `mul_clamped`, `mul_base_clamped`, and (behind `--cfg curve25519_dalek_bench_internals`) `mul`, `square`, `square2`, `add`, `sub`, `pow2k(k)` for k ∈ {1,2,4,10,50,100}, `invert`, and a synthetic "field mix of one `mul_clamped`" |
| `benches/harness/` | standalone crate with one set of kernels compiled for **both** x86_64 and wasm32; `src/main.rs` is the native driver, `run.mjs` the Node driver. Kernels: `fe_mul`, `fe_square`, `fe_pow2k50`, `fe_mul121666`, `fe_invert`, `edwards_mul_base` (and radix-32/64 variants), `edwards_to_montgomery`, `edwards_vartime_double_base`, `x25519_mul_clamped`, `x25519_mul_base_clamped`, and, added for §9, `ed25519_verify`, `edwards_table_create`, and `fe_invert_x{8,16}` against `fe_batch_invert_{8,16}`. The verification kernel makes `ed25519-dalek` a path dependency **of the harness only** — it is not a dependency of `curve25519-dalek`. `src/bin/cg.rs` is the `callgrind` entry point for §13.6. |
| `src/bench_internals.rs` | `--cfg`-gated, `#[doc(hidden)]` hook exposing `FieldElement` and the op counts to the benches. `FieldElement` is `pub(crate)` and Criterion benches are separate crates, so without this the only way to time `mul`/`square` is through a whole scalar multiplication — the exact confound these benches exist to avoid. Absent from every ordinary build; the public API is unchanged. |

Reproduce with:

```sh
# Criterion, field internals included
RUSTFLAGS='--cfg curve25519_dalek_bench_internals' \
    cargo bench --bench x25519_field

# min-of-repetitions driver, x86_64
cd curve25519-dalek/benches/harness
RUSTFLAGS='--cfg curve25519_dalek_bench_internals' cargo run --release --bin native -- --reps 15

# min-of-repetitions driver, wasm32
RUSTFLAGS='--cfg curve25519_dalek_bench_internals' \
    cargo build --release --target wasm32-unknown-unknown --lib
node run.mjs target/wasm32-unknown-unknown/release/x25519_field_harness.wasm 15
```

---

## 3. Front A — x86_64 with ADX/BMI2: **refused**

### 3.1 The gate: does LLVM already emit `mulx`/`adcx`/`adox`?

Counted before writing any code, from `--emit=asm` over the whole crate, then
per symbol for `<&FieldElement51 as Mul>::mul` and `FieldElement51::pow2k`.

Whole-crate totals:

| RUSTFLAGS | `mulx` | `adcx` | `adox` |
| --- | ---: | ---: | ---: |
| *(none)* | 0 | 0 | 0 |
| `-C target-feature=+adx,+bmi2` | 285 | **0** | **0** |
| `-C target-cpu=native` | 285 | **0** | **0** |

Per function:

| flags | fn | instrs | `mulx` | `mulq` | `adcx` | `adox` | `adc` | `add` |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| *(none)* | `mul` | 221 | 0 | 25 | 0 | 0 | 24 | 27 |
| *(none)* | `pow2k` | 152 | 0 | 15 | 0 | 0 | 14 | 21 |
| `+adx,+bmi2` | `mul` | 188 | **25** | 0 | **0** | **0** | 24 | 27 |
| `+adx,+bmi2` | `pow2k` | 129 | **15** | 0 | **0** | **0** | 14 | 21 |
| `target-cpu=native` | `mul` | 187 | 25 | 0 | 0 | 0 | 24 | 27 |
| `target-cpu=native` | `pow2k` | 129 | 15 | 0 | 0 | 0 | 14 | 21 |

Two readings, and both matter:

* **The multiply side is already optimal.** With `+bmi2`, LLVM emits exactly 25
  `mulx` in `mul` (= the 25 partial products of the 5×5 schoolbook) and exactly
  15 in `pow2k` (= the 15 distinct products of a 5-limb squaring). There is
  nothing left for `mulx` to win.
* **The carry side uses no ADX at all**: zero `adcx`, zero `adox`, in every
  configuration, replaced by 24 plain `adc` on a single flag chain.

So the gate as originally stated ("Front A dies if LLVM already emits
`mulx`/`adcx`/`adox`") came back split. The remaining question was whether
hand-written ADX could beat LLVM's `mulx` + single-`adc` output.

### 3.2 Why ADX has nothing to interleave in a 5×51 layout

`adcx`/`adox` exist to run **two independent carry chains** — one on CF, one on
OF — through a *saturated* multi-word accumulator, so that the flag dependency
between consecutive add-with-carry instructions stops serializing the
multiply-accumulate. That is the shape of the 4×64 schoolbook, and it is the
representation BoringSSL's `x25519_x86_64`, libsodium's ADX path and
fiat-crypto's ADX output all use to get the 20–40% these implementations are
cited for.

This crate's `FieldElement51` is **5 unsaturated 51-bit limbs**. Each output
limb is the sum of five ~102-bit partial products in a *two-word* `u128`
accumulator, which the asm above shows plainly: 27 `add` (low words) and 24
`adc` (high words) — one add/adc pair per product, exactly two words deep. A
two-word accumulator has one carry to propagate, so there is no second chain for
`adox` to carry; the dual-chain trick has nothing to interleave. The 5×51 layout
was chosen precisely so that carries stay cheap and rare, and that choice is
what removes ADX's advantage.

Getting the cited 20–40% would mean adopting the 4×64 saturated representation,
i.e. replacing the limb layout and every documented `bit excess` precondition
with a different backend. **That is explicitly out of scope for this work**, and
it is a new backend rather than a change to `mul`.

### 3.3 The measured caution: isolated `mul` gains need not transfer

Before refusing, one experiment tested whether a large isolated `mul`
improvement would even show up end to end. Adding `#[inline]` to
`<&FieldElement51 as Mul>::mul`:

| | isolated `fe_mul` | `mul_clamped` |
| --- | ---: | ---: |
| without `#[inline]` | 25.65 ns | 46.94 µs |
| with `#[inline]` | **18.08 ns (−29.5%)** | 46.71 µs (**−0.5%**) |
| without, `+adx,+bmi2` | 22.85 ns | 40.73 µs |
| with, `+adx,+bmi2` | **16.36 ns (−28.4%)** | 41.01 µs (**+0.7%**) |

A 30% isolated win produced nothing end to end (the isolated gain is largely
LLVM hoisting the loop-invariant `b[i]*19` precomputation out of the chained
kernel, which the ladder cannot do). The change was reverted. It is recorded
here because it is the calibration for how much an `asm!` rewrite would have to
deliver before it was worth its cost — and because the `mul_clamped` column is
the one that counts.

For balance: real `mul` improvements *do* transfer. `+bmi2` alone moves isolated
`fe_mul` by −10.9% and `mul_clamped` by −10.2% (§4.1).

### 3.4 Verdict on Front A

**Refused.** No `asm!`, no `core::arch` intrinsics, no new `unsafe`.

* LLVM already emits every available `mulx` (25/25 in `mul`, 15/15 in `pow2k`)
  under `+bmi2`.
* `adcx`/`adox` have no second carry chain to run in a 5×51 two-word
  accumulator; their payoff belongs to the 4×64 saturated layout, which would be
  a representation change and a new backend.
* An intervention that made isolated `mul` 30% faster moved `mul_clamped` by
  0.5%.

### 3.5 What Front A *did* produce, at zero code cost

Building with BMI2 enabled is worth **10–13% of X25519** with no source change
at all. On the default `x86-64` baseline ISA the crate gets `mulq`; the
consumer only gets `mulx` if they ask for it:

| build flags | `mul_clamped`, before this PR | after this PR |
| --- | ---: | ---: |
| *(default `x86-64` baseline)* | 59.84 µs | 46.76 µs |
| `-C target-feature=+adx,+bmi2` | 53.93 µs (−9.9%) | 41.24 µs (−11.8%) |
| `-C target-cpu=native` | 54.19 µs (−9.4%) | 40.40 µs (−13.6%) |

This is a documentation/deployment finding for consumers, not a code change:
`-C target-feature` cannot be set for a downstream binary from this crate's
`build.rs`, and doing so would be wrong.

---

## 4. What *was* changed: `square` no longer goes through `pow2k(1)`

This came out of the "`pow2k` versus repeated `square`" item in the inventory,
and it is the first of the four code changes in this work. §5 and §6 are the
other two covered by the combined certification in §6.4; §7 came later and is
measured separately, on top of that certified state.

### 4.1 The observation

`pow2k(k)` exists to amortize a squaring's call and loop setup over `k`
squarings. Before this change, `serial::u64`'s `square()` was literally
`self.pow2k(1)`, so it paid for all of that machinery and got none of the
amortization back:

| | before |
| --- | ---: |
| `pow2k(50)`, per squaring | 14.37 ns |
| `square()` = `pow2k(1)` | 22.51 ns |
| overhead per squaring | **+8.1 ns (+57%)** |

The ladder performs **1020** of those `pow2k(1)` calls per `mul_clamped`, and
the inversion — the one place `pow2k(k>1)` is actually used — only 254 squarings
in 12 calls. So `pow2k` *is* used where it should be, and its amortization is
real and still worth having on this target; the defect was that the common case
was routed through it.

Two independent signals said the same thing: `fiat_u64`, which has a dedicated
`square`, measured 15.29 ns against this crate's 22.51 ns, and `serial::u32`
(used by wasm32) already has a dedicated `square()` with `pow2k` written as
repeated squaring — so the defect was specific to `serial::u64`.

### 4.2 The change

`src/backend/serial/u64/field.rs`: the body of `pow2k`'s loop is factored into
`#[inline(always)] fn square_limbs(a: [u64; 5]) -> [u64; 5]`. `pow2k(k)` is now
a loop over it and `square()`/`square2()` call it directly.

The arithmetic is moved verbatim — same formulas, same carry chain, same
`debug_assert!` preconditions, same 5×51 limb representation, same documented
bit excess. No `unsafe`, no new dependency, no public API change, and no change
to `serial::u32`, `fiat_u32` or `fiat_u64`.

No `#[inline]` attribute is placed on `square`. Three variants were measured
across four codegen profiles; `#[inline]` (the hint) produced a **+9.2%
regression** under `lto=off, codegen-units=1`, and `#[inline(always)]` was
statistically indistinguishable from no attribute at all (57.2–58.0 µs for both
across two repeated runs of each of two profiles), so the plain form was kept.

### 4.3 Results

All numbers are ns, minimum of 15 repetitions, harness profile
(`lto=fat, codegen-units=1`), `taskset -c 2`.

**x86_64, `serial::u64` (the default backend's field arithmetic)**

| flags | kernel | before | after | change |
| --- | --- | ---: | ---: | ---: |
| default | `fe_mul` | 26.11 | 25.26 | −3.3% |
| default | `fe_square` | 22.51 | **15.05** | **−33.1%** |
| default | `fe_pow2k50` (per sq.) | 14.37 | 14.09 | −1.9% |
| default | `x25519_mul_clamped` | 59 841 | **46 755** | **−21.9%** |
| default | `x25519_mul_base_clamped` | 18 765 | 18 938 | +0.9% |
| `+adx,+bmi2` | `fe_square` | 23.07 | 14.01 | −39.3% |
| `+adx,+bmi2` | `x25519_mul_clamped` | 53 930 | **41 242** | **−23.5%** |
| `+adx,+bmi2` | `x25519_mul_base_clamped` | 17 352 | 17 427 | +0.4% |
| `target-cpu=native` | `x25519_mul_clamped` | 54 186 | **40 398** | **−25.4%** |

Run-to-run spread (min→max over the 15 repetitions) was 2–11% on the
`mul_clamped` rows above; the `mul_base_clamped` deltas are inside it and should
be read as "unchanged", which is expected — fixed-base multiplication is Edwards
arithmetic and barely touches `square` on the field element type in question.

`mul_base_clamped` was also measured with `curve25519_dalek_backend="serial"`
(18 812 ns) against the default `"simd"` (18 765 ns) to confirm §1's claim that
the vector backend does not reach it.

**Independent confirmation with Criterion** (bench profile forced to
`lto=fat, codegen-units=1`, default RUSTFLAGS, 4 s measurement):

| bench | before | after | change |
| --- | --- | --- | ---: |
| `x25519/MontgomeryPoint::mul_clamped` | 67.93 µs [67.59, 68.32] | **51.73 µs [51.58, 51.89]** | **−23.8%** |
| `field/FieldElement::square` | 15.35 ns [15.29, 15.41] | 14.42 ns [14.33, 14.51] | −6.1% |
| `field/field mix of one mul_clamped` | 48.95 µs | 47.22 µs | −3.5% |

Criterion's intervals here are ±0.3%, and it agrees with the harness on the
headline (−23.8% vs −21.9%). Criterion's absolute `mul_clamped` is higher than
the harness's because of per-iteration harness overhead; the ratio is the
comparable quantity.

**Codegen sensitivity — read this before quoting the headline.** The size of the
win depends strongly on the build profile, because it depends on whether LLVM
chooses to inline the now-loop-free squaring into `differential_add_and_double`:

| profile | before | after | change |
| --- | ---: | ---: | ---: |
| `lto=fat`, `codegen-units=1` | 59 841 | 46 834 | **−21.7%** |
| `lto=thin`, `codegen-units=16` | 59 574 | 57 372 | −3.7% |
| `lto=off`, `codegen-units=1` | 59 049 | 57 432 | −2.7% |
| `lto=off`, `codegen-units=16` (cargo's `--release` default) | 59 447 | 57 380 | −3.5% |

So: **−3.5% for a stock `cargo build --release`, and −22% for a consumer who
enables fat LTO.** No profile regresses. (These are the numbers for this change
alone, measured against the tree immediately before it. For the cumulative
figure across all three changes, certified against the base branch in a single
session, see §6.4.) The pre-change code was flat at
59.0–59.8 µs across all four, i.e. it got nothing from LTO; the refactor is what
makes LTO worth something here. A consumer who cares about X25519 throughput
should set `lto = "fat"` and `-C target-feature=+bmi2`: together those take
`mul_clamped` from 59.8 µs to 41.2 µs on this host, **−31%**.

---

## 5. Second change: the ladder's multiplication by 121666

### 5.1 The observation

`differential_add_and_double` computed `t13` as

```rust
let t13 = &APLUS2_OVER_FOUR * &t6;   // (A + 2)/4 = 121666
```

`APLUS2_OVER_FOUR` is `FieldElement51([121666, 0, 0, 0, 0])` — a 17-bit value in
limb 0 and four zero limbs. Run through the general 5x5 schoolbook that is
**twenty of the twenty-five partial products multiplying by zero**, plus four
`b[i] * 19` precomputations of `0 * 19`.

The compiler does not fix this for you. Disassembling the ladder in a stock
release build (`+adx,+bmi2`):

```text
differential_add_and_double   instrs=774  mulx=60  calls=6
```

The 60 inlined `mulx` are the four squarings (15 each, from §4); all six
multiplications — including the one whose operand is a compile-time constant
with four zero limbs — are real calls into the general `mul`, where the zeros
are invisible. This is exactly what `fe_mul121666` exists for in ref10, donna
and BoringSSL.

### 5.2 The change

A `mul121666` method on each backend's field element, called by the ladder in
place of the general multiplication:

* **`serial::u64`** — the five surviving products and the same carry chain
  `mul` uses. `c[i] = a[i] * 121666 < 2^54 * 2^16.9 = 2^70.9`, so the carries
  are below `2^20` and every bound in `mul`'s commentary holds with room to
  spare.
* **`serial::u32`** — the ten surviving products fed to the existing shared
  `reduce`. Limb 0 of the constant is even-indexed, so none of the
  radix-2^25.5 doubling factors that apply to odd-by-odd products come into
  play.
* **`fiat_u64` / `fiat_u32`** — `fiat_25519_carry_scmul_121666`, which
  fiat-crypto already generates for exactly this constant. The verified
  backends get a verified fast path rather than hand-written arithmetic.

Both hand-written versions are bit-for-bit identical to the general
multiplication by construction, and the differential tests in §12 check it.
`APLUS2_OVER_FOUR` is now `#[cfg(test)]`: the only remaining use is the test
that compares the two against each other.

After the change, the same disassembly:

```text
differential_add_and_double   instrs=824  mulx=65  calls=5  branches=[]
```

One fewer call, and the 25-`mulx` general multiply is replaced by five inlined
`mulx` — precisely the intended transformation. Still no branches.

### 5.3 Results

Isolated, ns per operation, minimum of 15 repetitions:

| target | `mul` | `mul121666` | ratio |
| --- | ---: | ---: | ---: |
| x86_64, `serial::u64` | 25.24 | **9.91** | 2.5x cheaper |
| wasm32, `serial::u32` | 53.74 | **14.04** | 3.8x cheaper |
| wasm32, `fiat_u32` | 53.45 | **16.33** | 3.3x cheaper |

End-to-end `mul_clamped`, ns, minimum of 15 repetitions:

| target / profile | before 5.2 | after 5.2 | change |
| --- | ---: | ---: | ---: |
| x86_64, `lto=off`, `cgu=16` (cargo `--release` default) | 57 380 | **53 935** | **-6.0%** |
| x86_64, `lto=thin`, `cgu=16` | 57 372 | **53 376** | **-7.0%** |
| x86_64, `lto=fat`, `cgu=1` | 46 755 | 46 860 | +0.2% (nothing) |
| x86_64, `lto=fat`, `cgu=1`, `+adx,+bmi2` | 41 242 | **40 857** | -0.9% |
| wasm32, `serial::u32` | 149 607 | **139 100** | **-7.0%** |
| wasm32, `fiat_u32` | 139 752 | **127 942** | **-8.5%** |

**This change and the one in §4 are complementary, and that is the useful part.**
Under fat LTO the inliner already inlines `mul` into the ladder and folds the
zero limbs itself, so §5 buys nothing there — but that is exactly the profile
where §4 is worth 22%. Under a stock release build LTO does neither, so §4 is
worth only 3.5% and §5 picks up another 6%. Neither profile regresses under
either change, and every profile is now meaningfully faster than the base:

wasm32 gains here where it gained nothing from §4, because `serial::u32`
already had a dedicated `square` and never had the `pow2k(1)` defect.

---

## 6. Third change: the ladder's conditional swap

### 6.1 The observation

`mul_bits_be` conditionally swaps its two working points once per scalar bit —
256 times per `mul_clamped`, and the swap must be constant time, so it is a
masked exchange rather than a branch.

Every field backend implements `ConditionallySelectable::conditional_swap`
directly, as a masked exchange of the limbs. `ProjectivePoint`, which is two
field elements, implemented only `conditional_select`, so it inherited
`subtle`'s default:

```rust
fn conditional_swap(a: &mut Self, b: &mut Self, choice: Choice) {
    let t: Self = *a;
    a.conditional_assign(&b, choice);
    b.conditional_assign(&t, choice);
}
```

That is a whole-struct copy plus two conditional assignments — twice the
per-limb work of the exchange the backends already provide, for the one
operation the ladder performs on every single bit. It reads like an oversight
rather than a decision: `FieldElement51` has the override, `ProjectivePoint`
was the type that did not.

### 6.2 The change

Forward `conditional_swap` and `conditional_assign` to the field element's own
implementations. Eight lines, no arithmetic touched, and the constant-time
property is unchanged — a masked exchange is what it was before and what it is
now, only once instead of twice.

### 6.3 Results

The instruction count is the noise-free part of the evidence. The ladder driver
(`<&Scalar as Mul<&MontgomeryPoint>>::mul`, which is where `mul_bits_be` and its
255-iteration loop end up):

| | instructions |
| --- | ---: |
| before | 324 |
| after | **294** |

Thirty fewer per iteration.

End to end, the two targets disagree, and the disagreement is the interesting
part:

| target | before | after | change |
| --- | ---: | ---: | ---: |
| x86_64, `lto=fat` | 46 249 | 45 987 | −0.6%, **not resolvable** |
| wasm32, `serial::u32` | 138 228 | **135 217** | **−2.2%**, 5/5 paired runs |

On x86_64 the difference is below the measurement floor. That floor was
established directly rather than guessed: two *identical* binaries measured
alternately four times varied by 2.5% in their minima, so a 0.6% difference on
this host means nothing. Thirty instructions removed from a loop whose critical
path is a chain of dependent field multiplications is work a wide
out-of-order core absorbs.

On wasm32 it shows up: five alternating paired runs, `with` faster in every one,
best-of-five −2.2%. `mul_base_clamped` is unchanged there (51 836 → 52 166,
52 124 → 51 151), which is the expected control — it does not use the ladder.

Kept on that basis: a measured gain on one target this fork's consumer ships to,
no regression on the other, and it brings `ProjectivePoint` in line with what
every field backend already does.

---

### 6.4 End-to-end certification of changes §4, §5 and §6

**Scope: this certifies §4 + §5 + §6 against `origin/main`. It does not include
§7**, which was measured later and separately, against the tree this section
leaves behind — so the two results compose rather than overlap.

The per-change numbers in §4, §5 and §6 were each measured against the state of
the tree immediately before that change, across several sessions. That is
exactly the kind of bookkeeping that accumulates errors — and this document has
already had to retract one figure (§13.5) — so the cumulative claim is certified
directly instead of being assembled: `origin/main` checked out into a second
worktree, the same harness copied into it, both trees measured **in one session,
alternating, on the same pinned core**. Only the end-to-end kernels are used, so
the base tree builds unmodified and no benchmark hook is involved.

`x25519_mul_clamped`, ns, minimum of 15 repetitions:

| target / profile | `origin/main` | this branch | change |
| --- | ---: | ---: | ---: |
| x86_64, `lto=off`, `cgu=16` (cargo `--release` default) | 58 778 | **53 071** | **−9.7%** |
| x86_64, `lto=off`, `cgu=16`, `+adx,+bmi2` | 54 625 | **49 182** | **−10.0%** |
| x86_64, `lto=fat`, `cgu=1` | 61 547 | **46 188** | **−25.0%** |
| x86_64, `lto=fat`, `cgu=1`, `+adx,+bmi2` | 54 656 | **40 816** | **−25.3%** |
| wasm32, `serial::u32` | 148 046 | **135 081** | **−8.8%** |

The wasm32 cell is best-of-three alternating runs, with this branch faster in
all three. An earlier certification pass — run before §6 landed, so covering
only §4 and §5 — gave −9.5%, −8.8%, −24.2% and −25.6% for the four x86_64 cells,
agreeing with the above to within about a percentage point.

`x25519_mul_base_clamped` is unchanged in every cell (x86_64 19 599 → 19 422,
18 359 → 18 584, 19 223 → 18 684, 17 572 → 17 677; wasm32 51 525 → 51 535). That
is the expected result — none of the three changes touches the fixed-base path —
and it is also the check that nothing regressed there.

Two things earlier passes got wrong, recorded because they are the reason this
section exists. The base branch under fat LTO measures 61.5–62.9 µs across
several runs, not the 59.8 µs recorded during the §4 work, so the fat-LTO
improvement is a quarter rather than the −21.7% assembled from separate
sessions. And host noise varied a great deal between passes — per-run spreads
from 2% to 130% — which is why the minimum is the reported statistic: it needs
one clean repetition, not a clean run.

---

## 7. Fourth change: squaring doubles its 64-bit inputs, not its 128-bit products

### 7.1 The observation

Squaring is the single largest line in the profile that motivated this work —
`pow2k` at 30.8%, which upstream is where every `square()` lands because
`square()` *was* `pow2k(1)` (§4). So it is worth looking at the squaring kernel
itself and not only at how it is reached.

`square_limbs` exploits the symmetry of the 5×5 product: of the twenty-five
partial products, the five squares `a[i]²` appear once and the other twenty
appear in ten mirror pairs, so ten products are computed once and doubled. The
code doubled them like this:

```rust
let c0: u128 = m(a[0], a[0]) + 2*( m(a[1], a4_19) + m(a[2], a3_19) );
```

That is a doubling of a **128-bit** value, five of them, one per output
coefficient — on x86-64 each is a `shld`/`add` pair rather than a single
instruction, and, more to the point, each sits on the dependency chain between
the multiplies and the carry chain.

Since `2*(x*y) == (2*x)*y` exactly over the integers, the doubling can move onto
one 64-bit operand instead, where four precomputed values cover all ten doubled
products:

```rust
let a0_2 = 2 * a[0];   // covers a0·a1, a0·a2, a0·a3, a0·a4
let a1_2 = 2 * a[1];   // covers a1·a4_19, a1·a2, a1·a3
let a2_2 = 2 * a[2];   // covers a2·a3_19, a2·a4_19
let a3_19_2 = 2 * a3_19; // covers a4·a3_19
```

The tell that this is the right shape is that **`serial::u32` already does it** —
`square_inner` has precomputed `x0_2 … x7_2` and has since the code was written.
The u64 backend was the outlier, carrying a comment saying the two forms "don't
seem any better or worse", which this section is the re-measurement of.

### 7.2 The change

`square_limbs` in `backend/serial/u64/field.rs` only. The coefficients are
identical integers, so the output is bit-for-bit identical by construction, and
the existing differential tests — which compare against a verbatim copy of the
pre-refactor code, including at the upper bit-excess bound — check it rather
than assume it.

The bit-excess preconditions are untouched and remain slacker than what the
function already requires: `2*a[i]` fits a u64 for b < 12 and `2*19*a[3]` for
b < 7.75, against the b < 3 the carry chain demands.

### 7.3 Results

Instructions, exact, `callgrind`, `lto=true`/`cgu=1`, baseline `x86-64` ISA:

| | baseline | this change | Δ |
| --- | ---: | ---: | ---: |
| `fe_square` ×20 000 | 3 024 484 | 3 002 201 | −0.74% |
| `mul_clamped` ×20 | 10 449 042 | 10 314 701 | **−1.29%** |
| `mul_base_clamped` ×20 | 4 435 161 | 4 405 500 | −0.67% |

Time, minimum of 15 repetitions, paired alternating runs of two binaries built
from the same tree with only this file differing:

| | baseline `x86-64` | | `+avx2,+adx,+bmi2` | |
| --- | ---: | ---: | ---: | ---: |
| | Δ | paired wins | Δ | paired wins |
| `fe_square` | **−5.4%** | 7/7 | **−10.5%** | 3/3 |
| `fe_invert` | — | — | **−5.8%** | 3/3 |
| `mul_clamped` | −0.7% | 6/7 | **−2.2%** | 3/3 |
| `mul_base_clamped` | −2.8% | — | −1.8% | 3/3 |

**The time win is several times the instruction win, which is the interesting
part.** −0.74% of instructions produced −5.4% of time on the isolated squaring,
and −10.5% once `+bmi2` makes the multiplies cheap enough for the shift to
matter relatively more. That is a latency effect: the 128-bit doubling was on
the critical path between the multiplies and the carry chain, and `callgrind`
counts instructions, not dependency chains. It is the mirror image of the §13.7
caution — there, removing instructions cost time; here, barely removing any
instructions saved a good deal of it. Neither instrument is sufficient alone.

Applying the brief's own rule — *the number that counts is `mul_clamped`'s* —
this lands at **−2.2%** on the configuration this crate recommends and the
profiled machine has (Zen 4 has BMI2), and at **−0.7%** on a stock baseline-ISA
build. Both are at or below this host's timing noise floor in magnitude, so the
case rests on the pairing (6/7 and 3/3, and 7/7 on the isolated kernel) plus an
exact instruction count that moves the same direction.

It was accepted rather than refused because it is the rare change with no other
side: the output is bit-for-bit identical, there is no `unsafe`, no
representation change, no API change, and it makes the u64 backend agree with
what the u32 backend has always done. It also speeds up every other caller of
squaring in the crate — Ed25519 verification, `invert`, `sqrt_ratio_i` — not
just X25519.

§13.7 warns that the ladder is register-pressure bound and that changes adding
live values lose. This one adds four live `u64`, so that risk was real; the
instruction count going *down* rather than up is the evidence it did not
materialize.

**wasm32 cannot be affected, and this was verified rather than argued:**
`serial::u64` is not compiled for `wasm32-unknown-unknown` (which takes
`bits="32"`, §10). Building the harness module from the same directory with only
`u64/field.rs` toggled between the two versions produces a byte-identical
`.wasm` (md5 `5364e78e…` both ways).

---

## 8. Fifth change: the ladder's subtractions skip their reduction

### 8.1 The observation

§13.10 established that the ladder is sensitive to instruction *count* and
indifferent to latency, because the core already overlaps its independent
operation pairs. That makes "remove instructions without adding live values" the
only shape worth looking for — and points at the four subtractions in the ladder
step, which are more expensive than they look.

`Sub` must accept anything satisfying the crate-wide bit excess `b < 3`, so it
adds `16p` to keep the difference positive. That leaves limbs around `2^55`,
outside the contract, so it has to finish with a `reduce`: five shifts, five
masks, five adds and the `×19` fold, about sixteen instructions on top of the
ten the subtraction itself needs. The source has always known this:

> If we could statically track the bitlengths of the limbs of every
> `FieldElement51`, we could choose a multiple of `p` just bigger than `_rhs`
> and avoid having to do a reduction.

In the ladder we *can* track them, by hand and locally. All four subtractions —
`t1 = U_P − W_P`, `t3 = U_Q − W_Q`, `t6 = t4 − t5`, `t10 = t7 − t8` — take
operands that are outputs of `mul` or `square`, and those are bounded by
`2^51 + 2^13`. Against operands that narrow, `2p` is a large enough offset, and
`2p` leaves the result below `2^53` — already inside `b < 3`, with no reduction
needed at all.

### 8.2 The change

A separate `FieldElement51::sub_unreduced`, **not** a change to `Sub`. This
matters: the brief protects the documented bit-excess preconditions as
correctness invariants, and narrowing `Sub`'s would break that contract for
every other caller in the crate. Adding an operation that states its own
stricter precondition does not, and it is the same shape as `mul121666` in §5.

```rust
pub(crate) fn sub_unreduced(&self, rhs: &FieldElement51) -> FieldElement51
```

- **Precondition:** every limb of both operands `< 2^52 − 38`, the smallest limb
  of `2p`. `mul`/`square` outputs satisfy it with eleven bits to spare, and
  `debug_assert!`s enforce it.
- **Postcondition:** limbs `< 2^53`, so `b < 3` holds and the result may be fed
  to anything in the module. The margin against the `2^54` the module allows is
  a full factor of two.

`fiat` gets a `sub_unreduced` that forwards to the ordinary `Sub`, so
`montgomery.rs` stays backend-agnostic and its behaviour is untouched.
`serial::u32` takes the offset too, with a much thinner margin and an assertion
holding it — §8.4.

### 8.3 Results

Instructions on `mul_clamped`, exact, over 20 calls: **10 314 701 → 9 636 449,
−6.58%.** That is larger than the arithmetic predicts (four reductions × ~16
instructions × 255 steps ≈ 16 000, or 3%) because dropping the reduction also
shortens live ranges and takes some of the step's spill traffic with it.

Time, minimum of 15, paired alternating runs:

| | `mul_clamped` | wins | `mul_base_clamped` |
| --- | ---: | :---: | ---: |
| baseline `x86-64`, §7 only | 43 035 | | 15 450 |
| baseline `x86-64`, + §8 | **40 428 (−6.1%)** | 3/3 | unchanged |
| `+avx2,+bmi2`, §7 only | 39 105 | | 15 450 |
| `+avx2,+bmi2`, + §8 | **37 641 (−3.7%)** | 3/3 | unchanged |

`mul_base_clamped` does not move, which is the expected result and the check
that nothing regressed: there is no Montgomery ladder in key generation.

This is the largest single-change improvement to `mul_clamped` in this document
after §4, and unlike §4 it does not depend on the build profile.

**Correctness.** The result is congruent to, but not limb-for-limb equal to,
`self - rhs` — it is a different representative. So the bit-for-bit argument
used everywhere else in this document does not apply, and the checks are
different in kind:

- `sub_unreduced_agrees_with_sub` builds operands the way the ladder does — by
  actually multiplying and squaring random elements rather than assuming what
  those outputs look like — and asserts both that `to_bytes()` matches the
  general subtraction and that every limb stays under `2^53`, the documented
  postcondition rather than the looser `2^54` the module allows. 1024 pairs, plus
  the worst case the precondition permits and the degenerate inputs.
- The `debug_assert!`s on the precondition were exercised by the entire test
  suite in a debug build, which includes the ladder. A violation anywhere would
  have panicked.
- `rfc7748_ladder_test2`, the 1000-iteration iterated-ladder vector that is
  `#[ignore]`d for cost, passes in release (45.2 s), as do the other RFC 7748
  vectors. The ladder's output is canonical bytes, so this is an exact check of
  the thing that actually ships.

**Constant time.** `sub_unreduced` is five adds and five subtracts, straight
line, no branch and no data-dependent index — strictly less work than the `Sub`
it replaces, which was itself branch-free.

### 8.4 `serial::u32` takes it too, on 12% of margin held by an assertion

**This section previously said the opposite.** It argued the margin was too thin
and recorded the analysis as a refusal. The refusal was reversed; what follows
is what ships, and §8.5 keeps the original reasoning because the objection was
sound on its own terms and only one thing about it changed.

The ten-limb layout documents `b < 1.75`, bounded by `19 * y[i]` having to fit in
a `u32`: `26 + b + lg(19) < 32`, i.e. `b < 1.752`. Offsetting by `2p` limb-wise
and dropping the reduction leaves each even limb below `2^26.007 + 2^27`, i.e.
`3.005 * 2^26`, and the same multiple of `2^25` for the odd ones — bit excess
**b = 1.5875**.

So it fits by **0.167 bits, about 12%**. The u64 case has a factor of two.

| | `serial::u64` | `serial::u32` |
| --- | ---: | ---: |
| bound that binds | limbs `< 2^54` | `19 * limb` fits a `u32` |
| what the offset produces | `b < 3` | `b = 1.5875` |
| headroom | 2× | **1.12×** |

**What changed is the failure mode, not the margin.** The objection was that a
`u32` wraparound here is silent: not a panic, not a failed bound assertion, but
a wrong field element inside a constant-time primitive. That is still true of an
unchecked implementation. So the bound is now *asserted* rather than argued.
`debug_assert!` checks the binding constraint itself — that `19 * limb` still
fits a `u32` — on every limb of every result, plus the no-underflow
precondition on the operands. A future change that widens an input past what
this absorbs fails loudly in debug and test builds instead of returning a wrong
answer quietly.

`sub_unreduced_agrees_with_sub` covers 1024 pairs of real `mul`/`square` outputs
and then the worst case the precondition admits — `self` at the top of what
`reduce` returns, `rhs` at zero, which maximises every output limb — asserting
agreement with `Sub`, that `19 * peak` still fits, and that the margin has not
moved from ~12%.

**Result.** `x25519_mul_clamped` on the u32 backend, x86_64 host,
`bits="32"`: **1 238 366 → 1 184 050, −4.39%.** That is the backend wasm32 uses,
which is the slowest target this crate has.

> **wasm32 figures elsewhere in this document predate this change.** §13.7's
> `bits="32"` column and §12's wasm ladder numbers were measured when
> `serial::u32::sub_unreduced` still forwarded to `Sub`, so they understate the
> current default there by roughly this 4.39% on the ladder. They have not been
> re-run under Node; the x86_64 `bits="32"` measurement above is the one that
> supports the claim. Flagged rather than silently left, because the earlier
> revision of this section asserted wasm32 was *unchanged*, and it is not.

### 8.5 The refusal this reversed, kept

The argument against was: twelve percent of headroom, on a silent-wraparound
failure mode, in a constant-time primitive, on the target where this document
has the least direct visibility, in exchange for a few percent of wasm32 DH.

That reasoning is still the reason the assertions exist. It was wrong only about
what to do next — the answer to an unverified bound is to verify it, not to
decline the change. It is kept because someone who wants to widen these inputs
later needs to meet this objection, not rediscover it.

---

## 9. The signature-verification path, and the inversions it does not do

Everything above came from a profile of direct messaging, where X25519 runs once
per message. A second profile of the same client — group messaging, Sender Key,
128 members — reports a different shape: `vector::avx2` at 40.8% against
`serial::u64` at 15.4%, where the direct-message profile had ~0% and ~64%. It
also reports a call count that looks like an invitation: eighteen field
inversions per message, sixteen of them from projective-to-affine conversion.

As before, those numbers come from the consumer and are the reason for the work,
not the specification. What follows is measured inside this repository, with a
new `ed25519_verify` kernel in the harness (`ed25519-dalek` is a bench-only path
dependency of the harness, not of `curve25519-dalek`).

### 9.1 What verification actually costs

> **Corrected.** The figures first published in this section were measured
> against the *wrong crate*, and the error and its consequences are recorded in
> §9.8 rather than quietly overwritten. Everything below is the re-measurement.

`callgrind`, AVX2 backend live, `lto=true`/`cgu=1`. Setup is removed exactly by
differencing 220 iterations against 20, so these are per-verify figures with no
amortized key generation or hashing in them:

**330 457 instructions per `ed25519_verify`.**

| | instructions | share |
| --- | ---: | ---: |
| `avx2::FieldElement2625x4::mul` | 170 560 | **51.61%** |
| `avx2::FieldElement2625x4::square_and_negate_D` | 64 250 | 19.44% |
| `vartime_double_base::spec_avx2::mul`, self | 53 938 | 16.32% |
| **`serial::u64::FieldElement51::pow2k` — the one inversion** | **30 724** | **9.30%** |
| `sha2::sha512` | 4 447 | 1.35% |
| `serial::u64::FieldElement51::mul` (the inversion's) | 2 860 | 0.87% |

**The vector backend is 87.4% of a verification** (the three AVX2 rows). That is
the same conclusion as §13.3, reached from the other end: the vector backend
delivers where it is wired in, and signature verification is where it is wired
in. It is also the in-repository counterpart of the group profile's 40.8%.

The runtime `NafLookupTable5` build no longer appears as its own row. With
`precomputed-tables` on — which is the default, and which the corrected harness
now actually links — the basepoint's table is the static `NafLookupTable8`
constant, so only **one** table is built at runtime (for `A`), and it inlines
into the caller. That is why the `vartime_double_base` self row grows from
15.08% to 17.09% while the total falls.

### 9.2 The sixteen inversions are not in verification. There is exactly one

This is the question the batch hinges on, and the answer is available from the
source before any measurement.

`as_affine` exists once in this crate — `montgomery::ProjectivePoint::as_affine`,
`montgomery.rs:427`. It is called from exactly one place, the last line of
`MontgomeryPoint::mul_bits_be`, so it runs **once per X25519 scalar
multiplication**. `curve_models::ProjectivePoint`, the Edwards projective type
the profile's symbol names, **has no `as_affine` at all**.

The verification path does not call it. Verification is

```rust
EdwardsPoint::vartime_double_scalar_mul_basepoint(&k, &minus_A, &s).compress()
```

and the inversions on that path are:

| site | inversions |
| --- | ---: |
| `NafLookupTable5<CachedPoint>::from(A)` — the vector table build | **0** |
| the wNAF loop — doublings and additions in extended coordinates | **0** |
| `EdwardsPoint::compress` — one `Z.invert()` | **1** |

Measured, not assumed: `callgrind` counts **12.0 `pow2k` calls per `compress`**,
which is one Fermat inversion (`pow22501`'s eleven plus `invert`'s own), and
**zero** `pow2k` anywhere else in the verify subtree. Grepping agrees from the
other direction: every `invert` in `backend/vector/` is inside a `#[cfg(test)]`
module. **The vector backend performs no field inversion at runtime at all.**

So: **the answer to "are the sixteen batchable" is that there are not sixteen.**
There is one, and one inversion cannot be batched with itself. The batch does not
apply to verification, and no amount of restructuring inside this crate changes
that.

**A cross-check on the premise, offered because it may be useful to whoever
owns the consumer.** If the sixteen were sixteen `mul_bits_be` calls, they would
cost sixteen X25519 scalar multiplications. One measures 32 µs here, about
90 000 cycles; sixteen is roughly 1.4 M cycles, which is four times the 358 906
cycles the same profile reports for the *entire* message. So the sixteen cannot
be sixteen X25519 operations either. The only construct in this crate that emits
eight inversions from a single call is `LookupTable<AffineNielsPoint>::from` —
via `as_affine_niels`, a name that truncates to something very close to
`as_affine` — and two of those is sixteen. That is §9.3, and it is where the
batch does apply.

### 9.3 Sixth change: basepoint-table construction batches its inversions

There are exactly three inverting conversions in the crate, all in `window.rs`,
all of the same shape:

| construction | inversions | callers |
| --- | ---: | --- |
| `LookupTable<AffineNielsPoint>::from` | 8 | `BasepointTable::create`, ×32 |
| `NafLookupTable5<AffineNielsPoint>::from` | 8 | **none** |
| `NafLookupTable8<AffineNielsPoint>::from` | 64 | **`serial::scalar_mul::precomputed_straus.rs:44`** — see the correction below |

Each is written as a chain — convert, add, convert, add — and in that form the
conversions really are sequential: `points[j+1]` is computed from `points[j]`,
which is already affine. Taken literally there is nothing to batch, which is the
trap. The dependency is on the *value* of the previous multiple, not on its
affine form, so the chain can run in extended coordinates, where no inversion is
needed, with all the conversions done together at the end.

That makes `EdwardsBasepointTable::create` the one place where Montgomery's
trick applies, and it is not a marginal one: it builds 32 tables of 8, so
**256 field inversions**, and at 3 685 ns each that is 943 µs of the 1 034 µs it
takes — **91% of building a table is field inversion.** The existing source
comment, `XXX batch inversion would be good if perf mattered here`, is now
answered with a number.

The change keeps the additions but defers the conversions, then inverts all
`Z` at once. The batch is a fixed-size array rather than a `Vec`, so this stays
available without `alloc`; `FieldElement::internal_invert_batch` became
`pub(crate)` for that, which is a visibility change and not an API change.

**Results.** Paired alternating runs of two binaries from the same tree with only
`window.rs` differing — seven alternating rounds, minimum of 15 repetitions each:

| | sequential | batched | change | paired wins |
| --- | ---: | ---: | ---: | :---: |
| `EdwardsBasepointTable::create` | 1 033 829 ns | **260 676 ns** | **−74.8% (3.97×)** | 7/7 |
| instructions per `create` (callgrind, exact) | 9 616 631 | **2 678 042** | **−72.1%** | — |
| `x25519_mul_base_clamped` | 15 555 | 15 755 | unchanged | — |
| ~~`ed25519_verify`~~ | ~~38 919~~ | ~~37 936~~ | **withdrawn, see §9.7** | — |

The unchanged `mul_base_clamped` row is the point of including it: the crate
ships its basepoint table as a constant, so nothing on either hot path builds
one. This helps a consumer that calls `BasepointTable::create` at runtime for a
point of its own, and it helps nobody else.

The `ed25519_verify` row is **withdrawn**. It was measured against a different
copy of this crate (§9.7), so it was guaranteed to report "unchanged" whatever
this change did — it was not evidence. The claim it was offered for is better
made without it, and was always available: verification never calls
`BasepointTable::create`, which is a fact about the call graph rather than a
measurement, so batching that constructor cannot move it in either direction.

> **Correction.** The figures above were re-measured. The `edwards_table_create`
> kernel called `EdwardsPoint::mul_base` and `compress` *inside* the timed loop,
> so every repetition charged `create` with a fixed-base scalar multiplication
> and an inversion that `create` does not perform. The kernel now builds its
> input point once outside the loop — `black_box` on that input is what keeps
> `create` from being hoisted as loop-invariant — and compresses once at the
> end, outside the timing boundary.
>
> The overhead was **constant**, so it entered numerator and denominator alike
> and pulled the ratio *toward 1*: the published 4.19× was an **under**statement,
> not an overstatement. Measured on this tree, removing it takes 15.8 µs off
> each arm and moves the ratio from 3.80× to 3.97×. The callgrind counts settle
> this independently, being deterministic and host-independent: they fall by
> 268 370 and 264 652 instructions, the same constant to within 1.4%, which is
> the `mul_base` and `compress` the kernel was charging to `create`.
>
> The published timings also came from a different host and an older tree, so
> only the same-tree pairing above should be read as the effect of the change.
> The direction, the conclusion, and the 91% inversion share are unaffected —
> the last one is now better supported, at 91.2% against the corrected
> denominator. Thanks to Codex for catching it.

**Correctness.** Montgomery's trick returns a different *representative* of each
inverse than `invert` does — both are weakly reduced, and weak reduction is not
unique — so unlike every other change in this document this one is **not**
limb-for-limb identical, and claiming otherwise would be wrong. The tests compare
canonical bytes instead, which is the strongest true statement and the one the
table's users depend on:

* `batched_lookup_table_matches_reference` — 32 distinct base points, all eight
  entries, all three coordinates of each, against a verbatim copy of the
  pre-batch construction.
* `batched_lookup_table_handles_identity` — the identity is the input whose `Z`
  handling differs most between the two (`internal_invert_batch` special-cases
  zeros), and a table built from it must still select the identity everywhere.
* `created_table_agrees_with_precomputed` — a table built the new way multiplies
  identically to the crate's own precomputed constant, over 16 scalars.

**Constant time.** `internal_invert_batch` was already the crate's constant-time
batch inversion, used by `Scalar::invert_batch`: it is straight-line, and skips
zero inputs with `conditional_assign` rather than a branch. Table construction
takes a public point, so this is not a secret-dependent path in the first place,
but the change does not introduce a branch either way.

### 9.4 What a batch is worth, priced

Since the batch has exactly one applicable site, it is worth recording what it
would have been worth elsewhere. Per element, minimum of 5 repetitions:

| n | n separate `invert` | one `invert_batch` | ratio |
| ---: | ---: | ---: | ---: |
| 1 | 3 528 ns | — | — |
| 8 | 3 534 ns/element | **578 ns/element** | **6.1×** |
| 16 | 3 535 ns/element | **355 ns/element** | **9.9×** |

The same comparison on wasm32 (32-bit limbs, serial backend) gives **5.6×** at
n = 8 (9 138 → 1 634 ns/element) and **8.5×** at n = 16 (9 062 → 1 065). The
ratio survives a target with no 64x64→128 multiply and no vector backend, which
is what makes it a property of the algorithm rather than of one machine. These
numbers only exist because the wasm driver&#39;s kernel list stopped being a
hand-maintained copy — see §14.

Montgomery's trick replaces \\(n\\) inversions with one inversion and
\\(3(n-1)\\) multiplications, so the ratio grows with \\(n\\) and is bounded by
the ratio of an inversion to three multiplications. These kernels ship even
though only one site could use them, because they are what makes the refusals in
§9.2 and §9.5 checkable rather than assertions.

### 9.5 safegcd, reassessed against verification. Still refused

§13.11 deferred safegcd with the inversion at 5.3% of `mul_clamped`. On the
verification path it is worth more: 9.30% in `pow2k` plus 0.87% in the
inversion's multiplications, so **10.2% of a verification** (§9.1 as corrected).
Against
the group profile's 40.8% for verification that is about 4% of the client.

Refused again, on three counts, and the first is new:

1. **The ceiling is low and the alternative is better where it exists.** A
   safegcd inversion is usually quoted at 2–4× Fermat, so 10.2% becomes perhaps
   5–7% of a verification, i.e. 2–3% of the group client. Meanwhile §9.4 shows
   the crate's *existing* batch inversion is worth 6–10× wherever several
   inversions co-occur. The algorithmic lever with real leverage is "find a set
   of inversions to batch", not "make one inversion faster" — and §9.2 shows
   verification has exactly one.
2. **It cannot be a variable-time inversion, even though verification is
   variable time.** `compress` is shared: it is called on secret-derived points
   elsewhere in the workspace, including in signing. A vartime inversion reached
   from `compress` would be a timing leak in those callers, and splitting
   `compress` into two variants is public API surface this fork should not add
   while it is trying to merge upstream.
3. **Unchanged from before:** a constant-time safegcd's constant-time property is
   global to the iteration count rather than local to each step, so it cannot be
   checked the way every change in this document has been checked, and
   fiat-crypto does not generate `divstep` for curve25519, so there is no
   verified primitive to build on.

### 9.6 Headroom in `vartime_double_scalar_mul_basepoint`: none in the tables

The obvious follow-up is whether the wNAF window or the table has slack, and the
first thing to check is whether §13.5's rejection of bigger tables transfers.
**It does not, and for a reason worth stating:** §13.5 rejected larger tables
because `LookupTable::select` is a *constant-time linear scan*, so scan cost
grows with table size. The vartime path's `NafLookupTable5::select` is

```rust
self.0[x / 2]
```

a direct index. Table size costs nothing to read here. The scan argument is
simply absent, and the trade is a different one: build cost against addition
count.

That trade is already at its optimum. For a 256-bit scalar, width-\\(w\\) wNAF
averages \\(256/(w+1)\\) additions and needs \\(2^{w-2}\\) table entries, which
cost \\(2^{w-2}-1\\) additions to build:

| \\(w\\) | entries | build additions | loop additions | total |
| ---: | ---: | ---: | ---: | ---: |
| 4 | 4 | 3 | 51.2 | 54.2 |
| **5 (what the crate uses)** | **8** | **7** | **42.7** | **49.7** |
| 6 | 16 | 15 | 36.6 | 51.6 |
| 7 | 32 | 31 | 32.0 | 63.0 |

Width 5 is the minimum, and the crate already uses it. Going to width 6 would
double the table build to gain about six additions out of a loop that also
performs 256 doublings. The conclusion rests on that addition count, which is
arithmetic and does not depend on any measurement; it was not attempted beyond
the arithmetic because the arithmetic is not close.

**One supporting number here was wrong** and is withdrawn: this section
originally said "the table build is 17 778 instructions per verification, 5.0%
of it". That was measured against a build *without* `precomputed-tables`
(§9.7), where the basepoint's table is also built at runtime — so it was two
builds, not one, and priced against an inflated total. On the corrected build
the single remaining build inlines into its caller and does not separate out at
all. The ranking of widths is unaffected.

Consistently, the crate uses width **8** for the *precomputed* basepoint table,
where the build cost is paid once at compile time and the only term left is the
addition count. That is the same trade resolved with one term deleted, and it is
resolved the other way, which is a good sign the model is right.

### 9.7 The verification kernel was measuring the wrong crate

This is recorded rather than quietly fixed, because the failure mode is worth
more than the numbers it cost.

`ed25519-dalek` depends on `curve25519-dalek` **by version**, not by path. The
benchmark harness is its own workspace, and its `[patch.crates-io]` covered
`curve25519-dalek-derive` but not `curve25519-dalek` itself. So Cargo resolved
`ed25519-dalek`'s dependency to the *published* 5.0.0 from crates.io, and the
harness binary linked **two copies** of the crate: the path copy under the field
and ladder kernels, the registry copy under `ed25519_verify`. The lockfile said
so plainly, with two `name = "curve25519-dalek"` entries, one carrying a
`source = "registry+…"` and a checksum.

**The consequence is worse than wrong numbers.** `ed25519_verify` was blind to
this repository. It would have reported "unchanged" for a real regression
exactly as readily as for a real improvement — which is the worst failure a
measurement tool can have, because it fails silently and in the direction of
reassurance. Every use of that kernel as a *control* was vacuous, including the
§9.3 row that offered it as evidence that batched table construction leaves the
hot paths alone.

What it cost, all now corrected in place:

| claim | published | corrected |
| --- | ---: | ---: |
| instructions per verification (§9.1) | 354 575 | **330 457** |
| AVX2 share of a verification (§9.1) | 87.8% | **87.4%** |
| the one inversion's share (§9.5) | 10.0% | **10.2%** |
| runtime NAF table build (§9.6) | 17 778 Ir, 5.0% | one build, inlined |

**The first attempt at this correction was itself wrong, and that is worth
recording too.** It published 315 697 rather than 330 457. The re-measurement
had been taken while a concurrent experiment sat uncommitted in the working
tree, so it was measuring a tree that was neither the old code nor the new. The
tell was visible and missed: 315 697 is within 0.3% of what the same kernel
reports *with* that experiment applied. The figure above was re-taken in a
pristine `git worktree` checked out at the commit in question, with nothing else
running, and reproduces exactly across two independent runs. **Check
`git status` immediately before a measurement, not merely at the end of the
session** — the same discipline §12 already applies to `cfg`s and features.

The registry build was 7.3% more expensive for two structural reasons, both
downstream of `precomputed-tables` being off in it: the basepoint's odd-multiples
table was built at runtime instead of being the static constant, so there were
two table builds per verification rather than one; and the basepoint used a
width-5 wNAF instead of width 8, costing roughly fourteen extra additions.

**No conclusion in §9 changes.** The load-bearing arguments there are facts
about the call graph — that `as_affine` is called from exactly one place, that
the vector backend contains no runtime `invert`, that `compress` performs one
inversion — and those were established by reading the source and by counting
`pow2k` calls, both of which are true of either copy. The safegcd refusal moves
from 10.0% to 10.6%, which strengthens it as an argument and does not come close
to reversing it.

The fix is one line of `[patch.crates-io]` in the harness manifest, with a
comment saying why it must stay. The general lesson is the one this document
already applies to `cfg`s and features in §12: **a measurement you have not
verified is measuring the thing you changed is not a measurement.** The cheap
check that would have caught it is to make a change you *know* is large, and
confirm the number moves.

---

### 9.8 Where a group message's `serial::u64` time goes

The group profile's 15.4% in `serial::u64`, against 87.4% of a *verification*
being AVX2, is consistent: whatever serial time a group message spends is not
inside verification, except for the single `compress` inversion. From this
crate's side there are exactly two sources it can be:

1. **Montgomery ladder work** — every X25519 operation is `serial::u64` in every
   backend (§13.3), and this is precisely the work §4 through §8 reduced by
   10.1% stock and 28.7% under fat LTO.
2. **`compress`'s inversion**, at 10.2% of each verification.

There is no third: `mul_base` is serial in every backend too, but it appears in
the group profile only through key generation, not per message. Nothing else in
the crate runs `serial::u64` on that path.

---

## 10. Three changes inside the AVX2 backend

§9 established that a verification is 87.4% AVX2 code, then spent its length on
things that were *not* the vector backend. That was backwards, and it took five
independent explorations of this repository — deliberately run without sight of
this document until they had formed their own candidate lists — to say so. The
first two changes below are in the two functions §9.1's table puts at the top;
§10.5, added later, is in the multiply both of them call.

### 10.1 Seventh change: `ExtendedPoint::double` fuses its two signed blends

`avx2/edwards.rs` accumulated `+S2` into lanes A,D and `-S2` into lanes B,C in
two separate passes:

```rust
tmp0 = tmp0 + zero.blend(S_2, Lanes::AD);
tmp0 = tmp0 + zero.blend(S_2.negate_lazy(), Lanes::BC);
```

The lane sets are disjoint and the quantity is the same, so one signed blend
carries both — `(S2, 2p-S2, 2p-S2, S2)` — for one blend and one add per limb
instead of two of each, 9 ops/limb down to 7.

> **This section claimed the wrong provenance.** It said "the `ifma` backend
> already did it this way; AVX2 was the odd one out." It does not.
> `ifma/edwards.rs` fuses a *different* pair — `-S2` with `-S4`, via
> `S2_S2_S2_S4` — in two blend-and-adds, and is still at **9 ops/limb**. The
> fusion taken here, `+S2` with `-S2` in one signed blend, is not present
> there. The change is right; the claim that someone else had already made it
> was not, and the practical cost of the error is that a live back-port to
> `ifma` went unnoticed for want of anyone checking. Comparing the two vector backends
against each other, rather than each against its own history, is what surfaced
it.

**Limb-identical**, not merely congruent: the two vectors are equal, and the
`2p` in lanes B,C is already present today via `negate_lazy`. Bounds unchanged
at `b < (1.01, 1.6, 2.33, 1.6)` against `mul`'s `(2.5, 1.75)`; `negate_lazy`'s
`b < 0.999` precondition still holds because `S2` leaves `square_and_negate_D`
at `b < 0.007`.

### 10.2 Eighth change: the field multiply is emitted column-major

`FieldElement2625x4::mul` is a 10x10 schoolbook written **row-major**: row `z0`
needs `y1_19..y9_19` while row `z9` needs `y0..y9`, so all ten `y_j`, all nine
`y_j_19`, all ten `x_i` and all five `x_odd_2` are live at once — about **44
values against 16 YMM registers**. The shipped function was 416 instructions, of
which roughly a quarter was spill traffic: 30 stack stores and 64 stack loads.

Emitting it **column-major over `rhs`** instead — column `j` contributing one
product to each of the ten accumulators — lets `y_j` and `y_j_19` be born and die
inside their column, and unpacks `rhs.0[j/2]` only when its column needs it.
Peak liveness drops to about 27.

> **The mechanism stated here was wrong, and the result is unaffected.** This
> section said the `x_i`, invariant across all ten columns, "fold into
> `vpmuludq mem, ymm, ymm` rather than being reloaded." Measured on the shipped
> binary, **36 of 417 `vpmuludq` take a memory operand — 8.6%.** The `x_i` are
> largely *not* folding. What does fold is the accumulator reload, into
> `vpaddq mem, ymm, ymm`: **122 of 508, 24%**.
>
> The −4.23% and the 416 → 388 instructions per call reproduce and are not in
> question; only the explanation of *why* was. It matters because "peak
> liveness drops to 27" reads as though the spilling problem were solved, and
> 27 against 16 registers is still eleven over. The successor question — why
> the `x_i` do not fold — is open, and the answer is not visible from the
> source: `vpunpckldq`, which builds them, requires its data operand in a
> register, so the unpacked form cannot be rematerialised from memory.

The mapping is mechanical: in column `j`, accumulator `z_k` takes `x_i` with
`i = (k - j) mod 10`; the multiplier is `y_j` when `k >= j` and `y_j_19` when
`k < j`; and `x_i` is doubled iff `i` and `j` are both odd. That rule was
**checked against the shipped code before anything was changed** — it reproduces
all 100 of the existing partial products exactly, term for term, which is what
makes this a reordering rather than a re-derivation.

**Bit-identical for every input.** These are the same 100 products in a
different order, and `u64` addition wraps, so it is associative and commutative;
the equality holds even outside the documented bounds. `mul_matches_row_major`
checks it against a verbatim copy of the old row-major code over 512 random
input pairs at the top of the documented ranges, comparing all five limbs across
all eight lanes — limb equality, which is strictly stronger than field equality.

### 10.3 What the two are worth

Paired A/B, same tree, `callgrind`, setup differenced out by the usual
220-against-20 differencing. `ed25519_verify` is the kernel §9.7 had to repair
before it could measure anything at all:

| | `ed25519_verify` | vs base |
| --- | ---: | ---: |
| base (§9.1) | 330 457 | — |
| + §10.1 blend fusion | 327 957 | **-0.76%** |
| + §10.2 column-major | **316 477** | **-4.23%** |

The blend fusion's -2 500 instructions per verification is exactly the
prediction from the disassembly: 5 `vpaddd` and 5 `vpblendd` removed, times 250
doublings. On `edwards_vartime_double_base` in isolation the two are -0.82% and
-4.64% cumulative, and the field multiply itself drops **416 -> 388 Ir per call,
-6.7%**.

Both move every AVX2 point operation — `variable_base`, `straus`,
`precomputed_straus`, `pippenger` and the NAF table builds — not just
verification. Neither touches wasm32, which has no vector backend, so §13.7's
warning that x86-64 instruction counts are not a proxy for wasm32 cannot bite
here.

### 10.4 Why the earlier register-pressure result did not predict this

§13.7 recorded three failed attempts at instruction-level parallelism in the
*serial* ladder and concluded it was register-pressure bound. Every one of those
attempts **added** live values. §10.2 removes seventeen: it is the same
diagnosis pointed the other way, which §13.7 itself half-anticipated — "the one
change in this area that did work went the other way."

The reason it was not found earlier is duller and worth naming: **§12 never
opened `avx2/field.rs`.** The measurement said 87.4% of a verification sat in
two vector functions, and the investigation went to the serial backend anyway,
because that is where the previous front had been. A profile is only useful if
you follow it to where it points.

---

### 10.5 Ninth change: the multiply is tagged per call site

`FieldElement2625x4`'s multiply is called three times on the verification path —
once in `ExtendedPoint::double`, twice in `Add<&CachedPoint>` — and LLVM emits
**one shared outlined body** for all three. Shared, it cannot specialise on the
operand shapes each site actually has.

`mul_tagged::<N>` takes an unused `const N: u8`, which is enough to give each
site its own instantiation. The `Mul` operator remains as the untagged spelling,
delegating to `mul_tagged::<0>`, so nothing outside the three hot sites changes.

| | instructions per `ed25519_verify` | `.text` |
| --- | ---: | ---: |
| before | 316 280 | — |
| tagged | **306 124** | +784 bytes |
| | **−3.21%** | +0.2% |

Wall-clock does not resolve this: the host's spread on `ed25519_verify` ran from
−16.5% to +18.4% across seven paired rounds, which is five times the effect. The
instruction count is the claim, and it is exact rather than statistical. The
`.text` figure is there because executed-instruction count is precisely the
wrong metric for a change that duplicates code — it would not show I-cache
pressure. +784 bytes is small enough that the question does not arise.

This is a compiler-behaviour workaround, and it is worth being clear about what
that means for its lifetime: `mul_tagged::<0>` through `::<3>` all compute the
same function, so if a future LLVM stops sharing the body the tags become inert
rather than wrong. Nothing depends on them for correctness.

**One half of this was refused.** The change arrived paired with replacing `-x`
by `x.negate_lazy()` in `From<ExtendedPoint> for CachedPoint`. Measured alone it
is worth **−47 instructions, −0.015%** — nothing — and it is not free: `neg`
returns `b < 0.0002` where `negate_lazy` returns `b < 1`, so the comment two
lines below it, "the coefficients of the output are bounded with `b < 0.007`",
would become false, and the resulting `CachedPoint` would sit at exactly the
`b < 1.0` the addition documents assuming. Correct today, on that reading, with
no margin and a stale comment, in exchange for nothing measurable. Dropped.

## 11. The constant-time window scan, and what its speed actually costs

`LookupTable::select` is **19.8% of a keygen** — the largest single item on the
fixed-base path. Two candidate fixes existed, each measured on its own, and the
question this section set out to answer was whether they add up. They do not,
and finding out why turned a 6.5% speedup into a 2.4% one.

### 11.1 Where the cost is

The scan keeps a live accumulator and conditionally assigns into it once per
entry:

```rust
let mut t = T::identity();
for j in $range {
    let c = (xabs as u16).ct_eq(&(j as u16));
    t.conditional_assign(&self.0[j - 1], c);
}
```

For `AffineNielsPoint` that accumulator is fifteen `u64` limbs, live across all
eight entries, each entry a read-modify-write per limb. The disassembly showed
78 stack accesses and — the part that matters — **a `call` to
`subtle::black_box`**. `subtle`'s barrier is `#[inline(never)]` around a
`read_volatile`, so every `Choice` is a real function call with a mandatory
memory round-trip, and every caller-saved register holding the accumulator has
to be spilled around it. Eight crossings per select, 64 selects per operation.

### 11.2 Four variants, measured

| variant | `subtle`'s barrier | keygen | Δ |
| --- | --- | ---: | ---: |
| baseline | intact | 202 172 | — |
| **A** `subtle/core_hint_black_box` | **weakened** | 194 812 | −3.64% |
| **B** hand-rolled comparison mask | **bypassed** | 189 116 | **−6.46%** |
| **B′** OR-accumulate, barrier per entry | intact | 202 940 | +0.38% |
| **B″** OR-accumulate, barrier hoisted | **intact** | 197 308 | **−2.41%** |
| A + B | weakened + bypassed | 189 308 | −6.36% |
| A + B″ | weakened | 194 364 | −3.86% |

**They are not additive.** Summing A and B predicts −10.10%; measured together
they are −6.36%, and A on top of B is *worse* than B alone. Every pairing lands
on the same ceiling because every one of them is attacking the same thing.

### 11.3 What the win actually was

B′ is the control, and it is the important row. It applies the whole algorithmic
change — OR-accumulate into a zeroed accumulator (`acc |= mask & v`, two
operations and an independent chain per limb) instead of read-modify-write into
a live one (`t ^= mask & (t ^ v)`, three operations) — while crossing `subtle`'s
barrier exactly where the original does. It is worth **+0.38%**: nothing, very
slightly negative.

So the algorithmic change is not where the speed was. LLVM was already hoisting
the mask across the fifteen limbs, and the extra entry that OR-accumulation
needs for the identity pays back what the cheaper inner operation saves. **The
entire measured win in A and B came from weakening or bypassing the constant-time
barrier**, and B is faster than A only because bypassing beats cheapening.

That reframes both candidates. Neither was an optimization that happened to
touch `subtle`; each *was* the barrier change, wearing different clothes.

### 11.4 Ninth change: hoisting the barrier crossings — **taken**

B″ keeps `ct_eq` and `Choice` exactly as they are and moves all nine crossings
**before** the accumulator exists:

```rust
let mut masks = [0u64; $size + 1];
for j in 0..=$size {
    let c = (xabs as u16).ct_eq(&(j as u16));
    masks[j] = (c.unwrap_u8() as u64).wrapping_neg();
}
// ... then accumulate, with no barrier in the way
```

Same barrier, same number of crossings, same semantics — but nothing expensive
is live when they happen, so the spills they force are cheap. That is worth
**−2.41%** on keygen with the crate's constant-time posture untouched, and it
scales with table size, because a bigger table interleaves more crossings with
the live accumulator:

| | baseline | with §11.4 | Δ |
| --- | ---: | ---: | ---: |
| `x25519_mul_base_clamped` | 202 172 | 197 308 | **−2.41%** |
| `edwards_mul_base` (radix-16) | 170 544 | 165 680 | −2.85% |
| `edwards_mul_base_radix32` | 170 917 | 159 165 | **−6.88%** |
| `edwards_mul_base_radix64` | 191 324 | 171 888 | **−10.16%** |
| `x25519_mul_clamped` (the ladder) | 463 460 | 463 460 | unchanged |
| `ed25519_verify` | 316 477 | 316 477 | unchanged |

The last two rows are the scope check: the Montgomery ladder uses no lookup
table, and verification's `NafLookupTable5::select` is a direct index, so
neither should move, and neither does.

> **The ladder row is a correct measurement behind an incorrect reason.** "The
> Montgomery ladder uses no lookup table" is true and is not why it did not
> move. This section's subject is *barrier crossings under high liveness*, and
> the ladder has **255 of them per DH** — more than any other kernel in the
> crate — at `montgomery.rs:241`, where `choice.into()` is a `Choice::from`
> around `subtle`'s `#[inline(never)]` volatile read. The profile prices it
> exactly: `subtle::black_box` self cost is 46 080 over ~60 operations, i.e.
> three instructions × 255 calls per `mul_clamped`.
>
> Scoping by *where the crossings came from* rather than by *whether there were
> any* is what let the ladder fall through. Hoisting them there is a real, if
> small, candidate — see §11.7, which also records why it is smaller than this
> section's numbers would suggest.

### 11.5 A and B — **refused**

**A, `subtle/core_hint_black_box`**, swaps the `read_volatile` barrier for
`core::hint::black_box`, which the standard library documents as best-effort and
explicitly permits the compiler to ignore. Worse, it is a *dependency feature*,
and Cargo unifies features across the whole graph: a library enabling it decides
the barrier for **every crate in its consumer's build that uses `subtle`**, not
just for itself. That is not a library's call to make. It remains available to a
consumer who wants it — worth a further −1.45 points on top of §11.4 — in the
same category as `+avx2` and fat LTO in §2: a build setting whose trade the
person assembling the binary should choose.

**B, the hand-rolled comparison mask**, is the faster one at −6.46%, and is
refused for a sharper reason. It computes the mask with plain arithmetic
(`((d | -d) >> 63) - 1`), written branchlessly — and the compiler emitted
`sete`, having recognised the value as a boolean. Nothing worse happened here:
the whole-binary branch census was unchanged and the emitted code was
straight-line. But recognising the boolean is precisely what `subtle`'s barrier
exists to prevent, and a scan that selects on a secret scalar digit is the last
place to trade a guarantee for two percent on the strength of one compiler
version's output.

### 11.6 Bigger tables: the answer depends on the build flags

§12 rejected larger basepoint tables, measuring radix-32 as a wash and radix-64
as clearly worse. §11.4 helps bigger tables most, which appeared to overturn
that for radix-32. **It does not, and the first version of this section said it
did — the claim was made without stating a build configuration, and it reverses
under one.**

| | radix-16 | radix-32 | |
| --- | ---: | ---: | --- |
| stock `--release` | 165 680 | 159 165 | **−3.93%** |
| `-C target-feature=+avx2` | 157 670 | **161 476** | **+2.41%** |

Both rows reproduce exactly, and they disagree. The scan grows linearly in table
size while the additions it saves shrink only logarithmically — radix-32 does 52
windows over 16 entries (832) against radix-16's 64 over 8 (512) — so the
trade-off turns on how expensive a scanned entry is relative to an addition, and
`+avx2` moves exactly that ratio. §2 already measured `+avx2` halving the
constant-time scan; what it also does is change which table size wins.

**Under the flags this document recommends for x86_64, radix-16 remains
optimal**, so §12's rejection stands for any build that follows the advice in
§2. Radix-32 is worth it only for a stock build that does not enable `+avx2`,
and it costs twice the table for under 4%. Radix-64 loses in both.

The lesson is the narrow one: a comparison between two implementations is only a
result *together with the flags it was taken under*, and a one-line claim that
omits them can be exactly backwards for the configuration a reader is actually
using. Found by an independent audit that measured the same comparison under
`+avx2` and got the opposite sign.

---

### 11.7 The same hoist on the *vector* scan — **refused, measured**

§13.12 moved the fixed base onto `LookupTable<CachedPoint>::select`, the
generic scan, which never received §11.4's fix — `select_or` is
`impl LookupTable<AffineNielsPoint>` only. An independent audit found this
three times over, from three different angles, and predicted **−5% to −8%** on
`edwards_mul_base` from porting §11.4 across. The diagnosis was right and the
prediction was wrong.

The scan really is ~21% of the vectorised `edwards_mul_base`, and the
disassembly really does show the five-`ymm` accumulator spilled around each of
the nine `subtle` crossings. Two implementations were written and measured:

| | `edwards_mul_base` | |
| --- | ---: | ---: |
| generic `select` (shipped) | 84 546 | — |
| hoisted, masks kept as `[u32x8; 9]` | 89 602 | **+5.98%** |
| hoisted, masks kept as `[u32; 9]`, broadcast at use | 86 658 | **+2.50%** |

Both are worse. The first is easy to explain — nine `u32x8` masks are 288 bytes
of stack, and their reloads cost more than the accumulator spills being removed.
The second removes that and is still worse, which is the interesting half.

The arithmetic says the hoisted version should win: `or_masked_assign` is two
vector ops per limb against `conditional_assign`'s three, so 9 entries × 5 limbs
× 2 = 90 against 8 × 5 × 3 = 120. What eats the 30-op saving is the bookkeeping
the hoist requires and the generic version does not: nine mask stores, nine
reloads, nine `vpbroadcastd`, and a ninth entry for the identity — which
`select` gets free by initialising its accumulator to it, and an OR-accumulator
cannot.

**This is §11.3 repeating.** There, the control (B′) showed the OR-accumulate
*algorithm* was worth nothing on its own and the entire win was the hoist. Here
the hoist is worth nothing on its own, because the serial case's premise does
not transfer: on the serial path the accumulator was fifteen GPRs the compiler
could not keep, and hoisting freed them; on the vector path five `ymm` out of
sixteen fit comfortably, and the crossings were not what was costing.

Recorded rather than dropped, because the analysis is most of the work and
because three independent readings of the disassembly reached the opposite
conclusion. **The spilling is visible and real; removing it this way costs more
than it saves.** Anyone attacking it again should start from the measurement
that the barrier crossings are not the binding constraint here, and should be
sceptical of any prediction — including this document's — that is not a paired
measurement.

## 12. Front B — wasm32 with 64-bit limbs: **refused**

`build.rs` picks `curve25519_dalek_bits` from `target_pointer_width`, so wasm32
gets `DalekBits::Dalek32` (`serial::u32::FieldElement2625`), carrying the note:

```rust
//TODO(Wasm32): Needs tests + benchmarks to back this up
```

The hypothesis was that because wasm32 has native `i64.mul`/`i64.add`, a `u64`
is not emulated the way it would be on a real 32-bit ARM, so `bits="64"` might
win.

### 12.1 Method

Criterion does not run on wasm32-unknown-unknown. Rather than introduce
`wasm-pack` and a second set of kernels, the **same** harness crate is compiled
as a `cdylib` for `wasm32-unknown-unknown` and driven from Node
(`benches/harness/run.mjs`), which instantiates the module and times
`bench_kernel(which, iters)` with `process.hrtime.bigint()`. The module has no
imports and no `wasm-bindgen` glue, so the only thing inside the timed region
besides the kernel is one cross-boundary call per repetition, amortized over a
calibrated iteration count (~20 ms). Same statistic as native: minimum of 15
repetitions. The wasm and native numbers therefore come from literally the same
Rust source, which is the point.

`node v22.22.2`, `--release` with `lto=true, codegen-units=1, panic=abort`.

### 12.2 Results

ns per operation, minimum of 15 repetitions:

| backend | `fe_mul` | `fe_square` | `fe_pow2k50` | `mul_clamped` | `mul_base_clamped` |
| --- | ---: | ---: | ---: | ---: | ---: |
| **`bits=32` (current default)** | **53.71** | **37.75** | **33.27** | **149 607** | **52 161** |
| `bits=64` | 148.70 | 93.10 | 82.90 | 346 175 | 103 752 |
| ratio 64/32 | 2.77× | 2.47× | 2.49× | **2.31× slower** | 1.99× slower |
| `bits=32`, fiat | 53.84 | 36.67 | 31.50 | 139 752 | 53 920 |
| `bits=64`, fiat | 147.99 | 93.04 | 81.91 | 353 674 | 104 166 |

Run-to-run spread was 2.0% on the `bits=32` `mul_clamped` row and 16.5% on the
`bits=64` one; the 2.3× gap is two orders of magnitude larger than the noise.

### 12.3 Verdict

**Refused. `curve25519_dalek_bits="64"` is 2.31× slower than the default on
wasm32**, and the current `build.rs` behaviour is correct.

The reason is that the relevant primitive is not `i64` arithmetic but the
**64×64→128 widening multiply**. `serial::u64` builds every partial product as
`(x as u128) * (y as u128)` and accumulates in `u128`; wasm has no 128-bit
integer type and no widening 64-bit multiply, so each `m(x, y)` becomes a
software 128-bit multiply out of four 32×32→64 pieces plus multi-word
accumulation. `serial::u32` needs only 32×32→64, which is a single `i64.mul`
after two extends. The wasm32 default is a good fit and the 64-bit backend is
not.

This is also a concrete instance of the warning that wasm32 profiles do not
generalize from native: the *sign* is inverted. On x86_64 the 64-bit backend is
the only sensible choice; on wasm32 it costs 2.3×.

The `TODO(Wasm32)` in `build.rs` is updated to point at this document. No
behaviour change: the code path it documents is the one that was already taken.

---

## 13. The rest of the inventory

### 13.1 `pow2k` versus repeated `square`

Answered in §4. Summary: `pow2k` is used where it should be (only the inversion
tail needs `k > 1`), its amortization is real on this target — after the change
`pow2k(50)` is 14.09 ns/squaring against 15.05 ns for a standalone `square`, a
6% edge that justifies keeping it — and the bug was that the *common* case was
routed through it and paid 8.1 ns per squaring for the privilege.

Criterion's `pow2k(k)` ladder after the change (default flags, no forced LTO):
`pow2k(1)` 17.09 ns, `pow2k(2)` 30.64, `pow2k(4)` 60.26, `pow2k(10)` 149.33,
`pow2k(50)` 739.45 (14.79/sq), `pow2k(100)` 1478 (14.78/sq) — i.e. the
per-squaring cost flattens by about k = 10.

### 13.2 Cost of the `fiat` backend

The formally verified backend was measured as the cheap control on both targets.
`mul_clamped`, ns, minimum of 15:

| target | dalek before | dalek after | fiat | fiat `+adx,+bmi2` |
| --- | ---: | ---: | ---: | ---: |
| x86_64, default flags | 59 841 | **46 755** | 53 192 | 51 439 |
| x86_64, `+adx,+bmi2` | 53 930 | **41 242** | — | 51 492 |
| wasm32 | 149 607 (u32, unaffected) | 149 607 | **139 752** | n/a |

Two things worth recording for a consumer:

* **Before this change, `fiat` was the faster backend on x86_64** — 53.19 µs
  against 53.93 µs even with BMI2, and 53.19 against 59.84 without. The cheap
  control was beating the default. After the change the default backend is ahead
  by 11–20%, so the change had to clear `fiat`, not just the old default, and it
  does.
* **On wasm32 `fiat` is still ~7% faster** (139.8 µs vs 149.6 µs) and the
  difference survives the noise. For a wasm consumer,
  `--cfg curve25519_dalek_backend="fiat"` is a real, free improvement on top of
  formally verified arithmetic. `fiat` is *slower* for `mul_base_clamped` on
  x86_64 (21.9 µs vs 18.9 µs), because selecting it also disables the vector
  backend for Edwards work; on wasm32 there is no such trade-off (53.9 vs 52.2 µs,
  inside noise).

`fiat`'s profile is different in shape, not just in magnitude: its `square` is a
dedicated routine (15.29 ns, competitive) but its `pow2k` is repeated squaring
with no amortization (16.65 ns/sq, worse than this crate's 14.09).

### 13.3 Could the vector backend cover the Montgomery ladder?

**Viable, but not worth it for a single X25519, and it is a new backend rather
than an extension of the existing one.**

The AVX2 backend's `FieldElement2625x4` packs four field elements into 32-bit
lanes to evaluate the *four independent coordinate products* of one Edwards
point operation in parallel. `differential_add_and_double` does expose 4-way
width in places — `t4`/`t5` are two independent squarings, `t7`/`t8` two
independent multiplications — but the step as a whole is a serial chain of six
dependent stages, so a lane-parallel ladder would idle most lanes most of the
time and pay pack/unpack on every stage. The published wins for vectorized
X25519 come from batching **independent** scalar multiplications (4 or 8 key
exchanges at once), which is a different API than `mul_clamped` and not what a
request-response protocol like the motivating one has available.

Cost estimate if someone wanted it anyway: a `MontgomeryPoint` ladder over
`FieldElement2625x4` plus a 4-way batched entry point, i.e. new formulas, new
constant-time swap over packed lanes, its own test vectors, and a public API
addition for the batched form — comparable in size to `docs/parallel-formulas.md`'s
Edwards work, and it would not speed up the single-exchange case this work is
about.

There is a related and more tractable gap worth recording, though. The vector
backend covers `variable_base`, `straus`, `precomputed_straus` and `pippenger`,
but **not** fixed-base multiplication: `EdwardsPoint::mul_base` is serial in
every backend. That is 80% of `mul_base_clamped` (§13.4), i.e. of every ephemeral
key generation. Whether vectorizing it would pay is genuinely unclear — §13.5
shows the constant-time window scan, not the point additions, is what dominates
there, and a scan is a different thing to vectorize than an addition chain — but
it is the one place where the existing vector backend has an obvious hole on
this hot path.

For what it is worth, the vector backend demonstrably works where it is wired
in: §1 measures `vartime_double_scalar_mul_basepoint` at 38 491 ns under
`"simd"` against 47 949 ns under `"serial"`. So this is a coverage gap, not a
dead backend.

---

### 13.4 Where `mul_base_clamped` spends its time

The other X25519 operation on a libsignal-style hot path is ephemeral key
generation, which is `MontgomeryPoint::mul_base_clamped`. It splits cleanly in
two, and the parts add up to the whole (min of 15, ns):

| | x86_64 | wasm32 |
| --- | ---: | ---: |
| `EdwardsPoint::mul_base` (fixed-base, radix-16 table) | 15 012 | 44 080 |
| `EdwardsPoint::to_montgomery` (the conversion) | 4 038 | 9 990 |
| — of which `FieldElement::invert` | 3 766 (93%) | 9 280 (93%) |
| sum of the parts | 19 050 | 54 070 |
| measured `mul_base_clamped` | **18 789** | **52 811** |

So ~80% is the fixed-base Edwards multiplication and ~21% is a single field
inversion, `to_montgomery` being `(Z+Y) * (Z-Y).invert()`.

The inversion is not doing anything wasteful: `invert` is Fermat, 254 squarings
and 11 multiplications, and 254 x 14.10 + 11 x 25.36 = 3 860 ns predicts the
measured 3 766 to within 2.5%. It is optimal *as an exponentiation*.

### 13.5 Bigger basepoint tables do not help, on either target

`EdwardsPoint::mul_base` uses the 30 KB radix-16 table. The crate also exposes
radix-32/64/128/256 tables as public API, documented as needing fewer additions
(64 → 52 → 43), so a consumer chasing fixed-base throughput would reasonably try
one. Best of three runs, minimum of 15 repetitions each, ns:

| table | size | additions | x86_64 | wasm32 |
| --- | ---: | ---: | ---: | ---: |
| radix-16 (what `mul_base` uses) | 30 KB | 64 | **15 229** | **45 039** |
| radix-32 | 60 KB | 52 | 15 369 (+0.9%) | 45 655 (+1.4%) |
| radix-64 | 120 KB | 43 | 17 302 (+13.6%) | 54 066 (+20.0%) |

Doubling the table to radix-32 is a wash — the difference is inside this host's
run-to-run spread — and quadrupling it to radix-64 is clearly worse. So there is
no reason to move off the default, and a consumer who does pays in binary size
and, at radix-64, in time as well.

The addition count is the wrong thing to optimize, because the lookup has to be
constant time: selecting one entry from a window is a linear scan over all of
them. Radix-64 cuts additions by a third but takes the scan from 64 × 8 = 512
conditional selects to 43 × 32 = 1376, each over a three-field-element
`AffineNielsPoint`. The scan wins the trade.

> **Correction.** An earlier revision of this document reported these
> differences as +8.3%/+30.7% on x86_64 and +37%/+115% on wasm32. Those figures
> were an artifact of the harness: the radix-32 and radix-64 tables were being
> constructed *inside* the timed kernel, so a table build — on the order of a
> thousand point operations — was amortized into every repetition, and worst
> where the calibrated iteration count was smallest. The tables are now built
> once behind a `OnceLock`, outside the timing boundary. The ordering was
> unaffected; the magnitudes were not. Thanks to CodeRabbit for catching it.

### 13.6 Exact instruction profile, and what it says is left

Everything above is wall-clock on a shared virtual machine, which is why the
minimum over repetitions is the reported statistic. Instruction counts have no
such problem, so the question "where does the work actually go" is answered
exactly rather than statistically, with `callgrind` (`benches/harness/src/bin/cg.rs`
is the entry point; the command is in its header).

Instructions per operation, stock release build (`lto=off`, `cgu=16`, baseline
x86-64 ISA — so `mulq` rather than `mulx`), averaged over 20 calls:

**`mul_clamped`** — 600 000 instructions:

| | share |
| --- | ---: |
| `<&FieldElement51 as Mul>::mul` | **46.4%** |
| `differential_add_and_double` (self: the four inlined squarings, `mul121666`, the eight add/sub) | **40.2%** |
| `pow2k` (the inversion's 254 squarings) | 5.2% |
| ladder driver, self (conditional swap, bit extraction) | 4.9% |
| `pow22501` | 0.1% |

**`mul_base_clamped`** — 226 000 instructions:

| | share |
| --- | ---: |
| `<&FieldElement51 as Mul>::mul` | **46.2%** |
| `AffineNielsPoint::conditional_assign` (the constant-time window scan) | **15.0%** |
| `pow2k` (the inversion's squarings) | 14.1% |
| `EdwardsPoint + AffineNielsPoint` (the mixed addition) | 6.3% |
| `LookupTable::select` (self: the scan's loop) | 4.8% |
| `ProjectivePoint::double` | 1.3% |

Three things worth pulling out.

**The op-count model checks out exactly.** `mul` costs 284 427 instructions per
`mul_clamped` over 1287 calls = 221 each, which is precisely the 221-instruction
`mul` counted in the disassembly in §3.1. That is an independent confirmation of
the operation counts in §1 from a completely different instrument.

**The inversion is smaller than the wall-clock suggested.** It is 5.3% of
`mul_clamped`'s instructions, against the ~8% estimated from timing, and about
20% of `mul_base_clamped` once its share of `mul` is included. That is still the
largest single algorithmic lever, but it is worth recording that the timing
estimate flattered it.

**The constant-time window scan is a fifth of key generation** — 15.0% in
`conditional_assign` plus 4.8% in `select`, against 6.3% for the point addition
the scan exists to feed. This is the quantitative version of §13.5: the scan, not
the addition count, is what dominates fixed-base multiplication, which is why
larger tables lose.

### 13.7 The ladder is register-pressure bound, not schedule bound

§13.6 shows `differential_add_and_double`'s self cost is 966 instructions per
step while the arithmetic in it accounts for about 820. The rest is register
pressure: the step holds up to six live field elements — thirty `u64` — against
sixteen general-purpose registers, so the allocator spills. And the step runs in
181 ns against a 197 ns sum of isolated operation latencies, so only about 8% of
the available instruction-level parallelism is being recovered.

Both of those look like an invitation. Three separate attempts to take it up all
failed, and they failed the same way, which is the useful part.

**Attempt 1: reschedule the step.** Its statements are pure dataflow, so any
topological order computes the same values — reordering is risk-free and
`callgrind` measures it exactly. Instructions for 20 `mul_clamped`:

| schedule | `differential_add_and_double` | program total |
| --- | ---: | ---: |
| as written | 4 926 600 | 12 269 825 |
| A: doubling half completed first | 4 834 800 | 12 176 959 |
| **B: doubling half first *and* `t2`/`t3` sunk to their use** | **4 753 200** | **12 095 359** |
| C: both multiplications first, doubling tail last | 4 819 500 | 12 161 659 |

B removes **173 400 instructions, −1.42% of the program**, 34 per step. (A is
worse than B because hoisting `t2`/`t3` only trades `t4`/`t5`'s live range for
theirs.) On wasm32 it then measured **124 359 → 126 167 ns, +1.5%, slower in all
five paired runs**; on x86_64 the reduction is below the 2.5% noise floor and
unmeasurable. Reverted.

**Attempt 2: fuse the independent pairs.** The step has three pairs of mutually
independent operations — `t4`/`t5`, `t7`/`t8`, `t11`/`t12` — and `mul` is an
outlined call, so each pair is separated by a call boundary that an
out-of-order core can only see across as far as its reorder buffer reaches.
Factoring `mul`'s body into an inlinable `mul_limbs` and adding `mul_pair` /
`square_pair`, which inline two bodies into one function so the scheduler can
interleave them freely, is bit-for-bit identical by construction. Result:
**+1.0% instructions** (the fused bodies spill more), **−0.36% on x86_64** —
inside the noise floor — and **+3.6% on wasm32, slower in all four paired
runs**. Reverted.

**Attempt 3: inline `mul`.** Recorded in §3.3: isolated `fe_mul` −29.5%,
`mul_clamped` −0.5%.

The common thread is that every one of these *adds* live values, and the step is
already spilling. Scheduling freedom that costs register pressure is not a trade
this loop can afford, and the effect is worst on wasm32, where `serial::u32`
holds ten limbs per element instead of five and V8 does its own allocation on
top. The one change in this area that did work — forwarding the conditional swap
in §6 — went the other way: it *removed* work without adding a single live
value, and it is the only one of the four that measured a win.

That also settles the ILP question from the other side. The gap between 181 ns
and the 197 ns serial sum is not scheduling slack waiting to be claimed; it is
what the core already recovers on its own, and there is no cheap way to get more.

**A methodological note worth keeping:** x86-64 instruction count is not a proxy
for wasm32 time. Attempt 1 removed 1.42% of instructions and lost 1.5% on
wasm32. The targets run different backends with different liveness and different
compilers downstream. This is the same caution as §10's `bits="64"` result,
reached from the opposite direction, and it is why §6 was accepted on a measured
wasm32 win rather than on its instruction count.

### 13.8 wasm32: `simd128` is worth 4.6% on DH and 15.7% on key generation

`simd128` is a stable wasm feature, but it is **not** enabled by default for
`wasm32-unknown-unknown`. Turning it on costs nothing but a flag:

```sh
RUSTFLAGS='-C target-feature=+simd128'
```

Six alternating paired runs, minimum of 11-13 repetitions, ns:

| kernel | default | `+simd128` | change |
| --- | ---: | ---: | ---: |
| `fe_mul` | 50.55 | 50.94 | none |
| `fe_square` | 36.90 | 38.13 | none |
| `fe_invert` | 8 941 | 8 909 | none |
| **`x25519_mul_clamped`** | 126 265 | **120 446** | **-4.6%** |
| **`x25519_mul_base_clamped`** | 49 777 | **41 965** | **-15.7%** |

Faster in every paired run for both X25519 operations. The module also gets
*smaller*, 118 754 to 110 732 bytes.

**The field arithmetic does not move at all**, which is what makes this
interesting. `fe_mul`, `fe_square` and `fe_invert` are unchanged — LLVM does not
vectorize the radix-2^25.5 schoolbook, and that is where one might have expected
a SIMD flag to pay. The gain is entirely in the **constant-time selection
code**, and §13.6 says exactly why: `AffineNielsPoint::conditional_assign` is 15%
of key generation and `LookupTable::select` another 4.8%, and a masked select
over ten 32-bit limbs is the most vectorizable thing in the crate — four lanes
to a `v128`. That the gain is 15.7% on key generation, which is scan-dominated,
and 4.6% on the ladder, which does 256 conditional swaps but is otherwise
multiplication-bound, follows directly from that split.

**Verified before recommending it.** A build flag that changes generated code is
a correctness question, not just a speed one, so the two modules were run
side by side over all twelve kernels at three iteration counts each and their
returned checksums compared: identical in every case. The kernels chain their
output back into the next iteration, so a checksum match covers the full
sequence of field and point operations, not just a final value.

So this is the wasm32 counterpart of the `+bmi2` result in §3.5: a consumer-side
build flag worth more than any source change measured here, and for the same
reason — the target has a capability the code can use that the default feature
set does not expose. `simd128` shipped in Chrome 91, Firefox 89, Safari 16.4 and
Node 16, so for most deployments it is free; a consumer targeting older engines
should check their floor first.

### 13.9 x86_64: `+avx2` vectorizes the same scan, and is complementary to `+bmi2`

The wasm32 result in §13.8 raises the obvious question for the other target: the
constant-time window scan is 19.8% of key generation there too, so does x86_64
already vectorize it?

Not on the baseline ISA. `callgrind` on `mul_base_clamped`, 20 calls:

| | `AffineNielsPoint::conditional_assign` | program total |
| --- | ---: | ---: |
| baseline `x86-64` | 679 680 (15.0%) | 4 520 254 |
| `-C target-feature=+avx2` | **322 560 (8.1%)** | **3 996 414 (−11.6%)** |

AVX2 halves the scan. Baseline `x86-64` already has SSE2 and its two `u64`
lanes, but LLVM does not vectorize there and does at four lanes.

In time, best of three alternating runs, fat LTO, minimum of 15 repetitions:

| flags | `mul_clamped` | `mul_base_clamped` |
| --- | ---: | ---: |
| baseline `x86-64` | 43 595 | 16 729 |
| `+avx2` | 43 942 (none) | **15 996 (−4.4%)** |
| `+avx2,+adx,+bmi2` | **39 594 (−9.2%)** | **15 729 (−6.0%)** |

**The two flags help different operations, and neither substitutes for the
other.** `+bmi2` is what gets `mulx` into the field multiplication, so it moves
the ladder; `+avx2` is what vectorizes the constant-time lookup, so it moves key
generation. `+avx2` alone does nothing for `mul_clamped` — the ladder's swap is
two field elements, too little to vectorize profitably — and `+bmi2` alone does
comparatively little for key generation, which is scan-bound.

**`+adx` is not one of them, and was dropped from the recommendation.** The
measurements above were taken with `+avx2,+adx,+bmi2` because that is the
conventional incantation; a later check found the `+adx` half contributes
nothing. A whole-binary opcode census of the two flag sets is *identical* — 882
`mulx`, **0** `adcx`, **0** `adox` either way — and the times are
indistinguishable on every kernel across three paired rounds (`mul_clamped` min
37 698 with `+adx` against 37 814 without, `fe_mul` 26.05 against 26.10,
`mul_base_clamped` 15 464 against 15 515; the without-`+adx` build wins two of
three rounds). This is §3 confirmed from the deployment side rather than the
disassembly: 5×51 unsaturated limbs give each output limb a single carry, so
there is no second carry chain for ADX to interleave, and LLVM correctly emits
none. Since ADX shipped a generation after AVX2/BMI2 — Broadwell rather than
Haswell on Intel, Zen rather than Excavator on AMD — including it narrows the
deployable CPU set and buys nothing. The recommendation is `+avx2,+bmi2`.
Thanks to Codex for raising it.

Note that the −11.6% instruction reduction becomes −4.4% in time. That is the
expected direction: the scan streams 30 KB of table and is not purely
instruction-bound. It is also another instance of the §13.7 caution against
reading instruction counts as time.

One caveat that is a deployment decision rather than a measurement, and it
applies to **every** flag recommended here, not only `+avx2`. `-C
target-feature` lets the compiler emit those instructions unconditionally, with
no runtime check, so the binary's CPU floor rises to whatever it was told it
had: Haswell (2013) for `+avx2`, Broadwell (2014) for `+adx`, and for the
combination the later of the two. `-C target-cpu=native` is the same hazard in
sharper form, since it targets the build machine rather than a stated baseline —
fine for something compiled where it runs, wrong for a redistributed artifact.
On an older CPU the result is a fault, not a slowdown.

This is unlike the crate's own vector backend, which detects AVX2 at runtime and
falls back. A consumer who cannot raise the floor keeps the runtime-detected
backend for Edwards work, keeps `lto = "fat"` — which costs nothing in
portability and is the larger half of the win anyway — and simply does not get
the rest.

**Reproduction note:** `-C target-cpu=native` on this host emits AVX-512 that
valgrind 3.22 rejects with SIGILL, so the instruction counts above use explicit
`+avx2` rather than `native`. The timing runs are unaffected.

### 13.10 `mul`'s carry chain: isolated −16%, `mul_clamped` +2%. Refused

§7 succeeded by taking a dependency chain off the critical path, so the obvious
follow-up is the longest chain left in the hottest function. `mul`'s carry
propagation is strictly serial and six deep:

```
c0 → c1 → c2 → c3 → c4 → out[0] → out[1]
```

Because the limbs form a *ring* — `c4` folds back into `c0` with a factor of 19 —
the carries can be reordered into four levels with two independent carries each,
at the cost of one extra carry (seven instead of six):

| level | carries (independent within a level) |
| --- | --- |
| 1 | `0 → 1` and `2 → 3` |
| 2 | `1 → 2` and `3 → 4` |
| 3 | `4 → 0` (×19) and `2 → 3` again |
| 4 | `0 → 1` |

The bound analysis works out with room to spare — the final limbs are
`c0 < 2^51`, `c1 < 2^51 + 2^12.7`, `c2 < 2^51`, `c3 < 2^51 + 2^6.4`,
`c4 < 2^51`, all far inside the `b < 3` the next operation requires — and the
full suite passes, RFC 7748 vectors included.

**It makes isolated `mul` dramatically faster and `mul_clamped` slower.**

| | `fe_mul` | `mul_clamped` |
| --- | ---: | ---: |
| §7 only, baseline `x86-64` | 28.70 | 43 082 |
| §7 + ring carry | **23.97 (−16.5%)** | 42 547 (−1.2%, 2/3 runs) |
| §7 only, `+avx2,+adx,+bmi2` | 26.05 | 39 116 |
| §7 + ring carry | **23.85 (−8.4%)** | **39 928 (+2.0%, slower 3/3)** |

`fe_invert`, `fe_square` and `mul_base_clamped` do not move at all.
Instructions go **up**: +3.6% on the isolated kernel, +0.76% on `mul_clamped`.

**Refused**, by the brief's own rule that the number that counts is
`mul_clamped`'s: +2.0% on the configuration the crate recommends, losing all
three paired runs.

The reason is worth more than the result, because it is the sharpest form of the
§3.3 calibration and it points the opposite way from §7. **The isolated `fe_mul`
kernel is a serial chain — each multiplication consumes the previous one's
output — so it is latency-bound, and shortening the carry chain is exactly what
it wants. The ladder is not.** `differential_add_and_double` has three pairs of
mutually independent operations (§13.7), so an out-of-order core already overlaps
one multiplication's carry chain with the next multiplication's products.
Latency there is already hidden; the extra carry is not. The same intervention
is therefore worth −16% in one place and +2% in the other.

So `fe_mul` and `mul_clamped` are not merely weakly correlated, as §3.3 found —
here they are **anti-correlated**, and an optimizer trusting the isolated
microbenchmark would have shipped a regression. This is also why §7 was accepted
on `mul_clamped` rather than on its far more impressive isolated `fe_square`
number.

Not implemented; the working tree was reverted to the serial chain.

### 13.11 What is left, and why it was not attempted (one item since taken)



After §4 through §8, both hot loops are at their published operation counts and
the field operations are at the floor for this representation. Three separate
checks say so (the isolated `mul`/`square` timings quoted below predate §7, so
the *ratio* is the claim here, not the absolute figures):

* **The field operations themselves.** `mul` is 25.36 ns for 25 partial products
  and `square` is 15.03 ns for 15; the ratio 1.687 is within 1.2% of the
  25/15 = 1.667 the product counts predict. Neither has meaningful slack short
  of changing the limb layout.
* **The ladder step.** `differential_add_and_double` now performs
  5M + 4S + 1×a24 — five general multiplications (`t7`, `t8`, `t14`, `t16`,
  `t17`), four squarings, and the specialized multiply by `(A+2)/4` from §5.
  That is the textbook cost of a differential add-and-double; before §5 it was
  6M + 4S, paying a general multiplication for the constant.
* **The fixed-base loop.** Each of the 64 windows costs one mixed addition
  (`&EdwardsPoint + &AffineNielsPoint`, 3M) plus `CompletedPoint::as_extended`
  (4M) = 7M, which is the standard cost for a mixed addition against
  Niels-form precomputed points. 64 × 7 × 25.36 ns = 11.4 µs against a measured
  15.0 µs for `mul_base`, the remainder being the constant-time window scan and
  the four doublings — consistent with §13.5, where the scan is what makes the
  larger tables lose.

  > **This bullet was wrong, and §13.12 is what it missed.** Every operation
  > count in it is correct, and the conclusion drawn from them — that the
  > fixed-base loop is at its floor — was off by **45%**. The error is in the
  > `25.36 ns`: that is the cost of a `serial::u64` multiplication, and the
  > analysis never asked whether that was the multiplication being used. It was,
  > and it did not have to be. Pricing 7M per window at the serial rate assumes
  > the backend the loop happens to run on, and an operation count cannot see an
  > assumption like that — it prices the operations it is given.
  >
  > Worth keeping as written, because "both hot loops are at their published
  > operation counts" was true the whole time and is exactly why this went
  > unexamined for eleven sections.

The one substantial remaining lever is the **field inversion**: 3 766 ns on
x86_64 and 9 280 ns on wasm32, which §13.6 pins down exactly as 5.3% of a
`mul_clamped`'s instructions and about 20% of a `mul_base_clamped`'s. Fermat's little theorem is optimal as an exponentiation, but
it is not the only algorithm: a constant-time binary GCD in the style of
Bernstein–Yang "safegcd" typically runs 2–4x faster than Fermat for a 255-bit
field. At 2.5x that would be worth roughly −5% on `mul_clamped` and −12% on
`mul_base_clamped`, or about −6% of the per-message X25519 cost.

It is also worth recording that fiat-crypto, which this crate already carries
behind a cfg, does **not** generate `divstep` for curve25519 — its requested
operation list is `carry_mul, carry_square, carry, add, sub, opp, selectznz,
to_bytes, from_bytes, relax, carry_scmul121666`. So there are no verified
primitives to build on here; a safegcd would be hand-written from scratch.

**Not attempted here.** safegcd is subtle, its constant-time property is a
property of the divstep bound and the whole iteration count rather than of any
line in isolation, and getting it wrong in a way tests do not catch is exactly
the failure mode this kind of code cannot afford. It would need its own
correctness argument, its own differential testing against Fermat across the
full input range, and its own constant-time review — a change of a different
character from the two in §4 and §5, both of which move existing arithmetic and
are bit-for-bit verifiable against what they replace. Recorded as the top
candidate for anyone who wants to take it on, with the numbers that say what it
is worth.

Two smaller candidates that the instruction profile surfaced, and which were
measured and left alone rather than merely skipped:

* **XOR-accumulating the window scan.** `AffineNielsPoint::conditional_assign`
  is `a ^ (mask & (a ^ b))` per limb — three operations. Accumulating into a
  zeroed register instead would be `acc ^= mask & b`, two. Against 15.0% of
  key generation that is worth roughly 3%, but it needs an extra pass to
  reinstate the identity when the window digit is zero, and it is shared
  constant-time code on every scalar multiplication in the crate, not just the
  X25519 path. Poor trade.
* **An unreduced subtraction.** ~~Left alone~~ — **this was subsequently done,
  in §8.** The estimate here was that dropping the reduction was worth about 3%
  of `mul_clamped`'s instructions; it measured **−6.58%**, and −6.1% of the time
  on a baseline `x86-64` build. The reason for leaving it alone was that it
  needs a narrower bit-excess precondition, and that contract is a correctness
  invariant rather than a detail. That objection was answered by *adding* an
  operation with its own stated precondition instead of narrowing `Sub`'s, which
  leaves the existing contract untouched. The entry is kept rather than deleted
  because the reasoning that deferred it, and what changed to unblock it, is the
  useful part.

Neither was refused on principle; both were recorded with what they were worth
so the trade could be re-made — and one of them since has been.

---

### 13.12 `mul_base` now reaches the vector backend: −45.1%

`edwards_mul_base` used to cost an identical 153 913 instructions under
`curve25519_dalek_backend = "serial"` and `= "simd"` — not close, identical. The
fixed-base ladder ran entirely in `serial::u64` even when AVX2 was live, because
the vector backend had `variable_base` and `vartime_double_base` and no
fixed-base path at all. The same pair of builds moved `vartime_double_base` by
−47.3% and `ed25519_verify` by −43.9%, so the backend was not short of headroom
here; it simply was not wired in.

Where the 153 913 went:

| | share |
| --- | ---: |
| `serial::u64::FieldElement51::mul` | **63.92%** |
| `LookupTable<AffineNielsPoint>::select_or` (the constant-time scan) | 23.03% |
| `EdwardsPoint + AffineNielsPoint` | 5.75% |
| `mul_base` itself | 2.47% |
| `ProjectivePoint::double` | 1.90% |
| `subtle::black_box` | 1.23% |

About three quarters was field arithmetic the vector types replace.

**The change.** `backend/vector/scalar_mul/fixed_base.rs` runs the same radix-16
ladder — 64 additions and one `mul_by_pow_2(4)`, indexed identically — over
`ExtendedPoint` and `CachedPoint`, reusing the `LookupTable<CachedPoint>` and
constant-time `select` that vector `variable_base` already had. The only new
thing needed was the table: `BASEPOINT_TABLE`, 32 sub-tables of eight
`CachedPoint`s, **40 960 bytes exactly**, generated and then checked against the
basepoint rather than trusted.

`EdwardsPoint::mul_base` now goes through `backend::mul_base`, which dispatches
on `get_selected_backend()` — the same runtime `cpuid` check every other vector
entry point uses. That matters: `curve25519_dalek_backend = "simd"` means "a
runtime-dispatched backend was compiled in", **not** "this CPU has AVX2", and
calling the AVX2 ladder on the strength of the cfg alone would be a `SIGILL` on
an x86_64 machine without it.

| | instructions | wall clock |
| --- | ---: | ---: |
| `edwards_mul_base`, serial ladder | 153 910 | 12 579 ns |
| `edwards_mul_base`, vector ladder | **84 546** | **8 444 ns** |
| | **−45.1%** | **−32.9%** |
| `x25519_mul_base_clamped` | | 16 256 → **12 349 ns (−24.0%)** |

The wall-clock gain is smaller than the instruction gain, which is what a shift
to 256-bit lanes looks like: fewer instructions, lower throughput each.

**What it costs.** 40 KB of `.rodata` on AVX2 builds, on top of the serial
table's 30 KB, which every non-AVX2 build and the public
`EdwardsBasepointTable` API still need. Both tables ship in an AVX2 build; this
is a size-for-speed trade, not a replacement.

**`avx512` keeps the serial ladder.** `CachedPoint` has a different limb layout
under `ifma`, so it would need a second generated constant with different
contents. Under `curve25519_dalek_backend = "avx512"` the dispatcher's AVX2 arm
is not compiled in the first place, so that build is unaffected either way.

**Correctness.** Three tests, in `fixed_base.rs`:

* `vector_fixed_base_matches_serial` — 64 consecutive scalars against
  `EdwardsPoint::mul_base`. Unlike §9.3's batch inversion this is expected to
  agree exactly: no inverse representatives are involved, only the same group
  operations in the same order.
* `vector_fixed_base_handles_edge_scalars` — zero, one, `-1`, and `[0xff; 32]`,
  whose top radix-16 digit is the one `as_radix_2w(4)` has to carry into.
* `basepoint_table_constant_holds_the_right_multiples` — the generated constant
  is checked entry by entry: sub-table `i` must hold \( k (2^{8})^{i} B \) for
  every `k` it can select, all 256 of them. A generated table that is wrong in
  one entry would still pass a handful of scalar comparisons.

Who benefits is worth stating plainly, because it is not verification:
`mul_base` is Ed25519 key generation and **signing**, and X25519 base-point
clamping. Verification goes through `vartime_double_scalar_mul_basepoint`, which
was already vectorised. A workload dominated by verification does not notice
this.

### 13.13 The whole branch, end to end

§6.4 certified §4 + §5 + §6 and said so explicitly; §7 onwards were each
measured against the tree in front of them. Nine more changes landed after that,
so the cumulative claim is worth taking once more, directly rather than
assembled from per-change deltas.

Both arms use the **same harness**, differing only in the library: the base
restored with `git checkout a3549477 -- curve25519-dalek/src`, `HEAD` for the
branch, `-C target-feature=+avx2` on both. Instruction counts are exact, with
fixed setup differenced out by the usual long-run-against-short-run subtraction.

| kernel | `main` (a3549477) | this branch | change |
| --- | ---: | ---: | ---: |
| `x25519_mul_clamped` — DH | 662 659 | **449 330** | **−32.2%** |
| `x25519_mul_base_clamped` — keygen | 179 696 | **116 266** | **−35.3%** |
| `edwards_mul_base` — fixed base | 146 601 | **84 546** | **−42.3%** |
| `ed25519_verify` — verification | 331 696 | **306 124** | **−7.7%** |
| `edwards_vartime_double_base` | 289 626 | **265 538** | **−8.3%** |
| `edwards_table_create` | 10 011 269 | **2 641 096** | **−73.6%** |
| `edwards_to_montgomery` | 37 789 | **36 312** | −3.9% |

Wall clock, minimum of 15 repetitions, best of five alternating rounds:

| | `main` | this branch | change |
| --- | ---: | ---: | ---: |
| DH | 58 188 ns | **38 095 ns** | **−34.5%** |
| keygen | 16 222 ns | **12 276 ns** | **−24.3%** |
| verification | 33 230 ns | **31 727 ns** | −4.5% |

**Verification moves least, and that is the honest shape of the result.** It was
already the best-served path in the crate — 87.4% AVX2, with the vector backend
wired in — so what was left there was the backend's own inefficiencies (§10) and
nothing structural. DH, keygen and fixed base moved most because they were not
served at all: the ladder never reaches the vector backend by construction
(§13.3), and `mul_base` did not reach it for no reason other than that nobody
had written the path (§13.12).

These numbers do not compose with the per-change figures elsewhere in this
document and should not be added to them. Several were measured on different
hosts across several sessions; this table is one host, one session, two binaries.

> **Two properties of the build these were taken under, both worth stating.**
>
> It carries `+avx2` but **not `+bmi2`**, so the binary contains **zero `mulx`**
> — 1002 plain `mul`. §3.1 measures `+bmi2` at `mul` 221 → 188 and `pow2k`
> 152 → 129, and every AVX2-capable CPU has BMI2, so the absolute figures here
> are roughly 15% above what the build this document *recommends* produces. The
> relative comparisons are unaffected — both arms carry the same flags.
>
> More consequentially, compile-time `+avx2` makes `cpufeatures`' check
> const-fold, so `backend::mul_base`'s serial arm is dead-code-eliminated:
> `nm` on the binary shows the vector `BASEPOINT_TABLE` and **no**
> `ED25519_BASEPOINT_TABLE`, and `mul_base` contains no `cpuid`. The −45.1% in
> §13.12 was therefore measured with the dispatch free. A consumer building for
> runtime dispatch — which is the point of the vector backend — pays for the
> check and links both tables. That cuts both ways: it means the speed figure
> is a mild overstatement and the **size** figure is a large one, since `+avx2`
> drops the 30 720-byte serial table outright.

### 13.14 `MontgomeryPoint::to_edwards` ran two Fermat exponentiations

§9.2 established that a *verification* performs exactly one field inversion and
concluded that "no amount of restructuring inside this crate changes that." The
inventory was right about verification. It never audited the Montgomery-to-
Edwards conversion, which is the XEdDSA entry point, and that function was
running the 250-squaring chain **twice**:

```rust
let y = &(&u - &one) * &(&u + &one).invert();   // chain #1
CompressedEdwardsY(y.to_bytes()).decompress()   // sqrt_ratio_i: chain #2
```

Decompression uses `y` only through a ratio — it wants `x` with
`x² = (y²−1)/(dy²+1)` — so substituting `y = yn/yd` and scaling numerator and
denominator by `yd²` gives `x² = (yn²−yd²)/(d·yn²+yd²)` with no division. And
`EdwardsPoint` is projective, so `yn` and `yd` are its `Y` and `Z`.

| | instructions |
| --- | ---: |
| before | 103 616 |
| after | **71 191** |
| | **−31.29%** |

Nothing had measured this: `edwards_to_montgomery` is the *other* direction, and
the profile set contained no kernel for this one. `montgomery_to_edwards` now
exists so it cannot drift unwatched.

**Why the earlier audit missed it.** §9.2's sentence generalises past its
evidence. It asks whether the inversion can be *batched* — correctly answering
no, one inversion cannot be batched with itself — and §9.5 and §13.11 then only
ever ask whether it can be made *faster*. The option of removing one was never
on the table, and here one was simply redundant. "Inside this crate" was also
doing work in that sentence that it should not have: the redundancy is visible
only when you look at the caller and the callee together.

**Correctness.** The result is congruent but not limb-for-limb identical — `Z`
is now `u+1` rather than `1` — so it is checked on canonical bytes against a
verbatim copy of the previous implementation over 256 `u` values and both signs,
asserting agreement where both accept, agreement on rejection where both reject,
and that both outcomes occur at all. `d·yn² + yd²` cannot vanish: it would need
`(yn/yd)² = −1/d`, and `−1/d` is a nonsquare because `−1` is square mod p and `d`
is not; `yd = 0` means `u = −1`, which the function already rejected.

## 14. Validation

* **A 12-cell feature x backend matrix, under `-D warnings`.** Three feature
  sets (default, `--no-default-features --features alloc`, `--all-features`)
  against four backend selections (default, `serial`, `fiat`, `bits="32"`), plus
  `+avx2` and `no_std` thumbv7em builds of all three crates.

  **That matrix exists because this class of bug has now bitten three times, and
  the third time it was the same mistake as the first.** §11.4 was green locally
  and broke seven CI jobs on dead code: `select_or` is called only by the
  basepoint tables, which are `precomputed-tables`-gated, so with that feature
  off nothing referenced it and `-D warnings` promoted the dead-code warning to
  an error. Every local run had the feature on — which is verbatim what change 6
  did, and the note written after change 6, quoted just below, was correct and
  got ignored.

  A subtler variant hid in the same fix: the local command was `cargo test`,
  CI's is effectively `cargo test` under `-D warnings`, so a duplicated
  `use super::*;` was invisible here and fatal there. **Match CI's strictness,
  not only its targets.** The widened matrix also caught a wart of this branch's
  own making: after change 2 moved the ladder off `APLUS2_OVER_FOUR`, the
  constant's only surviving reference was a test in `serial/u64/field.rs`, which
  the `fiat` build does not compile — dead there, and tolerated only because
  CI's `fiat` job does not pass `-D warnings`. Its `cfg` is now precise.

* **Full test suite**, `--all-features`, on every backend path — with and
  without the new cfg, since a path only tested when enabled is not tested:

  | configuration | result |
  | --- | --- |
  | default (`simd` / `bits=64`) | 156 passed, 21 doctests |
  | `curve25519_dalek_backend="serial"` | 156 passed, 21 doctests |
  | `curve25519_dalek_backend="fiat"` | 148 passed, 21 doctests |
  | `curve25519_dalek_bits="32"` | 151 passed, 21 doctests |
  | `curve25519_dalek_bits="32"` + `fiat` | 148 passed, 21 doctests |
  | `-C target-feature=+avx2,+adx,+bmi2` | 166 passed, 21 doctests |
  | `--cfg curve25519_dalek_bench_internals` | 156 passed, 21 doctests |
  | `--release` (debug assertions off) | 155 passed, 21 doctests |

  The counts differ by configuration because each selects a different backend
  module, and the backends carry different numbers of their own unit tests;
  `+avx2` additionally enables the AVX2 backend's tests, which the baseline ISA
  cannot run.

* **Differential tests** (`src/backend/serial/u64/field.rs`, `mod test`). The
  pre-refactor `pow2k` is kept verbatim as `reference_pow2k` and the new code is
  compared against it **limb for limb**, not merely as a field element — the
  limbs feed the next operation's bit-excess accounting, so field equality would
  not be a strong enough check.

  * `square_matches_reference_bit_for_bit` — 512 random elements at each of
    51/52/53/54 bits per limb, plus all-limbs-`2^54 − 1` (the top of the
    documented bit excess `b < 3`), plus zero and unit limb patterns.
  * `pow2k_matches_reference_bit_for_bit` — k ∈ 1..=8 and k ∈ {10, 20, 50, 100}
    (the values the inversion actually uses), at 51 and 54 bits.
  * `chained_square_matches_pow2k` — k chained `square()`s equal one `pow2k(k)`
    for k ∈ 1..=16.
  * `square_agrees_with_mul_at_bit_excess_bound` — an independent cross-check
    that does not use the reference copy: `x.square() == &x * &x` and
    `x.square2() == x² + x²`.
  * `mul_and_square_are_consistent` — `(x + y)² == x² + 2xy + y²`, exercising
    `mul` and `square` against each other. `mul` is unchanged here, so this is a
    regression guard for any future change to the multiplication.

  These same tests are what certify the squaring change in §7: the reference
  copy still contains the original `2*(u128)` form, so
  `square_matches_reference_bit_for_bit` is a genuine differential check of the
  new doubling against the old one, at every width up to the top of the
  documented bit excess. `mul121666`'s tests below are the model for why the
  reference is kept verbatim rather than rewritten.

  For `sub_unreduced` (§8) the check is different in kind, because the result is
  a different representative rather than the same limbs:

  * `sub_unreduced_agrees_with_sub` — 1024 operand pairs *built the way the
    ladder builds them*, by multiplying and squaring random elements rather than
    assuming what those outputs look like, asserting both `to_bytes()` equality
    with the general subtraction and that every limb stays under `2^53` — the
    documented postcondition, not the looser `2^54` that `b < 3` would accept.
    Plus the worst case the precondition permits, the degenerate inputs, and the
    ladder's *first* iteration, whose operands are the identity `(1, 0)` and
    `(from_bytes(u), 1)` rather than products — the one case the "both operands
    are `mul`/`square` outputs" phrasing does not literally cover.
  * The precondition's `debug_assert!`s were exercised by the entire suite in a
    debug build, ladder included; a violation anywhere would have panicked.

  For `mul121666` (§5), in both `serial::u64` and `serial::u32`:

  * `mul121666_matches_general_mul` — limb for limb against
    `&x * &APLUS2_OVER_FOUR` on 512 random elements at each width up to the top
    of the documented bit excess, plus the all-ones and degenerate limb
    patterns. This is the specialization's correctness condition.
  * `mul121666_is_multiplication_by_121666` — double-and-add over the bits of
    121666 using only `add` and `reduce`, which share no code with either
    `mul121666` or `mul`. This is what would catch the two being wrong the same
    way.

* **The harness links this crate, checked rather than assumed.** After §9.7,
  the harness manifest patches `curve25519-dalek` itself, and the lockfile
  carries exactly one entry for it with no `source = "registry+…"`. Worth
  re-checking whenever a bench-only dependency is added: any crate that depends
  on `curve25519-dalek` by version will silently pull a second copy.

* **Ed25519 reference vectors.** `window.rs` changed, and `LookupTable` is what
  every fixed-base multiplication selects from, so `ed25519-dalek`'s
  `against_reference_implementation` was run explicitly: the 128 vectors of
  `sign.input` from `ed25519.cr.yp.to`, sign and verify, all pass. Its full suite
  passes too, `--all-features` included, as does the `ed25519ph` RFC 8032
  prehash vector.

* **Batched table construction** (`src/window.rs`, `mod test`), three tests
  described in §9.3: `batched_lookup_table_matches_reference` against a verbatim
  copy of the pre-batch construction over 32 base points and all three
  coordinates of all eight entries, `batched_lookup_table_handles_identity`, and
  `created_table_agrees_with_precomputed`. The comparison is on canonical bytes
  rather than limbs, because Montgomery's trick returns a different
  weakly-reduced representative; §9.3 says why that is the strongest true
  statement here.

* **RFC 7748 known-answer tests.** `montgomery.rs` changed, so the ladder's own
  vectors were run from the dependent crate: `x25519-dalek`'s
  `rfc7748_ladder_test1_vectorset1`/`vectorset2` and the four
  `rfc7748_diffie_hellman` vectors pass, and so does `rfc7748_ladder_test2` —
  the 1000-iteration iterated-ladder vector that is `#[ignore]`d for cost — run
  explicitly in release (56.8 s, and again at 45.2 s after §8, which changes the
  ladder's arithmetic representative and so leans on this vector hardest).
  `ed25519-dalek`'s suite passes too.

* **Constant time.** `square_limbs`, `mul121666` and `sub_unreduced` are
  straight-line code: no branch, no data-dependent memory index, no early
  return. `sub_unreduced` in particular is strictly less work than the `Sub` it
  replaces, which was itself branch-free. Verified in the generated
  assembly — the emitted `mul` contains zero `j*`/`cmov`/`set*`; `pow2k` retains
  exactly one `jne`, which is its `k` loop, and `k` is a public exponent-chain
  constant, never secret. `square` no longer appears as an outlined symbol
  because it is a single squaring the optimizer folds into its callers. The
  refactor moves existing arithmetic without adding a conditional. After §5 the
  whole of `differential_add_and_double` still disassembles with an empty
  branch/`cmov`/`set*` set.

  For §7 this was checked as a whole-binary census across the two harness
  binaries, which differ only in that file:

  | | `shld` | `shrd` | branches | `cmov` | `set*` |
  | --- | ---: | ---: | ---: | ---: | ---: |
  | baseline | 353 | 54 | 8 912 | 573 | 446 |
  | §7 | **268** | 54 | 8 912 | 573 | 446 |

  Every column that could carry a secret-dependent decision is identical, and
  the one that moved is the intended one: 85 fewer `shld`, which is exactly 5 per
  inlined copy of `square_limbs` across the 17 copies in the binary — the five
  128-bit doublings, and nothing else.

* **wasm32 builds**: `cargo build --target wasm32-unknown-unknown` succeeds for
  the default, for `--cfg curve25519_dalek_bits="64"`, and for
  `--no-default-features`. For §7, the stronger statement: building the harness
  module from the same directory with only `u64/field.rs` swapped between the
  two versions gives a **byte-identical** `.wasm`, so that change cannot have
  moved wasm32 in either direction.

* **`cargo fmt --all`** clean.

* **`cargo clippy --all-targets -- -D warnings`**: 11 errors in
  `curve25519-dalek`, all pre-existing on the base branch (`edwards.rs` ×6,
  `ristretto.rs` ×3, `montgomery.rs` ×2 — all in `mod test`), plus 2 in
  `curve25519-dalek-derive`'s own tests. Re-verified after §7 by running the
  identical command against a pristine worktree of the base commit: same counts,
  same files. Note that at the workspace root the derive crate's failure can
  stop the build before `curve25519-dalek`'s test targets are checked, so the
  per-package form is the one that actually exercises them. The new harness
  crate is clippy-clean, `--all-targets` and default alike.

---

## 15. Summary

### 15.1 The whole branch, certified against `origin/main`

§6.4 certified §4+§5+§6 when those were all there was. With §7 and §8 landed,
the cumulative claim was re-measured the same way and end to end: `origin/main`
checked out into a second worktree, the same harness copied in, both trees
measured **in one session, alternating, on the same pinned core**. Minimum of 15
repetitions, three alternating rounds; `mul_clamped` in ns.

| profile | `origin/main` | this branch | change | paired wins |
| --- | ---: | ---: | ---: | :---: |
| x86_64, `lto=off`, `cgu=16` (cargo `--release` default) | 54 815 | **49 306** | **−10.1%** | 3/3 |
| x86_64, `lto=off`, `cgu=16`, `+avx2,+bmi2` | 54 037 | 48 967 | −9.4% | 3/3 |
| x86_64, `lto=fat`, `cgu=1` | 56 428 | **40 234** | **−28.7%** | 3/3 |
| x86_64, `lto=fat`, `cgu=1`, `+avx2,+bmi2` | 53 388 | **37 235** | **−30.3%** | 3/3 |
| wasm32, `serial::u32` | 137 213 | **124 191** | **−9.5%** | 3/3 |

`mul_base_clamped` over the same runs: 17 870 → 17 573, 17 021 → 16 620,
16 872 → 16 393, 15 694 → 15 510, and on wasm32 49 317 → 49 331. So key
generation is now slightly *better* rather than flat — between 1% and 3% on
x86_64 — which is §7 reaching the field inversion in `to_montgomery`; nothing
here touches the fixed-base path itself.

Two notes on reading this table. The `origin/main` column varies by about 6%
across profiles measured minutes apart on the same host, which is the honest
size of the run-to-run environment on a shared virtual machine and the reason
every figure in this document is a minimum over repetitions taken from an
alternating pair. And the `lto=fat` rows improve far more than the `lto=off`
rows for the reason §4 gives: only fat LTO inlines the ladder's squaring into
the step, which is what makes §4 and §7 pay.

### 15.2 Front by front

| front | outcome |
| --- | --- |
| **A — ADX/BMI2 asm on x86_64** | **Refused.** LLVM emits all 25/25 and 15/15 available `mulx` under `+bmi2`; `adcx`/`adox` have no second carry chain to run in a 5×51 two-word accumulator, and their payoff belongs to a 4×64 saturated layout that is out of scope. Calibration: an intervention giving isolated `mul` −30% moved `mul_clamped` by 0.5%. |
| **wasm32 `+simd128` (no code)** | `-C target-feature=+simd128` is worth **−4.6%** on `mul_clamped` and **−15.7%** on `mul_base_clamped`, and shrinks the module. Not on by default for `wasm32-unknown-unknown`. The field arithmetic does not vectorize; the whole gain is in the constant-time selection code (§13.8). |
| **A′ — build flags (no code)** | `+bmi2` moves the ladder (`mulx`), `+avx2` moves key generation (it halves the constant-time scan); together **−9.2%** on `mul_clamped` and **−6.0%** on `mul_base_clamped`, and with fat LTO ~−31% against stock `main`. **`+adx` was dropped from the recommendation**: identical instruction mix with and without it (882 `mulx`, 0 `adcx`, 0 `adox` either way), indistinguishable times, and it costs a CPU generation. Deployment findings for consumers (§13.9). |
| **B — wasm32 `bits="64"`** | **Refused.** 2.31× slower than the current default. wasm has no 64×64→128 multiply, so `u128` products are emulated. `build.rs`'s `TODO(Wasm32)` closed with evidence; behaviour unchanged. |
| **C — `square` via `pow2k(1)`** | **Changed.** `mul_clamped` −3.5% on a stock release build, −21.9% with fat LTO, −23.5% with fat LTO and `+adx,+bmi2`. Bit-for-bit identical output, no `unsafe`, no representation change, no API change, `serial::u32`/`fiat` untouched. |
| **D — ladder's multiply by 121666** | **Changed.** A specialized `mul121666` (fiat's verified `carry_scmul_121666` on the fiat backends) replaces a general multiplication whose operand had four zero limbs. Isolated 2.5x cheaper on x86_64, 3.8x on wasm32. `mul_clamped` −6.0% on a stock release build and −7.0% on wasm32; nothing under fat LTO, where the inliner already folded it. Complementary to C: between them every build profile improves. |
| **E — ladder's conditional swap** | **Changed.** `ProjectivePoint` inherited `subtle`'s default `conditional_swap` — a struct copy plus two conditional assignments — instead of forwarding to the masked exchange every field backend already implements. The ladder driver drops from 324 to 294 instructions per iteration. wasm32 **−2.2%** across five paired runs; on x86_64 the 0.6% difference is below a 2.5% noise floor measured from identical binaries. |
| **F — squaring's 128-bit doublings** | **Changed.** `square_limbs` doubled five 128-bit coefficients; `2*(x*y) == (2*x)*y`, so four precomputed 64-bit doublings cover all ten mirror-pair products instead — which is what `serial::u32` has always done. Isolated `fe_square` **−5.4%** (7/7 paired runs) and **−10.5%** with `+bmi2`; `mul_clamped` **−2.2%** with `+bmi2`, −0.7% stock; `fe_invert` −5.8%. The time win is several times the −1.29% instruction win because the 128-bit shift was on the dependency chain (§7). Bit-for-bit identical; wasm32 module byte-identical. |
| **G — the ladder's subtractions** | **Changed.** `Sub` adds `16p` and must then `reduce`, because it has to accept anything at the crate-wide `b < 3`. The ladder's four subtractions all take `mul`/`square` outputs, which are far narrower, so a separate `sub_unreduced` offsets by `2p` and needs no reduction — a new operation with its own stated precondition, not a change to `Sub`'s contract. `mul_clamped` **−6.1%** on baseline `x86-64` and **−3.7%** with `+avx2,+bmi2`, 3/3 paired runs each, −6.58% instructions; `mul_base_clamped` unchanged. Not limb-for-limb identical — a different representative — so it is checked on field equality, the limb bound, debug-assertions across the whole suite, and the 1000-iteration RFC 7748 ladder vector. `serial::u32` takes the same offset, where the bound closes with only **0.167 bits** (12%) of margin against a silent `u32` overflow, so `debug_assert!` checks that `19 * limb` still fits on every limb of every result rather than arguing it: **−4.39%** on `mul_clamped` there, the backend wasm32 uses (§8.4). |
| **H — batching the affine table conversions** | **Changed.** `LookupTable<AffineNielsPoint>::from` converted eight multiples one at a time, one field inversion each, so `EdwardsBasepointTable::create` did **256** — 91% of its cost. The chain depends on the previous multiple's *value*, not its affine form, so it runs in extended coordinates and converts all eight at the end with Montgomery's trick. `create` **−74.8% (3.97×)**, 7/7 paired runs, −72.1% instructions. Both hot paths unchanged: the crate ships its table as a constant. Not limb-for-limb identical — the batch returns a different weakly-reduced representative — so it is checked on canonical bytes against a verbatim copy of the old code, on the identity, and against the precomputed table (§9.3). |
| **I — tagging the AVX2 multiply per call site** | **Changed.** LLVM emitted one shared outlined body for the three `FieldElement2625x4` multiplies on the verification path; an unused `const N: u8` gives each site its own instantiation. `ed25519_verify` **−3.21% instructions** (316 280 → 306 124, callgrind, exact) for **+784 bytes** of `.text`. Wall-clock cannot resolve it — the host's spread was five times the effect — so the instruction count is the claim. The `Mul` operator stays as the untagged spelling; the tags are inert if a future LLVM stops sharing. A paired `negate_lazy` substitution was **refused** separately: −0.015%, and it would invalidate a documented `b < 0.007` bound (§10.5). |
| **J — vectorising `mul_base`** | **Changed.** `edwards_mul_base` cost an identical 153 913 instructions under the serial and simd backends: the fixed-base ladder never reached the vector backend, which had no fixed-base path. It now runs the same radix-16 ladder over `ExtendedPoint`/`CachedPoint`, against a new 40 960-byte `BASEPOINT_TABLE`. **−45.1% instructions** (153 910 → 84 546), **−32.9%** wall clock (12 579 → 8 444 ns), and `x25519_mul_base_clamped` **−24.0%**. Dispatched through `get_selected_backend()`'s runtime `cpuid`, not the backend cfg, which denotes a dispatched backend and not AVX2 hardware. Costs 40 KB of `.rodata` alongside the serial table, which the public API still needs; `avx512` keeps the serial ladder, its `CachedPoint` layout differing. The generated table is checked entry by entry against the basepoint (§13.12). |
| **K — the basepoint-table API missed the vector ladder** | **Changed.** `fc70cc9` wired the AVX2 fixed base into `backend::mul_base`, which only `EdwardsPoint::mul_base` calls; `&scalar * ED25519_BASEPOINT_TABLE`, the spelling this crate's own docs and `wacore-libsignal` use, kept the serial ladder at **153 910** instructions against **84 546**. §13.12's claim that signing benefits was false for the consumer that motivated it. `mul_base` now hands over when the table is the crate's own, recognised by address. `edwards_mul_base_table` is a permanent harness kernel so the two spellings cannot diverge silently again. |
| **L — operand-scanning the serial squaring** | **Changed.** `714dee6` did this for `mul` and left `square_limbs` grouped by output coefficient, where sixteen of eighteen spill stores are raw products. `fe_invert` **−3.83%**, `x25519_mul_base_clamped` −1.11%, `edwards_table_create` −1.75%, and **the ladder unchanged** — inlined there, the squaring was already well scheduled. Bit-for-bit identical; the differential test that certifies it already existed. |
| **M — hoisting the barrier on the *vector* scan** | **Refused, measured.** The vector scan never received §11.4's fix and is ~21% of `edwards_mul_base`; three independent readings of the disassembly predicted −5% to −8%. Two implementations measured **+5.98%** and **+2.50%**. The serial premise does not transfer: five `ymm` out of sixteen fit, so the crossings were not the binding constraint, and the hoist's bookkeeping costs more than it saves (§11.7). |
| **N — `to_edwards`'s second exponentiation** | **Changed.** `MontgomeryPoint::to_edwards` formed the Edwards `y` with a full Fermat inversion and then handed it to `decompress`, whose `sqrt_ratio_i` runs the same chain again. Decompression needs only the ratio, and `EdwardsPoint` is projective, so `yn`/`yd` become `Y`/`Z` and two squarings plus three multiplications replace the inversion: **103 616 → 71 191, −31.29%** on the XEdDSA verification entry point. Congruent, not limb-identical; checked against a verbatim copy of the old code over 256 inputs and both signs (§13.14). |
| **The verify kernel measured the wrong crate** | **Found and fixed (§9.7).** `ed25519-dalek` depends on `curve25519-dalek` by version, and the harness patched only `curve25519-dalek-derive`, so the binary linked two copies and `ed25519_verify` profiled the published 5.0.0 rather than this tree. It was blind to every change here — reporting "unchanged" for a regression as readily as for a win. Corrected: **354 575 → 330 457 Ir**, AVX2 87.8% → **87.4%**, the inversion 10.0% → **10.2%**. No conclusion in §9 changes, because they rest on the call graph rather than on these totals. The first attempt at the correction was itself contaminated by an uncommitted experiment and had to be re-taken in a pristine worktree; §9.7 records that too. |
| **Batching verification's inversions** | **Refused: there is nothing to batch.** An Ed25519 verification performs **exactly one** field inversion, in `compress`; the vector table build and the wNAF loop perform none, and the vector backend contains no runtime `invert` at all. The profile's sixteen `as_affine` cannot be in verification, and cannot be sixteen X25519 operations either — that would be four times the whole message's cycle budget (§9.2). |
| **safegcd, second look** | **Refused again.** Worth more here than in §13.11 — 10.0% of a verification, ~4% of the group client — but a 2–4× inversion caps the win at 2–3% of the client, while the crate's *existing* batch inversion is worth 6–10× wherever inversions co-occur. It also cannot be the variable-time kind, because `compress` is shared with secret-derived callers (§9.5). |
| **wNAF width / vartime tables** | **Refused analytically.** §13.5's rejection does *not* transfer — `NafLookupTable5::select` is a direct index, not a constant-time scan — but width 5 is already the minimum of build-plus-loop additions (49.7 against 51.6 at width 6), and the build is only 5.0% of a verification (§9.6). |
| **Combined C + D + E** | Superseded by the whole-branch certification in §15.1; kept because it is the only figure isolating these three. Certified against `origin/main` in one alternating session (§6.4) — F and G are *not* in it: x86_64 stock release **−9.7%**, stock + `+adx,+bmi2` **−10.0%**, fat LTO **−25.0%**, fat LTO + `+adx,+bmi2` **−25.3%**; wasm32 **−8.8%**. `mul_base_clamped` unchanged in every cell. |
| **Vector backend for Montgomery** | Viable but not worthwhile for single exchanges; the win would require a batched multi-exchange API. Not implemented. |
| **Bigger basepoint tables** | **Measured and rejected.** radix-32 is a wash against the default radix-16 (+0.9% x86_64, +1.4% wasm32) for twice the table size; radix-64 is +13.6% and +20.0% for four times. The constant-time window scan grows faster than the addition count falls. The crate's default is already right. |
| **Exploiting the ladder's spare ILP** | **Three attempts, all measured and reverted (§13.7).** Rescheduling the step removes 1.42% of x86-64 instructions but is 1.5% *slower* on wasm32; fusing the three independent operation pairs into single function bodies is +1.0% instructions, noise on x86_64 and 3.6% slower on wasm32; inlining `mul` gives −29.5% isolated and −0.5% end to end. All three add live values, and the step already spills — it is register-pressure bound, not schedule bound. |
| **Faster field inversion (safegcd)** | **Not attempted.** The top remaining lever — the inversion is 8% of `mul_clamped` and 20% of `mul_base_clamped`, and a constant-time binary GCD would plausibly be 2–4x faster than Fermat. Deferred because its constant-time property is global to the iteration rather than local, unlike the two changes above, which are bit-for-bit verifiable against the code they replace. |
