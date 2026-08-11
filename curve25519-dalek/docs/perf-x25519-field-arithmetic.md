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

**Conclusion: on x86_64, the AVX2 backend accelerates ed25519/ristretto
variable-base and multiscalar work, and does nothing at all for either X25519
operation.** A binary with AVX2 compiled in spending its time in `serial::u64`
is the expected behaviour, not a misconfiguration.

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
| `benches/harness/` | standalone crate with one set of kernels compiled for **both** x86_64 and wasm32; `src/main.rs` is the native driver, `run.mjs` the Node driver |
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
and it is the only code change in this work.

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
enables fat LTO.** No profile regresses. The pre-change code was flat at
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
multiplication by construction, and the differential tests in §8 check it.
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

| profile | base | after §4 + §5 | total |
| --- | ---: | ---: | ---: |
| x86_64, cargo `--release` default | 59 447 | **53 935** | **-9.3%** |
| x86_64, `lto=fat` | 59 841 | **46 860** | **-21.7%** |
| x86_64, `lto=fat`, `+adx,+bmi2` | 53 930 | **40 857** | **-24.2%** |
| wasm32, `serial::u32` | 149 607 | **139 100** | **-7.0%** |
| wasm32, `fiat_u32` | 139 752 | **127 942** | **-8.5%** |

wasm32 gains here where it gained nothing from §4, because `serial::u32`
already had a dedicated `square` and never had the `pow2k(1)` defect.

---

## 6. Front B — wasm32 with 64-bit limbs: **refused**

`build.rs` picks `curve25519_dalek_bits` from `target_pointer_width`, so wasm32
gets `DalekBits::Dalek32` (`serial::u32::FieldElement2625`), carrying the note:

```rust
//TODO(Wasm32): Needs tests + benchmarks to back this up
```

The hypothesis was that because wasm32 has native `i64.mul`/`i64.add`, a `u64`
is not emulated the way it would be on a real 32-bit ARM, so `bits="64"` might
win.

### 5.1 Method

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

### 5.2 Results

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

### 5.3 Verdict

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

## 7. The rest of the inventory

### 6.1 `pow2k` versus repeated `square`

Answered in §4. Summary: `pow2k` is used where it should be (only the inversion
tail needs `k > 1`), its amortization is real on this target — after the change
`pow2k(50)` is 14.09 ns/squaring against 15.05 ns for a standalone `square`, a
6% edge that justifies keeping it — and the bug was that the *common* case was
routed through it and paid 8.1 ns per squaring for the privilege.

Criterion's `pow2k(k)` ladder after the change (default flags, no forced LTO):
`pow2k(1)` 17.09 ns, `pow2k(2)` 30.64, `pow2k(4)` 60.26, `pow2k(10)` 149.33,
`pow2k(50)` 739.45 (14.79/sq), `pow2k(100)` 1478 (14.78/sq) — i.e. the
per-squaring cost flattens by about k = 10.

### 6.2 Cost of the `fiat` backend

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

### 6.3 Could the vector backend cover the Montgomery ladder?

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

---

## 8. Validation

* **Full test suite**, `--all-features`, on every backend path — with and
  without the new cfg, since a path only tested when enabled is not tested:

  | configuration | result |
  | --- | --- |
  | default (`simd` / `bits=64`) | 154 passed, 21 doctests |
  | `curve25519_dalek_backend="serial"` | 154 passed, 21 doctests |
  | `curve25519_dalek_backend="fiat"` | 148 passed, 21 doctests |
  | `curve25519_dalek_bits="32"` | 149 passed, 21 doctests |
  | `curve25519_dalek_bits="32"` + `fiat` | 148 passed, 21 doctests |
  | `-C target-feature=+adx,+bmi2` | 154 passed, 21 doctests |
  | `-C target-cpu=native` | 164 passed, 21 doctests |
  | `--cfg curve25519_dalek_bench_internals` | 154 passed, 21 doctests |
  | `--release` (debug assertions off) | 153 passed, 21 doctests |

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

  For `mul121666` (§5), in both `serial::u64` and `serial::u32`:

  * `mul121666_matches_general_mul` — limb for limb against
    `&x * &APLUS2_OVER_FOUR` on 512 random elements at each width up to the top
    of the documented bit excess, plus the all-ones and degenerate limb
    patterns. This is the specialization's correctness condition.
  * `mul121666_is_multiplication_by_121666` — double-and-add over the bits of
    121666 using only `add` and `reduce`, which share no code with either
    `mul121666` or `mul`. This is what would catch the two being wrong the same
    way.

* **RFC 7748 known-answer tests.** `montgomery.rs` changed, so the ladder's own
  vectors were run from the dependent crate: `x25519-dalek`'s
  `rfc7748_ladder_test1_vectorset1`/`vectorset2` and the four
  `rfc7748_diffie_hellman` vectors pass, and so does `rfc7748_ladder_test2` —
  the 1000-iteration iterated-ladder vector that is `#[ignore]`d for cost — run
  explicitly in release (56.8 s). `ed25519-dalek`'s suite passes too.

* **Constant time.** `square_limbs` and `mul121666` are straight-line code: no
  branch, no data-dependent memory index, no early return. Verified in the generated
  assembly — the emitted `mul` contains zero `j*`/`cmov`/`set*`; `pow2k` retains
  exactly one `jne`, which is its `k` loop, and `k` is a public exponent-chain
  constant, never secret. `square` no longer appears as an outlined symbol
  because it is a single squaring the optimizer folds into its callers. The
  refactor moves existing arithmetic without adding a conditional. After §5 the
  whole of `differential_add_and_double` still disassembles with an empty
  branch/`cmov`/`set*` set.

* **wasm32 builds**: `cargo build --target wasm32-unknown-unknown` succeeds for
  the default, for `--cfg curve25519_dalek_bits="64"`, and for
  `--no-default-features`.

* **`cargo fmt --all`** clean.

* **`cargo clippy --all-targets -- -D warnings`**: 11 errors, all pre-existing on
  the base branch (`edwards.rs` ×6, `ristretto.rs` ×3, `montgomery.rs` ×2).
  Verified identical before and after by stashing the change and re-running
  per-package. The new harness crate is clippy-clean.

---

## 9. Summary

| front | outcome |
| --- | --- |
| **A — ADX/BMI2 asm on x86_64** | **Refused.** LLVM emits all 25/25 and 15/15 available `mulx` under `+bmi2`; `adcx`/`adox` have no second carry chain to run in a 5×51 two-word accumulator, and their payoff belongs to a 4×64 saturated layout that is out of scope. Calibration: an intervention giving isolated `mul` −30% moved `mul_clamped` by 0.5%. |
| **A′ — build flags (no code)** | `-C target-feature=+adx,+bmi2` is worth −10% to −12% on `mul_clamped`. Deployment finding for consumers. |
| **B — wasm32 `bits="64"`** | **Refused.** 2.31× slower than the current default. wasm has no 64×64→128 multiply, so `u128` products are emulated. `build.rs`'s `TODO(Wasm32)` closed with evidence; behaviour unchanged. |
| **C — `square` via `pow2k(1)`** | **Changed.** `mul_clamped` −3.5% on a stock release build, −21.9% with fat LTO, −23.5% with fat LTO and `+adx,+bmi2`. Bit-for-bit identical output, no `unsafe`, no representation change, no API change, `serial::u32`/`fiat` untouched. |
| **D — ladder's multiply by 121666** | **Changed.** A specialized `mul121666` (fiat's verified `carry_scmul_121666` on the fiat backends) replaces a general multiplication whose operand had four zero limbs. Isolated 2.5x cheaper on x86_64, 3.8x on wasm32. `mul_clamped` −6.0% on a stock release build and −7.0% on wasm32; nothing under fat LTO, where the inliner already folded it. Complementary to C: between them every build profile improves. |
| **Combined C + D** | x86_64 stock release **−9.3%**, fat LTO **−21.7%**, fat LTO + `+adx,+bmi2` **−24.2%**; wasm32 **−7.0%** (`serial::u32`) and **−8.5%** (`fiat_u32`). |
| **Vector backend for Montgomery** | Viable but not worthwhile for single exchanges; the win would require a batched multi-exchange API. Not implemented. |
