# Changelog

Entries are listed in reverse chronological order per undeprecated
major series.

# 5.x series

## Unreleased

* Add `EdwardsPoint::vartime_triple_scalar_mul_basepoint`, computing `a1*A1 + a2*A2 + b*B` in roughly half the doublings of the naive approach when `a1` and `a2` are less than 2^128 ([#858](https://github.com/dalek-cryptography/curve25519-dalek/pull/858))
* Add `HalfWidthScalar`, a `Scalar` that is known to be less than 2^128 ([#858](https://github.com/dalek-cryptography/curve25519-dalek/pull/858))
* Perf: Square `FieldElement51` directly rather than through `pow2k(1)` ([#922](https://github.com/dalek-cryptography/curve25519-dalek/pull/922))

### Other Changes

* Perf: `MontgomeryPoint::to_edwards` no longer performs two field
  exponentiations. It formed the Edwards `y` with a full inversion and then
  called `decompress`, whose `sqrt_ratio_i` runs the same 250-squaring chain
  again. Decompression uses `y` only as a ratio, so substituting `y = yn/yd`
  removes the division, and `EdwardsPoint` is projective, so `yn` and `yd`
  become its `Y` and `Z`: two squarings and three multiplications replace an
  inversion and a byte round trip, -46.05%. This is the XEdDSA verification
  entry point. The result is equal as a point but not limb for limb -- `Z` is
  now `u+1` rather than `1` -- and is tested against the previous
  implementation on canonical bytes across 256 inputs and both signs.
* Perf: `&scalar * ED25519_BASEPOINT_TABLE` reaches the vector fixed-base ladder
  too. Only `EdwardsPoint::mul_base` did; the table spelling -- which this
  crate's own README and docs use -- went through the macro-generated serial
  ladder, at 153,910 instructions against 84,546 for the identical result.
  `BasepointTable::mul_base` now hands over to the dispatcher when the table is
  the crate's own basepoint table, recognised by address; a table a caller built
  for another point cannot alias a `static` and takes the serial path as before.
  `RistrettoBasepointTable` wraps the same static and is caught too.
* Perf: `serial::u64` field squaring is emitted by operand scanning, as the
  field multiply already was. Written as five sums of three products, every
  partial product had to exist before any accumulator could retire, and sixteen
  of the eighteen spill stores in the emitted code were raw products. Regrouped
  by first operand: `fe_invert` -3.83%, `x25519_mul_base_clamped` -1.11%,
  `EdwardsBasepointTable::create` -1.75%. The Montgomery ladder is unchanged --
  inlined there, the squaring was already scheduled well. Output is bit-for-bit
  identical, including outside the documented bounds.
* Perf: `EdwardsPoint::mul_base` runs on the vector backend when AVX2 is live.
  The fixed-base ladder ran entirely in `serial::u64` even with the vector
  backend selected -- `edwards_mul_base` cost an identical instruction count
  under `curve25519_dalek_backend="serial"` and `="simd"` -- because the vector
  backend had a variable-base path and no fixed-base one. The same radix-16
  ladder now runs over the AVX2 point types: `mul_base` -45.1% in instructions
  and -32.9% in wall clock, `x25519_mul_base_clamped` -24.0%. This adds a
  40,960-byte basepoint table to AVX2 builds, alongside the existing 30 KB
  serial table, which non-AVX2 builds and the public `EdwardsBasepointTable`
  API still need; it is a size-for-speed trade rather than a replacement.
  Dispatch is on the same runtime CPU check every other vector entry point
  uses, so a machine without AVX2 keeps the serial ladder. The `avx512` backend
  also keeps it, its `CachedPoint` having a different limb layout. The
  generated table is verified entry by entry against the basepoint.
* Perf: the AVX2 field multiply is instantiated once per call site. LLVM shared
  one outlined body across the three `FieldElement2625x4` multiplies on the
  verification path, where it cannot specialise on the operand shapes each site
  has; an unused const generic gives each its own copy. `ed25519_verify` -3.21%
  in instructions for +784 bytes of `.text`. Output is bit-for-bit identical --
  all instantiations compute the same function -- and the change is inert, not
  wrong, if a future compiler stops sharing the body.
* Perf: `serial::u64`'s field multiply is emitted by operand scanning. Without
  BMI2, `mul` writes `rdx:rax` and the register allocator spills: 234
  instructions, 89 of them touching `%rsp`. Streaming one `a[i]` across five
  long-lived accumulators -- the same 25 products and the same carry chain --
  gives 212 instructions and 45 stack touches. Isolated `fe_mul` -3.40%,
  `x25519_mul_clamped` -1.74%, `edwards_mul_base` 165,680 -> 161,920. Output is
  bit-for-bit identical.
* Perf: `serial::u32`'s ladder subtractions skip their reduction too, matching
  `serial::u64`. This is the backend wasm32 uses, and it closes the bit-excess
  bound with only 0.167 bits -- about 12% -- of margin, against a factor of two
  for the 64-bit case. Because the failure mode of that margin is a silent
  `u32` wraparound rather than a panic, the binding constraint is asserted
  rather than argued: `debug_assert!` checks that `19 * limb` still fits a
  `u32` on every limb of every result. `x25519_mul_clamped` -4.39% on the u32
  backend. Congruent to but not limb-for-limb identical with `Sub`.
* Perf: the Montgomery ladder indexes the scalar's bytes instead of walking an
  iterator chain, removing two bounds checks from the binary.
  `x25519_mul_clamped` -1.32%. Constant time is unaffected: the operation
  sequence is identical and the index is a public loop counter.
* Fix: `optional_multiscalar_mul` switched from Straus to Pippenger at 190
  points, below where the two actually cross. Pippenger carries a large fixed
  cost -- 43 digit columns, each summing 62 buckets regardless of input size --
  so Straus is still ahead there. Measured instruction counts put the crossover
  higher; at n = 200 the old threshold cost 3.63% more than Straus would have.
* Perf: `serial::u64` field squaring now doubles its 64-bit inputs rather than
  its 128-bit output coefficients. Ten of the twenty-five partial products
  appear twice, so they are computed once and doubled; since `2*(x*y)` equals
  `(2*x)*y`, four precomputed 64-bit doublings replace five 128-bit ones, which
  is what `serial::u32` has always done. Isolated squaring improves 5.4%, and
  10.5% with `-C target-feature=+bmi2`; field inversion 5.8%; X25519
  `mul_clamped` 2.2% with `+bmi2`. Output is bit-for-bit identical, the limb
  representation and bit-excess preconditions are unchanged, and the wasm32
  module is byte-identical since `serial::u64` is not compiled there.
* Perf: the Montgomery ladder's four subtractions no longer pay for a reduction.
  The general `Sub` must accept any input at the crate-wide bit excess, so it
  offsets by `16p` and has to reduce afterward; the ladder's subtractions all
  take `mul`/`square` outputs, which are narrow enough that an offset of `2p`
  leaves the result already in range. `serial::u64` gains a `sub_unreduced` with
  its own documented, stricter precondition — `Sub`'s contract is unchanged, and
  the other backends forward to it. X25519 `mul_clamped` improves 6.1% on a
  baseline x86-64 build and 3.7% with `-C target-feature=+avx2,+bmi2`;
  `mul_base_clamped` is unchanged. The result is congruent to but not
  limb-for-limb equal to `Sub`'s, so it is checked on field equality, on the
  limb bound, and against the RFC 7748 iterated-ladder vector.
* Perf: `EdwardsBasepointTable::create` is about 4x faster. Building a
  `LookupTable<AffineNielsPoint>` converted its eight multiples to affine one at
  a time, paying a field inversion for each, so creating a basepoint table did
  256 inversions — 91% of its cost. The multiples depend on each other's value,
  not on their affine form, so the chain now runs in extended coordinates and all
  eight conversions share a single inversion via Montgomery's trick. Table
  creation improves 74.8% (3.97x); the crate's own hot paths are unchanged, since it ships
  its basepoint table as a constant. Uses a fixed-size batch, so it needs no
  `alloc`. The batch returns a different weakly-reduced representative than
  `invert` does, so the result is equal as a field element but not limb for limb;
  it is tested against the previous construction on canonical bytes.
* Perf: the constant-time window scan crosses `subtle`'s optimisation barrier
  before its accumulator is live. `LookupTable::select` was 19.8% of a keygen.
  `subtle`'s barrier is `#[inline(never)]` around a `read_volatile`, so each
  `Choice` is a real call with a mandatory memory round-trip, and the original
  loop interleaved one crossing per table entry with a live 15-limb accumulator,
  spilling it every time. Computing all the comparison masks up front, then
  accumulating, keeps the same barrier and the same number of crossings but
  makes the spills cheap: `x25519_mul_base_clamped` -2.41%, `mul_base` -2.85%
  at radix-16, -6.88% at radix-32 and -10.16% at radix-64. The Montgomery ladder
  and signature verification are unchanged, as neither uses this scan. Constant
  time is unaffected: the comparison still goes through `ct_eq`, every branch in
  the emitted code is on the public loop index, and a test checks the new path
  against the generic one for every input on all five backend configurations.
* Perf: the AVX2 field multiply is emitted column-major. `FieldElement2625x4`'s
  10x10 schoolbook was written row-major, which needs all ten `y_j`, all nine
  `y_j_19`, all ten `x_i` and all five doubled `x_i` live simultaneously —
  roughly 44 values against 16 YMM registers — and about a quarter of the
  function's 416 instructions were spilling and reloading them. Emitting one
  column of `rhs` at a time drops peak liveness to about 27 and lets the
  loop-invariant `x_i` fold into memory operands. The multiply itself improves
  6.7% (416 to 388 instructions per call). Output is bit-identical: the same 100
  partial products summed in a different order, and `u64` addition wraps, so it
  is associative and commutative even outside the documented bounds; a test
  checks all five limbs and all eight lanes against a verbatim copy of the
  previous code.
* Perf: AVX2 `ExtendedPoint::double` fuses its two signed blends. Adding `+S2`
  into lanes A,D and `-S2` into lanes B,C used two blends and two adds per limb;
  the lane sets are disjoint and the quantity is the same, so one signed blend
  carries both. 9 operations per limb become 7. The `ifma` backend already did
  it this way. The result is limb-identical and the documented bounds are
  unchanged.
* Perf: together the two AVX2 changes take an Ed25519 verification down
  **4.23%** (330 457 to 316 477 instructions, paired measurement with setup
  differenced out), and `vartime_double_scalar_mul_basepoint` down 4.64%. They
  move every AVX2 point operation, not only verification.
* Fix: the benchmark harness linked two copies of `curve25519-dalek`. Because
  `ed25519-dalek` depends on it by version rather than by path, the harness's
  `ed25519_verify` kernel was profiling the published crate from crates.io
  instead of this repository, so it was blind to every change here. Harness-only;
  no library code is affected. See `docs/perf-x25519-field-arithmetic.md`
  section 9.7.
* Perf: together, the five X25519 changes listed here take `mul_clamped` down
  **10.1%** in a stock `cargo --release` build, **28.7%** with fat LTO and
  **30.3%** with fat LTO plus `-C target-feature=+avx2,+bmi2`, and **9.5%** on
  wasm32. Measured against the base branch in one alternating session on the
  same pinned core, minimum of 15 repetitions, winning all three paired rounds
  in every profile. `mul_base_clamped` improves 1-3% on x86_64 and is unchanged
  on wasm32. See `docs/perf-x25519-field-arithmetic.md` section 12.1.
* Perf: the Montgomery ladder's multiplication by `(A+2)/4 = 121666` no longer
  goes through the general field multiplication, whose second operand had four
  zero limbs. Each backend gains a `mul121666`; the fiat backends use
  fiat-crypto's own verified `carry_scmul_121666`. X25519 `mul_clamped` improves
  by ~6% in a stock release build and ~7-8.5% on wasm32, with bit-for-bit
  identical output on the hand-written backends.
* Perf: `ProjectivePoint` now implements `ConditionallySelectable::conditional_swap`
  and `conditional_assign` by forwarding to the field element's own masked
  exchange, instead of inheriting `subtle`'s default (a struct copy plus two
  conditional assignments). The Montgomery ladder performs this swap once per
  scalar bit; the ladder driver drops from 324 to 294 instructions per
  iteration, worth ~2.2% of X25519 on wasm32 and below the noise floor on
  x86_64.
* Docs: the README now records the build settings that matter for X25519-heavy
  workloads, which are worth more than any source change measured here: fat LTO
  plus `-C target-feature=+bmi2` takes `mul_clamped` from 58.8 us to
  40.8 us on x86_64; adding `+avx2` is worth a further 6% on `mul_base_clamped`,
  since it vectorizes the constant-time window scan that key generation is
  dominated by — the two halves are complementary, `+bmi2` moving the ladder and
  `+avx2` moving key generation; and `-C target-feature=+simd128` — stable but
  off by default on `wasm32-unknown-unknown` — is worth 4.6% on `mul_clamped`
  and 15.7% on `mul_base_clamped`. The README also records that these flags
  raise the binary's CPU requirement, unlike the crate's runtime-detected vector
  backends, and that `+adx` is deliberately absent: LLVM emits zero `adcx`/`adox`
  with or without it, so it costs a CPU generation and buys nothing.
* Docs: `docs/perf-x25519-field-arithmetic.md` records a measurement pass over
  the X25519 field arithmetic on x86_64 (ADX/BMI2) and wasm32, including why an
  ADX assembly path was not added and why wasm32 keeps the 32-bit backend. This
  closes the `TODO(Wasm32)` in `build.rs`.
* Add the `x25519_field` benchmark target and a wasm32-capable benchmark harness
  under `benches/harness/`.

## 5.0.0 - 2026-07-06

### Breaking Changes

* Update edition to 2024
* Update the MSRV from 1.60 to 1.85
* Remove `group-bits` feature due to soundness issues with underlying trait ([#909](https://github.com/dalek-cryptography/curve25519-dalek/pull/909))
* Re-export `rand_core` ([#908](https://github.com/dalek-cryptography/curve25519-dalek/pull/908))
* Rename `unstable_avx512` backend to `avx512`, and no longer require nightly for it ([#913](https://github.com/dalek-cryptography/curve25519-dalek/pull/913))
* Rename `Scalar::batch_invert` -> `Scalar::invert_batch` for consistency. Also make it no-alloc. ([#789](https://github.com/dalek-cryptography/curve25519-dalek/pull/789))
* Remove deprecated functions `FieldElement::as_bytes()` and `EdwardsPoint::nonspec_map_to_curve()` ([#778](https://github.com/dalek-cryptography/curve25519-dalek/pull/778))
* Upgrade `rand_core` dependency to v0.10.0
* Upgrade `digest` and `sha2` deps

### Other Changes

* Perf: Use maximum available NAF window size in `VartimePrecomputedStraus` ([#848](https://github.com/dalek-cryptography/curve25519-dalek/pull/848))
* Perf: Skip checking 8 candidate points in `RistrettoPoint::lizard_decode` ([#882](https://github.com/dalek-cryptography/curve25519-dalek/pull/882))
* Add Lizard bytes-to-point injection for Ristretto. Gated under `lizard` feature. ([#826](https://github.com/dalek-cryptography/curve25519-dalek/pull/826))
* Add an allocating batch inversion called `Scalar::invert_batch_alloc` ([#789](https://github.com/dalek-cryptography/curve25519-dalek/pull/789))
* Add `Scalar::div_by_2` ([#805](https://github.com/dalek-cryptography/curve25519-dalek/pull/805))
* Add `EdwardsPoint::hash_to_curve` ([#786](https://github.com/dalek-cryptography/curve25519-dalek/pull/786))
* Undeprecate `Scalar::from_bits()` ([#780](https://github.com/dalek-cryptography/curve25519-dalek/pull/780))
* Use constant-time equality testing for compressed Ristretto and Edwards points, rather than autoderived equality

# 4.x series

## 4.2.0 [YANKED] - 2025-07-08

NOTE: yanked because `hash_to_curve` was improperly implemented (#785)

* Move AVX-512 backend selection logic to a separate CFG flag that requires nightly
* Add Elligator2 hashing methods `EdwardsPoint::hash_to_curve()` and `FieldElement::hash_to_field()`
* Deprecate `FieldElement::as_bytes` in favor of `FieldElement::to_bytes`
* Remove deprecated `FieldElement::as_bytes`
* Add batch conversion function `EdwardsPoint::to_montgomery_batch()`
* Make `VartimePrecomputedStraus::optional_mixed_multiscalar_mul()` and `VartimeRistrettoPrecomputation::vartime_mixed_multiscalar_mul()` accept more points than static scalars

## 4.1.3 - 2024-06-18

* Security: Fix timing leak in Scalar subtraction on u32, u64, fiat_u32, and fiat_u64 backends
* Fix assorted new warnings and lints from rustc and clippy

## 4.1.2

* Fix nightly SIMD build

## 4.1.1

* Mark `constants::BASEPOINT_ORDER` deprecated from pub API
* Add implementation for `PrimeFieldBits`, behind the `group-bits` feature flag.

## 4.1.0

* Add arbitrary integer multiplication with `MontgomeryPoint::mul_bits_be`
* Add implementations of the `ff` and `group` traits, behind the `group` feature flag
* Adapt to new types introduced in `fiat-crypto` 0.2 in `fiat` backend
* Fix `no_std` for `fiat` backend
* Mark `Scalar::clamp_integer` as `#[must_use]`
* Various documentation fixes

## 4.0.0

### Breaking changes

* Update the MSRV from 1.41 to 1.60
* Provide SemVer policy
* Make `digest` an optional feature
* Make `rand_core` an optional feature
* Remove `std` feature flag
* Remove `nightly` feature flag
* Automatic serial backend selection between `u32` and `u64` over the default `u32`
* Backend `simd` is now automatically selected over `serial` when a supported CPU is detected
* Backend override is now via cfg(curve25519_dalek_backend) over additive features
* Provide override to select `u32` or `u64` backend via cfg(curve25519_dalek_bits)
* Replace methods `Scalar::{zero, one}` with constants `Scalar::{ZERO, ONE}`
* Deprecate `EdwardsPoint::hash_from_bytes` and rename it `EdwardsPoint::nonspec_map_to_curve`
* Require including a new trait, `use curve25519_dalek::traits::BasepointTable`
  whenever using `EdwardsBasepointTable` or `RistrettoBasepointTable`
* `Scalar::from_canonical_bytes` now returns `CtOption`
* `Scalar::is_canonical` now returns `Choice`
* Remove `Scalar::from_bytes_clamped` and `Scalar::reduce`
* Deprecate and feature-gate `Scalar::from_bits` behind `legacy_compatibility`

### Other changes

* Add `EdwardsPoint::{mul_base, mul_base_clamped}`, `MontgomeryPoint::{mul_base, mul_base_clamped}`, and `BasepointTable::mul_base_clamped`
* Add `precomputed-tables` feature
* Update Maintenance Policies for SemVer
* Migrate documentation to docs.rs hosted
* Fix backend documentation generation
* Fix panic when `Ristretto::double_and_compress_batch` receives the identity point
* Remove `byteorder` dependency
* Update the `criterion` dependency to 0.4.0
* Include README.md into crate Documentation
* Update the `rand_core` dependency version and the `rand` dev-dependency
  version.
* Relax the `zeroize` dependency to `^1`
* Update the edition from 2015 to 2021

# 3.x series

## 3.2.0

* Add support for getting the identity element for the Montgomery
  form of curve25519, which is useful in certain protocols for
  checking contributory behaviour in derivation of shared secrets.

## 3.1.2

* Revert a commit which mistakenly removed support for `zeroize` traits
  for some point types, as well as elligator2 support for Edwards points.

## 3.1.1

* Fix documentation builds on nightly due to syntax changes to
  `#![cfg_attr(feature = "nightly", doc = include_str!("../README.md"))]`.

## 3.1.0

* Add support for the Elligator2 encoding for Edwards points.
* Add two optional formally-verified field arithmetic backends which
  use the Fiat Crypto project's Rust code, which is generated from
  proofs of functional correctness checked by the Coq theorem proving
  system.
* Add support for additional sizes of precomputed tables for basepoint
  scalar multiplication.
* Fix an unused import.
* Add support for using the `zeroize` traits with all point types.
  Note that points are not automatically zeroized on Drop, but that
  consumers of `curve25519-dalek` should call these methods manually
  when needed.

## 3.0.3

* Fix documentation builds on nightly due to syntax changes to
  `#![cfg_attr(feature = "nightly", doc = include_str!("../README.md"))]`.

## 3.0.2

* Multiple documentation typo fixes.
* Fixes to make using `alloc`+`no_std` possible for stable Rust.

## 3.0.1

* Update the optional `packed-simd` dependency to rely on a newer,
  maintained version of the `packed-simd-2` crate.

## 3.0.0

### Breaking changes

* Update the `digest` dependency to `0.9`.  This requires a major version
  because the `digest` traits are part of the public API, but there are
  otherwise no changes to the API.

# 2.x series

## 2.1.3

* Fix documentation builds on nightly due to syntax changes to
  `#![fg_attr(feature = "nightly", doc = include_str!("../README.md"))]`.

## 2.1.2

* Multiple documentation typo fixes.
* Fix `alloc` feature working with stable rust.

## 2.1.1

* Update the optional `packed-simd` dependency to rely on a newer,
  maintained version of the `packed-simd-2` crate.

## 2.1.0

* Make `Scalar::from_bits` a `const fn`, allowing its use in `const` contexts.

## 2.0.0

The only significant change is the data model change to the `serde` feature;
besides the `rand_core` version bump, there are no other user-visible changes.

### Breaking changes

* Fix a data modeling error in the `serde` feature pointed out by Trevor Perrin
  which caused points and scalars to be serialized with length fields rather
  than as fixed-size 32-byte arrays.  This is a breaking change, but it fixes
  compatibility with `serde-json` and ensures that the `serde-bincode` encoding
  matches the conventional encoding for X/Ed25519.
* Update `rand_core` to `0.5`, allowing use with new `rand` versions.

### Other changes

* Switch from `clear_on_drop` to `zeroize` (by Tony Arcieri).
* Require `subtle = ^2.2.1` and remove the note advising nightly Rust, which is
  no longer required as of that version of `subtle`.  See the `subtle`
  changelog for more details.
* Update `README.md` for `2.x` series.
* Remove the `build.rs` hack which loaded the entire crate into its own
  `build.rs` to generate constants, and keep the constants in the source code.

# 1.x series

## 1.2.6

* Fixes to make using alloc+no_std possible for stable Rust.

## 1.2.5

* Update the optional `packed-simd` dependency to rely on a newer,
  maintained version of the `packed-simd-2` crate.

## 1.2.4

* Specify a semver bound for `clear_on_drop` rather than an exact version,
  addressing an issue where changes to inline assembly in rustc prevented
  `clear_on_drop` from working without an update.

## 1.2.3

* Fix an issue identified by a Quarkslab audit (and Jack Grigg), where manually
  constructing unreduced `Scalar` values, as needed for X/Ed25519, and then
  performing scalar/scalar arithmetic could compute incorrect results.
* Switch to upstream Rust intrinsics for the IFMA backend now that they exist in
  Rust and don't need to be defined locally.
* Ensure that the NAF computation works correctly, even for parameters never
  used elsewhere in the codebase.
* Minor refactoring to EdwardsPoint decompression.
* Fix broken links in documentation.
* Fix compilation on nightly broken due to changes to the `#[doc(include)]` path
  root (not quite correctly done in 1.2.2).

## 1.2.2

* Fix a typo in an internal doc-comment.
* Add the "crypto" tag to crate metadata.
* Fix compilation on nightly broken due to changes to the `#[doc(include)]` path
  root.

## 1.2.1

* Fix a bug in bucket index calculations in the Pippenger multiscalar algorithm
  for very large input sizes.
* Add a more extensive randomized multiscalar multiplication consistency check
  to the test suite to prevent regressions.
* Ensure that multiscalar and NAF computations work correctly on extremal
  `Scalar` values constructed via `from_bits`.

## 1.2.0

* New multiscalar multiplication algorithm with better performance for
  large problem sizes.  The backend algorithm is selected
  transparently using the size hints of the input iterators, so no
  changes are required for client crates to start using it.
* Equality of Edwards points is now checked in projective coordinates.
* Serde can now be used with `no_std`.

## 1.1.4

* Fix typos in documentation comments.
* Remove unnecessary `Default` bound on `Scalar::from_hash`.

## 1.1.3

* Reverts the change in 1.1.0 to allow owned and borrowed RNGs, which caused a breakage due to a subtle interaction with ownership rules.  (The `RngCore` change is retained).

## 1.1.2

* Disabled KaTeX on `docs.rs` pending proper [support upstream](https://github.com/rust-lang/docs.rs/issues/302).

## 1.1.1

* Fixed an issue related to `#[cfg(rustdoc)]` which prevented documenting multiple backends.

## 1.1.0

* Adds support for precomputation for multiscalar multiplication.
* Restructures the internal source tree into `serial` and `vector` backends (no change to external API).
* Adds a new IFMA backend which sets speed records.
* The `avx2_backend` feature is now an alias for the `simd_backend` feature, which autoselects an appropriate vector backend (currently AVX2 or IFMA).
* Replaces the `rand` dependency with `rand_core`.
* Generalizes trait bounds on `RistrettoPoint::random()` and `Scalar::random()` to allow owned and borrowed RNGs and to allow `RngCore` instead of `Rng`.

## 1.0.3

* Adds `ConstantTimeEq` implementation for compressed points.

## 1.0.2

* Fixes a typo in the naming of variables in Ristretto formulas (no change to functionality).

## 1.0.1

* Depends on the stable `2.0` version of `subtle` instead of `2.0.0-pre.0`.

## 1.0.0

Initial stable release.  Yanked due to a dependency mistake (see above).
