// -*- mode: rust; -*-
//
// This file is part of curve25519-dalek.
// See LICENSE for licensing information.

//! Hooks that let the crate's own benchmarks reach the internal field
//! arithmetic.
//!
//! **This module is not part of the public API of `curve25519-dalek`.** It only
//! exists when the crate is built with
//!
//! ```sh
//! RUSTFLAGS='--cfg curve25519_dalek_bench_internals' cargo bench
//! ```
//!
//! and it is absent from every ordinary build, so nothing here is covered by
//! the crate's semver guarantees. It exists because `FieldElement` is
//! `pub(crate)` while Criterion benchmarks are compiled as separate crates:
//! without a hook like this the only way to measure `mul`/`square`/`pow2k` is
//! through a full scalar multiplication, which is exactly the confound the
//! benchmarks are trying to avoid.
//!
//! See `benches/x25519_field.rs` and `docs/perf-x25519-field-arithmetic.md`.

/// The backend-selected field element type, i.e. whatever
/// `curve25519_dalek_bits` and `curve25519_dalek_backend` resolved to.
///
/// Note that this is always a *serial* type: no branch of the `cfg_if!` in
/// `crate::field` resolves to a vectorized field element, which is why the AVX2
/// backend does not accelerate the Montgomery ladder.
pub type FieldElement = crate::field::FieldElement;

/// Limb width of the compiled field backend, i.e. what `curve25519_dalek_bits`
/// resolved to. Reported by the benchmark drivers so a result record states the
/// backend that was actually built rather than the one that was requested.
pub const LIMB_BITS: u32 = if cfg!(curve25519_dalek_bits = "64") {
    64
} else {
    32
};

/// Name of the compiled backend, i.e. what `curve25519_dalek_backend` resolved
/// to. Note that `"simd"` and `"avx512"` only affect Edwards/Ristretto
/// arithmetic: the field arithmetic under X25519 is serial in every case.
pub const BACKEND: &str = if cfg!(curve25519_dalek_backend = "fiat") {
    "fiat"
} else if cfg!(curve25519_dalek_backend = "avx512") {
    "avx512"
} else if cfg!(curve25519_dalek_backend = "simd") {
    "simd"
} else {
    "serial"
};

/// `FieldElement::invert`, which is `pub(crate)` because it lives in
/// `crate::field` rather than in the backend.
pub fn invert(x: &FieldElement) -> FieldElement {
    x.invert()
}

/// Number of `differential_add_and_double` steps performed by one
/// `MontgomeryPoint::mul_clamped` (255 bits, the MSB is skipped).
pub const LADDER_STEPS: usize = 255;

/// Field multiplications performed by one `differential_add_and_double`.
pub const MULS_PER_LADDER_STEP: usize = 6;

/// Field squarings performed by one `differential_add_and_double`.
pub const SQUARES_PER_LADDER_STEP: usize = 4;

/// Field multiplications performed by `FieldElement::invert` plus the final
/// `ProjectivePoint::as_affine` multiplication.
pub const MULS_PER_FINAL_INVERSION: usize = 13;

/// Field squarings performed by `FieldElement::invert`.
pub const SQUARES_PER_FINAL_INVERSION: usize = 254;

/// Total field multiplications in one `MontgomeryPoint::mul_clamped`.
pub const MULS_PER_MUL_CLAMPED: usize =
    LADDER_STEPS * MULS_PER_LADDER_STEP + MULS_PER_FINAL_INVERSION;

/// Total field squarings in one `MontgomeryPoint::mul_clamped`.
pub const SQUARES_PER_MUL_CLAMPED: usize =
    LADDER_STEPS * SQUARES_PER_LADDER_STEP + SQUARES_PER_FINAL_INVERSION;
