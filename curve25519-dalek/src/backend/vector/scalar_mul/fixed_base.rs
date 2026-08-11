// -*- mode: rust; -*-
//
// This file is part of curve25519-dalek.
// See LICENSE for licensing information.

//! Constant-time fixed-base scalar multiplication over the AVX2 point types.
//!
//! The serial radix-16 ladder runs entirely in `serial::u64` even when the
//! vector backend is live: before this module, `edwards_mul_base` cost an
//! identical 153 913 instructions under `curve25519_dalek_backend="serial"`
//! and `="simd"`, because the vector backend had a variable-base path and no
//! fixed-base one. This runs the same ladder — 64 additions and one
//! `mul_by_pow_2(4)`, indexed identically — over `ExtendedPoint` and
//! `CachedPoint`. See §13.12 of `docs/perf-x25519-field-arithmetic.md`.
//!
//! **AVX2 only.** `CachedPoint` has a different limb layout in the `ifma`
//! backend, so a shipped table there would have to be a second constant with
//! different contents; the `avx512` build keeps using the serial ladder.

#![allow(non_snake_case)]

#[curve25519_dalek_derive::unsafe_target_feature_specialize("avx2")]
pub mod spec {

    #[for_target_feature("avx2")]
    use crate::backend::vector::avx2::ExtendedPoint;
    #[for_target_feature("avx2")]
    use crate::backend::vector::avx2::constants::BASEPOINT_TABLE;

    use crate::edwards::EdwardsPoint;
    use crate::scalar::Scalar;
    use crate::traits::Identity;

    /// Constant-time fixed-base scalar multiplication against the Ed25519
    /// basepoint, radix 16.
    ///
    /// Indexed exactly as `impl_basepoint_table!`'s ladder with `$radix = 4`:
    /// sub-table `i / 2` holds odd multiples of \\( (2^{8})^{i/2} B \\), the odd
    /// digits are summed, the accumulator is multiplied by \\( 2^4 \\), and the
    /// even digits are summed into it.
    pub fn mul_base(scalar: &Scalar) -> EdwardsPoint {
        let a = scalar.as_radix_2w(4);

        let tables = &BASEPOINT_TABLE.0;
        let mut P = ExtendedPoint::identity();

        for i in (0..64).filter(|x| x % 2 == 1) {
            P = &P + &tables[i / 2].select(a[i]);
        }

        P = P.mul_by_pow_2(4);

        for i in (0..64).filter(|x| x % 2 == 0) {
            P = &P + &tables[i / 2].select(a[i]);
        }

        P.into()
    }
}

#[cfg(test)]
mod test {
    use super::spec_avx2::mul_base;
    use crate::backend::vector::avx2::ExtendedPoint;
    use crate::constants::ED25519_BASEPOINT_POINT;
    use crate::edwards::EdwardsPoint;
    use crate::scalar::Scalar;
    use crate::traits::Identity;

    /// The vector ladder must agree with the serial one on canonical bytes.
    /// Unlike the batch inversion in §9.3 this *is* expected to agree exactly:
    /// no inverse representatives are involved, only a different order of the
    /// same group operations, and `compress` canonicalises regardless.
    #[test]
    fn vector_fixed_base_matches_serial() {
        let mut s = Scalar::ONE;
        for _ in 0..64 {
            assert_eq!(
                mul_base(&s).compress(),
                EdwardsPoint::mul_base(&s).compress()
            );
            s += Scalar::ONE;
        }
    }

    /// The edges the interior of the scalar range cannot reach: zero, one, and
    /// the largest scalar, whose top radix-16 digit is the one
    /// `as_radix_2w(4)` has to carry into.
    #[test]
    fn vector_fixed_base_handles_edge_scalars() {
        for s in [
            Scalar::ZERO,
            Scalar::ONE,
            -Scalar::ONE,
            Scalar::from_bytes_mod_order([0xff; 32]),
        ] {
            assert_eq!(
                mul_base(&s).compress(),
                EdwardsPoint::mul_base(&s).compress()
            );
        }
    }

    /// The table constant is generated, so it is checked against the basepoint
    /// it claims to hold multiples of rather than trusted: sub-table `i` must
    /// hold \\( k (2^{8})^{i} B \\) for every `k` it can select.
    #[test]
    fn basepoint_table_constant_holds_the_right_multiples() {
        use crate::backend::vector::avx2::constants::BASEPOINT_TABLE;

        let mut base = ED25519_BASEPOINT_POINT;
        for i in 0..32 {
            for k in 1..=8i8 {
                let want = base * Scalar::from(k as u64);
                // `CachedPoint` has no direct conversion; adding it to the
                // identity is the cheapest way back to an `EdwardsPoint`.
                let got: EdwardsPoint =
                    (&ExtendedPoint::identity() + &BASEPOINT_TABLE.0[i].select(k)).into();
                assert_eq!(
                    got.compress(),
                    want.compress(),
                    "sub-table {i}, multiple {k}"
                );
            }
            base = base.mul_by_pow_2(8);
        }
    }
}
