// -*- mode: rust; -*-
//
// This file is part of curve25519-dalek.
// Copyright (c) 2016-2021 isis lovecruft
// Copyright (c) 2016-2019 Henry de Valence
// See LICENSE for licensing information.

#![allow(non_snake_case)]

#[curve25519_dalek_derive::unsafe_target_feature_specialize(
    "avx2",
    conditional("avx512ifma,avx512vl", curve25519_dalek_backend = "avx512")
)]
pub mod spec {

    #[for_target_feature("avx2")]
    use crate::backend::vector::avx2::{CachedPoint, ExtendedPoint};

    #[for_target_feature("avx512ifma")]
    use crate::backend::vector::ifma::{CachedPoint, ExtendedPoint};

    #[cfg(feature = "precomputed-tables")]
    #[for_target_feature("avx2")]
    use crate::backend::vector::avx2::constants::BASEPOINT_ODD_LOOKUP_TABLE;

    #[for_target_feature("avx512ifma")]
    use crate::backend::vector::ifma::constants::BASEPOINT_ODD_LOOKUP_TABLE;

    use crate::backend::util::add_naf_digit;
    use crate::constants;
    use crate::edwards::EdwardsPoint;
    use crate::scalar::HEEA_MAX_INDEX;
    use crate::scalar::{HalfWidthScalar, Scalar};
    #[allow(unused_imports)]
    use crate::traits::Identity;
    use crate::window::NafLookupTable5;

    /// Compute \\(a_1 A_1 + a_2 A_2 + b B\\) in variable time, where \\(B\\) is the Ed25519 basepoint.
    ///
    /// \\(a_1\\) and \\(a_2\\) are [`HalfWidthScalar`]s, i.e. they are less than \\(2^{128}\\),
    /// while \\(b\\) is a full 256-bit scalar. That bound is what makes the strategy below sound;
    /// because it is carried by the type, there is nothing to check here and no input can make
    /// this function panic or return a wrong answer.
    ///
    /// # Optimization Strategy
    ///
    /// The function decomposes the 256-bit scalar \\(b\\) as \\(b = b_{lo} + b_{hi} \cdot 2^{128}\\),
    /// where both \\(b_{lo}\\) and \\(b_{hi}\\) are 128-bit values. This allows computing:
    ///
    /// \\[
    /// a_1 A_1 + a_2 A_2 + b_{lo} B + b_{hi} B'
    /// \\]
    ///
    /// where \\(B' = B \cdot 2^{128}\\) is a precomputed constant. Now all four scalars
    /// (\\(a_1, a_2, b_{lo}, b_{hi}\\)) are 128 bits, and two of the bases (\\(B\\) and \\(B'\\))
    /// use precomputed tables.
    ///
    /// # Implementation
    ///
    /// - For \\(A_1\\) and \\(A_2\\): NAF with window width 5 (8 precomputed points each)
    /// - For \\(B\\): NAF with window width 8 when precomputed tables available (64 points), otherwise width 5
    /// - For \\(B'\\): NAF with window width 5
    ///
    /// The algorithm shares doublings across all four scalar multiplications, processing
    /// only 128 bits instead of 256, providing approximately 2x speedup over the naive approach.
    ///
    /// This SIMD implementation uses vectorized point operations (AVX2 or AVX512-IFMA) for
    /// improved performance over the serial backend.
    pub fn mul_128_128_256(
        a1: &HalfWidthScalar,
        A1: &EdwardsPoint,
        a2: &HalfWidthScalar,
        A2: &EdwardsPoint,
        b: &Scalar,
    ) -> EdwardsPoint {
        // Decompose b into b_lo (lower 128 bits) and b_hi (upper 128 bits), so that
        // b = b_lo + b_hi * 2^128
        let (b_lo, b_hi) = b.split_at_128();

        // Compute NAF representations (all scalars are now ~128 bits)
        let a1_naf = a1.non_adjacent_form(5);
        let a2_naf = a2.non_adjacent_form(5);

        #[cfg(feature = "precomputed-tables")]
        let b_lo_naf = b_lo.non_adjacent_form(8);
        #[cfg(not(feature = "precomputed-tables"))]
        let b_lo_naf = b_lo.non_adjacent_form(5);

        let b_hi_naf = b_hi.non_adjacent_form(5);

        // Find starting index - check all NAFs up to bit 127
        // (with potential carry to bit 128 or 129)
        let mut i: usize = HEEA_MAX_INDEX;
        for j in (0..=HEEA_MAX_INDEX).rev() {
            i = j;
            if a1_naf[i] != 0 || a2_naf[i] != 0 || b_lo_naf[i] != 0 || b_hi_naf[i] != 0 {
                break;
            }
        }

        // Create lookup tables using SIMD-optimized CachedPoint
        let table_A1 = NafLookupTable5::<CachedPoint>::from(A1);
        let table_A2 = NafLookupTable5::<CachedPoint>::from(A2);

        #[cfg(feature = "precomputed-tables")]
        let table_B = &BASEPOINT_ODD_LOOKUP_TABLE;
        #[cfg(not(feature = "precomputed-tables"))]
        let table_B = &NafLookupTable5::<CachedPoint>::from(&constants::ED25519_BASEPOINT_POINT);

        // B' = B * 2^128 (precomputed constant point)
        // TODO: For optimal performance, this should also use the wider lookup table when precomputed-tables is enabled
        let table_B_128 =
            &NafLookupTable5::<CachedPoint>::from(&constants::ED25519_BASEPOINT_128_POINT);

        let mut Q = ExtendedPoint::identity();

        loop {
            Q = Q.double();

            add_naf_digit!(Q, a1_naf[i], table_A1);
            add_naf_digit!(Q, a2_naf[i], table_A2);
            add_naf_digit!(Q, b_lo_naf[i], table_B);
            // B' = B * 2^128
            add_naf_digit!(Q, b_hi_naf[i], table_B_128);

            if i == 0 {
                break;
            }
            i -= 1;
        }

        Q.into()
    }
}

#[cfg(test)]
mod test {

    use crate::constants;
    use crate::scalar::{HalfWidthScalar, Scalar};

    // Proptest for the SIMD `vartime_triple_scalar_mul_basepoint` equivalence.
    //
    // This mirrors the serial-backend proptest, but drives the vectorized
    // implementation via the runtime backend dispatcher, which safely selects
    // the AVX2/AVX512 path on SIMD-capable machines.
    proptest::proptest! {
        #[test]
        fn proptest_triple_scalar_mul_equivalence(
            a1_bytes_16 in proptest::array::uniform16(proptest::num::u8::ANY),
            a2_bytes_16 in proptest::array::uniform16(proptest::num::u8::ANY),
            b_bytes in proptest::array::uniform32(proptest::num::u8::ANY),
            A1_scalar_bytes in proptest::array::uniform32(proptest::num::u8::ANY),
            A2_scalar_bytes in proptest::array::uniform32(proptest::num::u8::ANY),
        ) {
            // Construct 128-bit scalars a1 and a2
            let a1 = HalfWidthScalar::from_bytes(a1_bytes_16);
            let a2 = HalfWidthScalar::from_bytes(a2_bytes_16);

            // Construct full 256-bit scalar b
            let b = Scalar::from_bytes_mod_order(b_bytes);

            // Generate random points A1 and A2 using scalar multiplication of basepoint
            let A1_scalar = Scalar::from_bytes_mod_order(A1_scalar_bytes);
            let A2_scalar = Scalar::from_bytes_mod_order(A2_scalar_bytes);
            let A1 = &constants::ED25519_BASEPOINT_POINT * &A1_scalar;
            let A2 = &constants::ED25519_BASEPOINT_POINT * &A2_scalar;

            // Compute using the optimized triple-base function (SIMD backend dispatch)
            let result_optimized =
                crate::backend::vartime_triple_base_mul_128_128_256(&a1, &A1, &a2, &A2, &b);

            // Compute using raw operations: a1*A1 + a2*A2 + b*B
            let expected = &(&(a1.as_scalar() * &A1) + &(a2.as_scalar() * &A2))
                + &(&b * &constants::ED25519_BASEPOINT_POINT);

            proptest::prop_assert_eq!(result_optimized, expected, "Optimized triple scalar mul should equal raw operations");
        }
    }
}
