// -*- mode: rust; -*-
//
// This file is part of curve25519-dalek.
// Copyright (c) 2016-2021 isis lovecruft
// Copyright (c) 2016-2019 Henry de Valence
// See LICENSE for licensing information.

#![allow(non_snake_case)]

use crate::backend::serial::curve_models::{ProjectiveNielsPoint, ProjectivePoint};
use crate::backend::util::add_naf_digit;
use crate::constants;
use crate::edwards::EdwardsPoint;
use crate::scalar::{HEEA_MAX_INDEX, HalfWidthScalar, Scalar};
use crate::traits::Identity;
use crate::window::NafLookupTable5;

/// Compute \\(a_1 A_1 + a_2 A_2 + b B\\) in variable time, where \\(B\\) is the Ed25519 basepoint.
///
/// \\(a_1\\) and \\(a_2\\) are [`HalfWidthScalar`]s, i.e. they are less than \\(2^{128}\\), while
/// \\(b\\) is a full 256-bit scalar. That bound is what makes the strategy below sound; because it
/// is carried by the type, there is nothing to check here and no input can make this function
/// panic or return a wrong answer.
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
/// - For \\(B\\): NAF with window width 8 when precomputed tables available (64 points)
/// - For \\(B'\\): NAF with window width 5 (could be optimized with precomputed table)
///
/// The algorithm shares doublings across all four scalar multiplications, processing
/// only 128 bits instead of 256, providing approximately 2x speedup over the naive approach.
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

    // Create lookup tables
    let table_A1 = NafLookupTable5::<ProjectiveNielsPoint>::from(A1);
    let table_A2 = NafLookupTable5::<ProjectiveNielsPoint>::from(A2);

    #[cfg(feature = "precomputed-tables")]
    let table_B = &constants::AFFINE_ODD_MULTIPLES_OF_BASEPOINT;
    #[cfg(not(feature = "precomputed-tables"))]
    let table_B =
        &NafLookupTable5::<ProjectiveNielsPoint>::from(&constants::ED25519_BASEPOINT_POINT);

    // B' = B * 2^128 (precomputed constant point)
    #[cfg(feature = "precomputed-tables")]
    let table_B_128 = &constants::AFFINE_ODD_MULTIPLES_OF_BASEPOINT_128;
    #[cfg(not(feature = "precomputed-tables"))]
    let table_B_128 =
        &NafLookupTable5::<ProjectiveNielsPoint>::from(&constants::ED25519_BASEPOINT_128_POINT);

    let mut r = ProjectivePoint::identity();

    loop {
        let mut t = r.double();

        add_naf_digit!(t.as_extended(), a1_naf[i], table_A1);
        add_naf_digit!(t.as_extended(), a2_naf[i], table_A2);
        add_naf_digit!(t.as_extended(), b_lo_naf[i], table_B);
        // B' = B * 2^128
        add_naf_digit!(t.as_extended(), b_hi_naf[i], table_B_128);

        r = t.as_projective();

        if i == 0 {
            break;
        }
        i -= 1;
    }

    r.as_extended()
}

#[cfg(test)]
mod test {

    use super::*;
    use crate::scalar::Scalar;

    fn random_scalar() -> Scalar {
        let mut wide = [0u8; 64];
        getrandom::fill(&mut wide).unwrap();
        Scalar::from_bytes_mod_order_wide(&wide)
    }

    /// Compute a1*A1 + a2*A2 + b*B with plain, unoptimized operations.
    fn naive(
        a1: &HalfWidthScalar,
        A1: &EdwardsPoint,
        a2: &HalfWidthScalar,
        A2: &EdwardsPoint,
        b: &Scalar,
    ) -> EdwardsPoint {
        &(&(a1.as_scalar() * A1) + &(a2.as_scalar() * A2))
            + &(b * &constants::ED25519_BASEPOINT_POINT)
    }

    #[test]
    fn test_triple_base_multiplication() {
        // Test vectors with random scalars
        let a1 = HalfWidthScalar::from(12345u64);
        let a2 = HalfWidthScalar::from(67890u64);
        let b = random_scalar();

        // Random points (using scalar multiplication of basepoint)
        let A1 = &constants::ED25519_BASEPOINT_POINT * &Scalar::from(2u64);
        let A2 = &constants::ED25519_BASEPOINT_POINT * &Scalar::from(3u64);

        let result = mul_128_128_256(&a1, &A1, &a2, &A2, &b);

        assert_eq!(result, naive(&a1, &A1, &a2, &A2, &b));
    }

    #[test]
    fn test_triple_base_multiplication_128() {
        // Test with full-width 128-bit scalars for a1 and a2
        let a1 = HalfWidthScalar::from_bytes([0xFF; 16]);
        let a2 = HalfWidthScalar::from_bytes([
            0x12, 0x34, 0x56, 0x78, 0x9A, 0xBC, 0xDE, 0xF0, 0x11, 0x22, 0x33, 0x44, 0x55, 0x66,
            0x77, 0x88,
        ]);

        // Full 256-bit scalar for b
        let b = random_scalar();

        // Test points
        let A1 = &constants::ED25519_BASEPOINT_POINT * &Scalar::from(5u64);
        let A2 = &constants::ED25519_BASEPOINT_POINT * &Scalar::from(7u64);

        let result_128 = mul_128_128_256(&a1, &A1, &a2, &A2, &b);

        assert_eq!(
            result_128,
            naive(&a1, &A1, &a2, &A2, &b),
            "Optimized 128-bit version failed"
        );
    }

    #[test]
    fn test_triple_base_with_zero_scalars() {
        let a1 = HalfWidthScalar::ZERO;
        let a2 = HalfWidthScalar::from(123u64);
        let b = random_scalar();

        let A1 = &constants::ED25519_BASEPOINT_POINT * &Scalar::from(2u64);
        let A2 = &constants::ED25519_BASEPOINT_POINT * &Scalar::from(3u64);

        let result = mul_128_128_256(&a1, &A1, &a2, &A2, &b);

        assert_eq!(result, naive(&a1, &A1, &a2, &A2, &b));
    }

    #[test]
    fn test_triple_base_with_identity_points() {
        let a1 = HalfWidthScalar::from(111u64);
        let a2 = HalfWidthScalar::from(222u64);
        let b = random_scalar();

        let A1 = EdwardsPoint::identity();
        let A2 = &constants::ED25519_BASEPOINT_POINT * &Scalar::from(3u64);

        let result = mul_128_128_256(&a1, &A1, &a2, &A2, &b);

        assert_eq!(result, naive(&a1, &A1, &a2, &A2, &b));
    }

    #[test]
    fn test_triple_base_consistency() {
        // Test that both functions give the same result for 128-bit inputs
        let a1 = HalfWidthScalar::from(0x123456789ABCDEFu64);
        let a2 = HalfWidthScalar::from(0xFEDCBA987654321u64);
        let b = random_scalar();

        let A1 = &constants::ED25519_BASEPOINT_POINT * &Scalar::from(11u64);
        let A2 = &constants::ED25519_BASEPOINT_POINT * &Scalar::from(13u64);

        let result_optimized = mul_128_128_256(&a1, &A1, &a2, &A2, &b);

        assert_eq!(result_optimized, naive(&a1, &A1, &a2, &A2, &b));
    }

    #[test]
    fn test_triple_base_large_scalars() {
        // Test with large scalars
        let a1 = HalfWidthScalar::from_bytes([0xFF; 16]);
        let a2 = HalfWidthScalar::from_bytes([0xAA; 16]);
        let b = random_scalar();

        let A1 = &constants::ED25519_BASEPOINT_POINT * &Scalar::from(17u64);
        let A2 = &constants::ED25519_BASEPOINT_POINT * &Scalar::from(19u64);

        let result = mul_128_128_256(&a1, &A1, &a2, &A2, &b);

        assert_eq!(result, naive(&a1, &A1, &a2, &A2, &b));
    }

    // Proptest for vartime_triple_scalar_mul_basepoint equivalence
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

            // Compute using the optimized triple-base function
            let result_optimized = mul_128_128_256(&a1, &A1, &a2, &A2, &b);

            proptest::prop_assert_eq!(
                result_optimized,
                naive(&a1, &A1, &a2, &A2, &b),
                "Optimized triple scalar mul should equal raw operations"
            );
        }
    }
}
