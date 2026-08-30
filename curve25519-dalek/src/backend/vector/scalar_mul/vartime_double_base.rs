// -*- mode: rust; -*-
//
// This file is part of curve25519-dalek.
// Copyright (c) 2016-2021 isis lovecruft
// Copyright (c) 2016-2019 Henry de Valence
// See LICENSE for licensing information.
//
// Authors:
// - isis agora lovecruft <isis@patternsinthevoid.net>
// - Henry de Valence <hdevalence@hdevalence.ca>

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

    #[cfg(feature = "precomputed-tables")]
    #[for_target_feature("avx512ifma")]
    use crate::backend::vector::ifma::constants::BASEPOINT_ODD_LOOKUP_TABLE;

    use crate::backend::util::add_naf_digit;
    use crate::edwards::EdwardsPoint;
    use crate::scalar::Scalar;
    use crate::traits::Identity;
    use crate::window::NafLookupTable5;

    /// Compute \\(aA + bB\\) in variable time, where \\(B\\) is the Ed25519 basepoint.
    pub fn mul(a: &Scalar, A: &EdwardsPoint, b: &Scalar) -> EdwardsPoint {
        let a_naf = a.non_adjacent_form(5);

        #[cfg(feature = "precomputed-tables")]
        let b_naf = b.non_adjacent_form(8);
        #[cfg(not(feature = "precomputed-tables"))]
        let b_naf = b.non_adjacent_form(5);

        // Find starting index
        let mut i: usize = 255;
        for j in (0..256).rev() {
            i = j;
            if a_naf[i] != 0 || b_naf[i] != 0 {
                break;
            }
        }

        let table_A = NafLookupTable5::<CachedPoint>::from(A);

        #[cfg(feature = "precomputed-tables")]
        let table_B = &BASEPOINT_ODD_LOOKUP_TABLE;

        #[cfg(not(feature = "precomputed-tables"))]
        let table_B =
            &NafLookupTable5::<CachedPoint>::from(&crate::constants::ED25519_BASEPOINT_POINT);

        let mut Q = ExtendedPoint::identity();

        loop {
            Q = Q.double();

            add_naf_digit!(Q, a_naf[i], table_A);
            add_naf_digit!(Q, b_naf[i], table_B);

            if i == 0 {
                break;
            }
            i -= 1;
        }

        Q.into()
    }
}
