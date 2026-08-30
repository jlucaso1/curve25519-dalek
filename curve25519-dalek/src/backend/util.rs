// -*- mode: rust; -*-
//
// This file is part of curve25519-dalek.
// Copyright (c) 2016-2021 isis lovecruft
// Copyright (c) 2016-2019 Henry de Valence
// See LICENSE for licensing information.

//! Helpers shared by the serial and SIMD scalar multiplication implementations.

/// Accumulate the NAF digit `$naf` of one scalar into the running sum `$acc`, looking the
/// corresponding odd multiple of that scalar's base point up in `$table`.
///
/// A zero digit contributes nothing, which is what lets a caller share a single doubling of
/// `$acc` across all of its scalars.
///
/// The two backends accumulate into different curve models, so `$acc` may be followed by a
/// method call that converts it into an addable form. The serial backends hold a
/// [`CompletedPoint`], which has to go through [`as_extended()`]:
///
/// ```ignore
/// add_naf_digit!(t.as_extended(), a_naf[i], table_A);
/// ```
///
/// while the SIMD backends hold an `ExtendedPoint`, which is already addable:
///
/// ```ignore
/// add_naf_digit!(Q, a_naf[i], table_A);
/// ```
///
/// Either way the sum is assigned back to `$acc` itself.
///
/// This is a macro rather than a function for two reasons. The callers' lookup tables have
/// several different types ([`NafLookupTable5`]/[`NafLookupTable8`], over
/// [`ProjectiveNielsPoint`], [`AffineNielsPoint`], or a backend's `CachedPoint`) with no common
/// trait to be generic over. And expanding in place keeps the SIMD backends' point arithmetic
/// inside the enclosing `#[target_feature]` function, which a closure would not inherit.
///
/// [`CompletedPoint`]: crate::backend::serial::curve_models::CompletedPoint
/// [`as_extended()`]: crate::backend::serial::curve_models::CompletedPoint::as_extended
/// [`ProjectiveNielsPoint`]: crate::backend::serial::curve_models::ProjectiveNielsPoint
/// [`AffineNielsPoint`]: crate::backend::serial::curve_models::AffineNielsPoint
/// [`NafLookupTable5`]: crate::window::NafLookupTable5
/// [`NafLookupTable8`]: crate::window::NafLookupTable8
macro_rules! add_naf_digit {
    ($acc:ident $(. $as_addend:ident ())?, $naf:expr, $table:expr) => {{
        let digit = $naf;
        match digit.cmp(&0) {
            core::cmp::Ordering::Greater => {
                $acc = &$acc $(. $as_addend())? + &$table.select(digit as usize)
            }
            core::cmp::Ordering::Less => {
                $acc = &$acc $(. $as_addend())? - &$table.select(-digit as usize)
            }
            core::cmp::Ordering::Equal => {}
        }
    }};
}

pub(crate) use add_naf_digit;
