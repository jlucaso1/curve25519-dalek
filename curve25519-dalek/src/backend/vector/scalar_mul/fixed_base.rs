// -*- mode: rust; -*-
//
// This file is part of curve25519-dalek.
// See LICENSE for licensing information.

//! **Prototype**, not wired into `mul_base`. Measures what a vectorised
//! fixed-base multiply would be worth; see `docs/perf-x25519-field-arithmetic.md`.
//!
//! The serial fixed-base ladder runs entirely in `serial::u64`, even when the
//! AVX2 backend is live — `edwards_mul_base` costs the same 153 910
//! instructions under `curve25519_dalek_backend="serial"` and `="simd"`. This
//! module runs the same radix-16 ladder over the vector types instead, so the
//! difference is the headroom.
//!
//! The table is built at construction time rather than being a static
//! constant. A shipped version would need the constant — 32 sub-tables of eight
//! `CachedPoint`s, about 40 KB — which is the cost the measurement is there to
//! justify.

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

    use crate::edwards::EdwardsPoint;
    use crate::scalar::Scalar;
    use crate::traits::Identity;
    use crate::window::LookupTable;

    /// The vector counterpart of `EdwardsBasepointTable`: 32 sub-tables, the
    /// `i`th holding odd multiples of \\( (2^{8})^{i} B \\).
    pub struct VectorBasepointTable(pub [LookupTable<CachedPoint>; 32]);

    impl VectorBasepointTable {
        /// Build the table for `basepoint`. Mirrors `impl_basepoint_table!`'s
        /// `create`, with `$radix = 4`, so the two tables index alike.
        pub fn create(basepoint: &EdwardsPoint) -> VectorBasepointTable {
            let mut table = VectorBasepointTable([LookupTable::default(); 32]);
            let mut P = *basepoint;
            for i in 0..32 {
                // P = (2^8)^i * B
                table.0[i] = LookupTable::<CachedPoint>::from(&P);
                P = P.mul_by_pow_2(8);
            }
            table
        }

        /// Constant-time fixed-base scalar multiplication, radix 16.
        ///
        /// The same 64 additions and one `mul_by_pow_2(4)` as the serial
        /// ladder; only the point arithmetic and the table scan differ.
        pub fn mul_base(&self, scalar: &Scalar) -> EdwardsPoint {
            let a = scalar.as_radix_2w(4);

            let tables = &self.0;
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
}
