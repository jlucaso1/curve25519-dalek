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

//! Field arithmetic modulo \\(p = 2\^{255} - 19\\), using \\(64\\)-bit
//! limbs with \\(128\\)-bit products.

use core::fmt::Debug;
use core::ops::Neg;
use core::ops::{Add, AddAssign};
use core::ops::{Mul, MulAssign};
use core::ops::{Sub, SubAssign};

use subtle::Choice;
use subtle::ConditionallySelectable;

#[cfg(feature = "zeroize")]
use zeroize::Zeroize;

/// A `FieldElement51` represents an element of the field
/// \\( \mathbb Z / (2\^{255} - 19)\\).
///
/// In the 64-bit implementation, a `FieldElement` is represented in
/// radix \\(2\^{51}\\) as five `u64`s; the coefficients are allowed to
/// grow up to \\(2\^{54}\\) between reductions modulo \\(p\\).
///
/// # Note
///
/// The `curve25519_dalek::field` module provides a type alias
/// `curve25519_dalek::field::FieldElement` to either `FieldElement51`
/// or `FieldElement2625`.
///
/// The backend-specific type `FieldElement51` should not be used
/// outside of the `curve25519_dalek::field` module.
#[derive(Copy, Clone)]
pub struct FieldElement51(pub(crate) [u64; 5]);

impl Debug for FieldElement51 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "FieldElement51({:?})", &self.0[..])
    }
}

#[cfg(feature = "zeroize")]
impl Zeroize for FieldElement51 {
    fn zeroize(&mut self) {
        self.0.zeroize();
    }
}

impl<'a> AddAssign<&'a FieldElement51> for FieldElement51 {
    fn add_assign(&mut self, _rhs: &'a FieldElement51) {
        for i in 0..5 {
            self.0[i] += _rhs.0[i];
        }
    }
}

impl<'a> Add<&'a FieldElement51> for &FieldElement51 {
    type Output = FieldElement51;
    fn add(self, _rhs: &'a FieldElement51) -> FieldElement51 {
        let mut output = *self;
        output += _rhs;
        output
    }
}

impl<'a> SubAssign<&'a FieldElement51> for FieldElement51 {
    fn sub_assign(&mut self, _rhs: &'a FieldElement51) {
        let result = (self as &FieldElement51) - _rhs;
        self.0 = result.0;
    }
}

impl<'a> Sub<&'a FieldElement51> for &FieldElement51 {
    type Output = FieldElement51;
    fn sub(self, _rhs: &'a FieldElement51) -> FieldElement51 {
        // To avoid underflow, first add a multiple of p.
        // Choose 16*p = p << 4 to be larger than 54-bit _rhs.
        //
        // If we could statically track the bitlengths of the limbs
        // of every FieldElement51, we could choose a multiple of p
        // just bigger than _rhs and avoid having to do a reduction.
        //
        // Since we don't yet have type-level integers to do this, we
        // have to add an explicit reduction call here.
        FieldElement51::reduce([
            (self.0[0] + 36028797018963664u64) - _rhs.0[0],
            (self.0[1] + 36028797018963952u64) - _rhs.0[1],
            (self.0[2] + 36028797018963952u64) - _rhs.0[2],
            (self.0[3] + 36028797018963952u64) - _rhs.0[3],
            (self.0[4] + 36028797018963952u64) - _rhs.0[4],
        ])
    }
}

impl<'a> MulAssign<&'a FieldElement51> for FieldElement51 {
    fn mul_assign(&mut self, _rhs: &'a FieldElement51) {
        let result = (self as &FieldElement51) * _rhs;
        self.0 = result.0;
    }
}

impl<'a> Mul<&'a FieldElement51> for &FieldElement51 {
    type Output = FieldElement51;

    #[rustfmt::skip] // keep alignment of c* calculations
    fn mul(self, _rhs: &'a FieldElement51) -> FieldElement51 {
        /// Helper function to multiply two 64-bit integers with 128
        /// bits of output.
        #[inline(always)]
        fn m(x: u64, y: u64) -> u128 { (x as u128) * (y as u128) }

        // Alias self, _rhs for more readable formulas
        let a: &[u64; 5] = &self.0;
        let b: &[u64; 5] = &_rhs.0;

        // Precondition: assume input limbs a[i], b[i] are bounded as
        //
        // a[i], b[i] < 2^(51 + b)
        //
        // where b is a real parameter measuring the "bit excess" of the limbs.

        // 64-bit precomputations to avoid 128-bit multiplications.
        //
        // This fits into a u64 whenever 51 + b + lg(19) < 64.
        //
        // Since 51 + b + lg(19) < 51 + 4.25 + b
        //                       = 55.25 + b,
        // this fits if b < 8.75.
        let b1_19 = b[1] * 19;
        let b2_19 = b[2] * 19;
        let b3_19 = b[3] * 19;
        let b4_19 = b[4] * 19;

        // Multiply to get 128-bit coefficients of output
        let     c0: u128 = m(a[0], b[0]) + m(a[4], b1_19) + m(a[3], b2_19) + m(a[2], b3_19) + m(a[1], b4_19);
        let mut c1: u128 = m(a[1], b[0]) + m(a[0],  b[1]) + m(a[4], b2_19) + m(a[3], b3_19) + m(a[2], b4_19);
        let mut c2: u128 = m(a[2], b[0]) + m(a[1],  b[1]) + m(a[0],  b[2]) + m(a[4], b3_19) + m(a[3], b4_19);
        let mut c3: u128 = m(a[3], b[0]) + m(a[2],  b[1]) + m(a[1],  b[2]) + m(a[0],  b[3]) + m(a[4], b4_19);
        let mut c4: u128 = m(a[4], b[0]) + m(a[3],  b[1]) + m(a[2],  b[2]) + m(a[1],  b[3]) + m(a[0] , b[4]);

        // How big are the c[i]? We have
        //
        //    c[i] < 2^(102 + 2*b) * (1+i + (4-i)*19)
        //         < 2^(102 + lg(1 + 4*19) + 2*b)
        //         < 2^(108.27 + 2*b)
        //
        // The carry (c[i] >> 51) fits into a u64 when
        //    108.27 + 2*b - 51 < 64
        //    2*b < 6.73
        //    b < 3.365.
        //
        // So we require b < 3 to ensure this fits.
        debug_assert!(a[0] < (1 << 54)); debug_assert!(b[0] < (1 << 54));
        debug_assert!(a[1] < (1 << 54)); debug_assert!(b[1] < (1 << 54));
        debug_assert!(a[2] < (1 << 54)); debug_assert!(b[2] < (1 << 54));
        debug_assert!(a[3] < (1 << 54)); debug_assert!(b[3] < (1 << 54));
        debug_assert!(a[4] < (1 << 54)); debug_assert!(b[4] < (1 << 54));

        // Casting to u64 and back tells the compiler that the carry is
        // bounded by 2^64, so that the addition is a u128 + u64 rather
        // than u128 + u128.

        const LOW_51_BIT_MASK: u64 = (1u64 << 51) - 1;
        let mut out = [0u64; 5];

        c1 += ((c0 >> 51) as u64) as u128;
        out[0] = (c0 as u64) & LOW_51_BIT_MASK;

        c2 += ((c1 >> 51) as u64) as u128;
        out[1] = (c1 as u64) & LOW_51_BIT_MASK;

        c3 += ((c2 >> 51) as u64) as u128;
        out[2] = (c2 as u64) & LOW_51_BIT_MASK;

        c4 += ((c3 >> 51) as u64) as u128;
        out[3] = (c3 as u64) & LOW_51_BIT_MASK;

        let carry: u64 = (c4 >> 51) as u64;
        out[4] = (c4 as u64) & LOW_51_BIT_MASK;

        // To see that this does not overflow, we need out[0] + carry * 19 < 2^64.
        //
        // c4 < a0*b4 + a1*b3 + a2*b2 + a3*b1 + a4*b0 + (carry from c3)
        //    < 5*(2^(51 + b) * 2^(51 + b)) + (carry from c3)
        //    < 2^(102 + 2*b + lg(5)) + 2^64.
        //
        // When b < 3 we get
        //
        // c4 < 2^110.33  so that carry < 2^59.33
        //
        // so that
        //
        // out[0] + carry * 19 < 2^51 + 19 * 2^59.33 < 2^63.58
        //
        // and there is no overflow.
        out[0] += carry * 19;

        // Now out[1] < 2^51 + 2^(64 -51) = 2^51 + 2^13 < 2^(51 + epsilon).
        out[1] += out[0] >> 51;
        out[0] &= LOW_51_BIT_MASK;

        // Now out[i] < 2^(51 + epsilon) for all i.
        FieldElement51(out)
    }
}

impl Neg for &FieldElement51 {
    type Output = FieldElement51;
    fn neg(self) -> FieldElement51 {
        let mut output = *self;
        output.negate();
        output
    }
}

impl ConditionallySelectable for FieldElement51 {
    fn conditional_select(
        a: &FieldElement51,
        b: &FieldElement51,
        choice: Choice,
    ) -> FieldElement51 {
        FieldElement51([
            u64::conditional_select(&a.0[0], &b.0[0], choice),
            u64::conditional_select(&a.0[1], &b.0[1], choice),
            u64::conditional_select(&a.0[2], &b.0[2], choice),
            u64::conditional_select(&a.0[3], &b.0[3], choice),
            u64::conditional_select(&a.0[4], &b.0[4], choice),
        ])
    }

    fn conditional_swap(a: &mut FieldElement51, b: &mut FieldElement51, choice: Choice) {
        u64::conditional_swap(&mut a.0[0], &mut b.0[0], choice);
        u64::conditional_swap(&mut a.0[1], &mut b.0[1], choice);
        u64::conditional_swap(&mut a.0[2], &mut b.0[2], choice);
        u64::conditional_swap(&mut a.0[3], &mut b.0[3], choice);
        u64::conditional_swap(&mut a.0[4], &mut b.0[4], choice);
    }

    fn conditional_assign(&mut self, other: &FieldElement51, choice: Choice) {
        self.0[0].conditional_assign(&other.0[0], choice);
        self.0[1].conditional_assign(&other.0[1], choice);
        self.0[2].conditional_assign(&other.0[2], choice);
        self.0[3].conditional_assign(&other.0[3], choice);
        self.0[4].conditional_assign(&other.0[4], choice);
    }
}

impl FieldElement51 {
    pub(crate) const fn from_limbs(limbs: [u64; 5]) -> FieldElement51 {
        FieldElement51(limbs)
    }

    /// The scalar \\( 0 \\).
    pub const ZERO: FieldElement51 = FieldElement51::from_limbs([0, 0, 0, 0, 0]);
    /// The scalar \\( 1 \\).
    pub const ONE: FieldElement51 = FieldElement51::from_limbs([1, 0, 0, 0, 0]);
    /// The scalar \\( -1 \\).
    pub const MINUS_ONE: FieldElement51 = FieldElement51::from_limbs([
        2251799813685228,
        2251799813685247,
        2251799813685247,
        2251799813685247,
        2251799813685247,
    ]);

    /// Invert the sign of this field element
    pub fn negate(&mut self) {
        // See commentary in the Sub impl
        let neg = FieldElement51::reduce([
            36028797018963664u64 - self.0[0],
            36028797018963952u64 - self.0[1],
            36028797018963952u64 - self.0[2],
            36028797018963952u64 - self.0[3],
            36028797018963952u64 - self.0[4],
        ]);
        self.0 = neg.0;
    }

    /// Given 64-bit input limbs, reduce to enforce the bound 2^(51 + epsilon).
    #[inline(always)]
    fn reduce(mut limbs: [u64; 5]) -> FieldElement51 {
        const LOW_51_BIT_MASK: u64 = (1u64 << 51) - 1;

        // Since the input limbs are bounded by 2^64, the biggest
        // carry-out is bounded by 2^13.
        //
        // The biggest carry-in is c4 * 19, resulting in
        //
        // 2^51 + 19*2^13 < 2^51.0000000001
        //
        // Because we don't need to canonicalize, only to reduce the
        // limb sizes, it's OK to do a "weak reduction", where we
        // compute the carry-outs in parallel.

        let c0 = limbs[0] >> 51;
        let c1 = limbs[1] >> 51;
        let c2 = limbs[2] >> 51;
        let c3 = limbs[3] >> 51;
        let c4 = limbs[4] >> 51;

        limbs[0] &= LOW_51_BIT_MASK;
        limbs[1] &= LOW_51_BIT_MASK;
        limbs[2] &= LOW_51_BIT_MASK;
        limbs[3] &= LOW_51_BIT_MASK;
        limbs[4] &= LOW_51_BIT_MASK;

        limbs[0] += c4 * 19;
        limbs[1] += c0;
        limbs[2] += c1;
        limbs[3] += c2;
        limbs[4] += c3;

        FieldElement51(limbs)
    }

    /// Load a `FieldElement51` from the low 255 bits of a 256-bit
    /// input.
    ///
    /// # Warning
    ///
    /// This function does not check that the input used the canonical
    /// representative.  It masks the high bit, but it will happily
    /// decode 2^255 - 18 to 1.  Applications that require a canonical
    /// encoding of every field element should decode, re-encode to
    /// the canonical encoding, and check that the input was
    /// canonical.
    ///
    #[rustfmt::skip] // keep alignment of bit shifts
    pub const fn from_bytes(bytes: &[u8; 32]) -> FieldElement51 {
        const fn load8_at(input: &[u8], i: usize) -> u64 {
               (input[i] as u64)
            | ((input[i + 1] as u64) << 8)
            | ((input[i + 2] as u64) << 16)
            | ((input[i + 3] as u64) << 24)
            | ((input[i + 4] as u64) << 32)
            | ((input[i + 5] as u64) << 40)
            | ((input[i + 6] as u64) << 48)
            | ((input[i + 7] as u64) << 56)
        }

        let low_51_bit_mask = (1u64 << 51) - 1;
        FieldElement51(
        // load bits [  0, 64), no shift
        [  load8_at(bytes,  0)        & low_51_bit_mask
        // load bits [ 48,112), shift to [ 51,112)
        , (load8_at(bytes,  6) >>  3) & low_51_bit_mask
        // load bits [ 96,160), shift to [102,160)
        , (load8_at(bytes, 12) >>  6) & low_51_bit_mask
        // load bits [152,216), shift to [153,216)
        , (load8_at(bytes, 19) >>  1) & low_51_bit_mask
        // load bits [192,256), shift to [204,112)
        , (load8_at(bytes, 24) >> 12) & low_51_bit_mask
        ])
    }

    /// Serialize this `FieldElement51` to a 32-byte array.  The
    /// encoding is canonical.
    #[rustfmt::skip] // keep alignment of s[*] calculations
    pub fn to_bytes(self) -> [u8; 32] {
        // Let h = limbs[0] + limbs[1]*2^51 + ... + limbs[4]*2^204.
        //
        // Write h = pq + r with 0 <= r < p.
        //
        // We want to compute r = h mod p.
        //
        // If h < 2*p = 2^256 - 38,
        // then q = 0 or 1,
        //
        // with q = 0 when h < p
        //  and q = 1 when h >= p.
        //
        // Notice that h >= p <==> h + 19 >= p + 19 <==> h + 19 >= 2^255.
        // Therefore q can be computed as the carry bit of h + 19.

        // First, reduce the limbs to ensure h < 2*p.
        let mut limbs = FieldElement51::reduce(self.0).0;

        let mut q = (limbs[0] + 19) >> 51;
        q = (limbs[1] + q) >> 51;
        q = (limbs[2] + q) >> 51;
        q = (limbs[3] + q) >> 51;
        q = (limbs[4] + q) >> 51;

        // Now we can compute r as r = h - pq = r - (2^255-19)q = r + 19q - 2^255q

        limbs[0] += 19 * q;

        // Now carry the result to compute r + 19q ...
        let low_51_bit_mask = (1u64 << 51) - 1;
        limbs[1] += limbs[0] >> 51;
        limbs[0] &= low_51_bit_mask;
        limbs[2] += limbs[1] >> 51;
        limbs[1] &= low_51_bit_mask;
        limbs[3] += limbs[2] >> 51;
        limbs[2] &= low_51_bit_mask;
        limbs[4] += limbs[3] >> 51;
        limbs[3] &= low_51_bit_mask;
        // ... but instead of carrying (limbs[4] >> 51) = 2^255q
        // into another limb, discard it, subtracting the value
        limbs[4] &= low_51_bit_mask;

        // Now arrange the bits of the limbs.
        let mut s = [0u8;32];
        s[ 0] =   limbs[0]                           as u8;
        s[ 1] =  (limbs[0] >>  8)                    as u8;
        s[ 2] =  (limbs[0] >> 16)                    as u8;
        s[ 3] =  (limbs[0] >> 24)                    as u8;
        s[ 4] =  (limbs[0] >> 32)                    as u8;
        s[ 5] =  (limbs[0] >> 40)                    as u8;
        s[ 6] = ((limbs[0] >> 48) | (limbs[1] << 3)) as u8;
        s[ 7] =  (limbs[1] >>  5)                    as u8;
        s[ 8] =  (limbs[1] >> 13)                    as u8;
        s[ 9] =  (limbs[1] >> 21)                    as u8;
        s[10] =  (limbs[1] >> 29)                    as u8;
        s[11] =  (limbs[1] >> 37)                    as u8;
        s[12] = ((limbs[1] >> 45) | (limbs[2] << 6)) as u8;
        s[13] =  (limbs[2] >>  2)                    as u8;
        s[14] =  (limbs[2] >> 10)                    as u8;
        s[15] =  (limbs[2] >> 18)                    as u8;
        s[16] =  (limbs[2] >> 26)                    as u8;
        s[17] =  (limbs[2] >> 34)                    as u8;
        s[18] =  (limbs[2] >> 42)                    as u8;
        s[19] = ((limbs[2] >> 50) | (limbs[3] << 1)) as u8;
        s[20] =  (limbs[3] >>  7)                    as u8;
        s[21] =  (limbs[3] >> 15)                    as u8;
        s[22] =  (limbs[3] >> 23)                    as u8;
        s[23] =  (limbs[3] >> 31)                    as u8;
        s[24] =  (limbs[3] >> 39)                    as u8;
        s[25] = ((limbs[3] >> 47) | (limbs[4] << 4)) as u8;
        s[26] =  (limbs[4] >>  4)                    as u8;
        s[27] =  (limbs[4] >> 12)                    as u8;
        s[28] =  (limbs[4] >> 20)                    as u8;
        s[29] =  (limbs[4] >> 28)                    as u8;
        s[30] =  (limbs[4] >> 36)                    as u8;
        s[31] =  (limbs[4] >> 44)                    as u8;

        // High bit should be zero.
        debug_assert!((s[31] & 0b1000_0000u8) == 0u8);

        s
    }

    /// Given `k > 0`, return `self^(2^k)`.
    pub fn pow2k(&self, mut k: u32) -> FieldElement51 {
        debug_assert!(k > 0);

        let mut a: [u64; 5] = self.0;

        loop {
            a = square_limbs(a);

            k -= 1;
            if k == 0 {
                break;
            }
        }

        FieldElement51(a)
    }

    /// Returns the square of this field element.
    pub fn square(&self) -> FieldElement51 {
        FieldElement51(square_limbs(self.0))
    }

    /// Compute `self - rhs` for operands narrow enough not to need a reduction.
    ///
    /// The general [`Sub`] must accept any input satisfying the crate-wide bit
    /// excess `b < 3`, so it adds `16p` to keep the difference positive — which
    /// leaves limbs around `2^55` and forces a [`FieldElement51::reduce`] to get
    /// back inside the contract. That reduction is about sixteen instructions,
    /// and it is pure overhead whenever the operands are already narrow.
    ///
    /// This subtracts against `2p` instead, which is enough whenever both
    /// operands are outputs of `mul`, `square` or `square2` — all of which
    /// produce limbs below `2^51 + 2^13` — and so needs no reduction at all.
    ///
    /// # Preconditions
    ///
    /// Every limb of **both** operands must be `< 2^52 - 38`, which is the
    /// smallest limb of `2p`. `mul`/`square` outputs satisfy this with eleven
    /// bits to spare, and so do the ladder's *initial* values, which are not
    /// products: `ProjectivePoint::identity()` is `(1, 0)` and the other point
    /// is `(from_bytes(u), 1)`, all limbs `< 2^51`. The `debug_assert!`s below
    /// enforce it in debug builds, so the whole test suite checks it.
    ///
    /// # Postcondition
    ///
    /// Limbs are `< 2^53`, so the result satisfies the documented `b < 3` and
    /// may be fed to any operation in this module.
    ///
    /// The result is congruent to, but generally *not* limb-for-limb equal to,
    /// `self - rhs`: it is a different representative of the same field
    /// element. `sub_unreduced_agrees_with_sub` checks the field equality and
    /// the limb bound; the ladder's callers are checked end to end against the
    /// byte output of `mul_clamped`, which is canonical.
    pub(crate) fn sub_unreduced(&self, rhs: &FieldElement51) -> FieldElement51 {
        /// `2p`, limb by limb: `2 * (2^51 - 19)` then four times `2 * (2^51 - 1)`.
        const TWO_P: [u64; 5] = [
            4503599627370458,
            4503599627370494,
            4503599627370494,
            4503599627370494,
            4503599627370494,
        ];

        // The bound that makes the reduction unnecessary. It is stricter than
        // the module-wide `b < 3`, which is why this is a separate operation
        // rather than a change to `Sub`.
        debug_assert!(self.0.iter().all(|&l| l < TWO_P[0]));
        debug_assert!(rhs.0.iter().all(|&l| l < TWO_P[0]));

        FieldElement51([
            (self.0[0] + TWO_P[0]) - rhs.0[0],
            (self.0[1] + TWO_P[1]) - rhs.0[1],
            (self.0[2] + TWO_P[2]) - rhs.0[2],
            (self.0[3] + TWO_P[3]) - rhs.0[3],
            (self.0[4] + TWO_P[4]) - rhs.0[4],
        ])
    }

    /// Multiply this field element by \\((A+2)/4 = 121666\\), the constant the
    /// Montgomery ladder needs once per step.
    ///
    /// `&x * &constants::APLUS2_OVER_FOUR` gives the same answer, but that
    /// constant is `[121666, 0, 0, 0, 0]`, so twenty of the twenty-five partial
    /// products in the general multiplication are multiplications by zero, and
    /// the four `b[i] * 19` precomputations are `0 * 19`. The compiler cannot
    /// fold them away because `mul` is not inlined into the ladder.
    ///
    /// Only the five surviving products are computed here. The carry chain is
    /// the same one `mul` uses, so the result is bit-for-bit identical to the
    /// general multiplication; `mul121666_matches_general_mul` checks that.
    #[rustfmt::skip] // keep alignment of c* calculations
    pub fn mul121666(&self) -> FieldElement51 {
        /// \\((A+2)/4\\), the only value this is ever called with.
        const APLUS2_OVER_FOUR: u128 = 121666;

        let a: &[u64; 5] = &self.0;

        // Precondition, as for `mul`: a[i] < 2^(51 + b) with b < 3.
        debug_assert!(a[0] < (1 << 54));
        debug_assert!(a[1] < (1 << 54));
        debug_assert!(a[2] < (1 << 54));
        debug_assert!(a[3] < (1 << 54));
        debug_assert!(a[4] < (1 << 54));

        // c[i] = a[i] * 121666 < 2^54 * 2^16.9 = 2^70.9, so the carries
        // c[i] >> 51 are below 2^20 and everything below stays far inside its
        // type. This is the same shape as `mul`'s coefficients, just smaller.
        let     c0: u128 = (a[0] as u128) * APLUS2_OVER_FOUR;
        let mut c1: u128 = (a[1] as u128) * APLUS2_OVER_FOUR;
        let mut c2: u128 = (a[2] as u128) * APLUS2_OVER_FOUR;
        let mut c3: u128 = (a[3] as u128) * APLUS2_OVER_FOUR;
        let mut c4: u128 = (a[4] as u128) * APLUS2_OVER_FOUR;

        const LOW_51_BIT_MASK: u64 = (1u64 << 51) - 1;
        let mut out = [0u64; 5];

        c1 += ((c0 >> 51) as u64) as u128;
        out[0] = (c0 as u64) & LOW_51_BIT_MASK;

        c2 += ((c1 >> 51) as u64) as u128;
        out[1] = (c1 as u64) & LOW_51_BIT_MASK;

        c3 += ((c2 >> 51) as u64) as u128;
        out[2] = (c2 as u64) & LOW_51_BIT_MASK;

        c4 += ((c3 >> 51) as u64) as u128;
        out[3] = (c3 as u64) & LOW_51_BIT_MASK;

        let carry: u64 = (c4 >> 51) as u64;
        out[4] = (c4 as u64) & LOW_51_BIT_MASK;

        // carry < 2^20, so out[0] + carry * 19 < 2^51 + 2^24.3, no overflow.
        out[0] += carry * 19;

        out[1] += out[0] >> 51;
        out[0] &= LOW_51_BIT_MASK;

        FieldElement51(out)
    }

    /// Returns 2 times the square of this field element.
    pub fn square2(&self) -> FieldElement51 {
        let mut square = square_limbs(self.0);
        for limb in &mut square {
            *limb *= 2;
        }

        FieldElement51(square)
    }
}

/// One squaring step in radix \\(2^{51}\\): given limbs `a` of \\(x\\), return
/// the limbs of \\(x^2\\).
///
/// This is the body of [`FieldElement51::pow2k`], factored out so that
/// [`FieldElement51::square`] can be a single squaring rather than a call into
/// a loop. Squaring is a third of the field operations on the X25519 ladder,
/// and it was previously spelled `pow2k(1)`, which paid for `pow2k`'s call and
/// loop machinery without getting any amortization back for it.
///
/// The arithmetic below is unchanged from the original `pow2k` loop body, so
/// `pow2k(k)` and `k` chained `square()`s remain bit-for-bit identical.
#[rustfmt::skip] // keep alignment of c* calculations
#[inline(always)]
fn square_limbs(mut a: [u64; 5]) -> [u64; 5] {
    /// Multiply two 64-bit integers with 128 bits of output.
    #[inline(always)]
    fn m(x: u64, y: u64) -> u128 {
        (x as u128) * (y as u128)
    }

    // Precondition: assume input limbs a[i] are bounded as
    //
    // a[i] < 2^(51 + b)
    //
    // where b is a real parameter measuring the "bit excess" of the limbs.

    // Precomputation: 64-bit multiply by 19.
    //
    // This fits into a u64 whenever 51 + b + lg(19) < 64.
    //
    // Since 51 + b + lg(19) < 51 + 4.25 + b
    //                       = 55.25 + b,
    // this fits if b < 8.75.
    let a3_19 = 19 * a[3];
    let a4_19 = 19 * a[4];

    // Precomputation: the doublings.
    //
    // Every product below other than the five squares a[i]^2 appears exactly
    // twice in the full 5x5 product, so it is computed once and doubled. The
    // doubling is done here, on one 64-bit operand, rather than on the 128-bit
    // product: `2*(x*y) == (2*x)*y`, so the coefficients are unchanged, but a
    // 64-bit shift is one instruction where a 128-bit one is two, and four
    // shifts here cover ten doubled products.
    //
    // `2*a[i]` fits into a u64 whenever 51 + b + 1 < 64, i.e. b < 12, and
    // `2*19*a[3]` whenever 51 + b + lg(19) + 1 < 64, i.e. b < 7.75. Both are
    // slacker than the b < 3 this function already requires.
    let a0_2 = 2 * a[0];
    let a1_2 = 2 * a[1];
    let a2_2 = 2 * a[2];
    let a3_19_2 = 2 * a3_19;

    // Multiply to get 128-bit coefficients of output.
    let     c0: u128 = m(a[0], a[0]) + m(a1_2, a4_19) + m(a2_2, a3_19);
    let mut c1: u128 = m(a[3], a3_19) + m(a0_2,  a[1]) + m(a2_2, a4_19);
    let mut c2: u128 = m(a[1], a[1])  + m(a0_2,  a[2]) + m(a[4], a3_19_2);
    let mut c3: u128 = m(a[4], a4_19) + m(a0_2,  a[3]) + m(a1_2, a[2]);
    let mut c4: u128 = m(a[2], a[2])  + m(a0_2,  a[4]) + m(a1_2, a[3]);

    // Same bound as in multiply:
    //    c[i] < 2^(102 + 2*b) * (1+i + (4-i)*19)
    //         < 2^(102 + lg(1 + 4*19) + 2*b)
    //         < 2^(108.27 + 2*b)
    //
    // The carry (c[i] >> 51) fits into a u64 when
    //    108.27 + 2*b - 51 < 64
    //    2*b < 6.73
    //    b < 3.365.
    //
    // So we require b < 3 to ensure this fits.
    debug_assert!(a[0] < (1 << 54));
    debug_assert!(a[1] < (1 << 54));
    debug_assert!(a[2] < (1 << 54));
    debug_assert!(a[3] < (1 << 54));
    debug_assert!(a[4] < (1 << 54));

    const LOW_51_BIT_MASK: u64 = (1u64 << 51) - 1;

    // Casting to u64 and back tells the compiler that the carry is bounded by 2^64, so
    // that the addition is a u128 + u64 rather than u128 + u128.
    c1 += ((c0 >> 51) as u64) as u128;
    a[0] = (c0 as u64) & LOW_51_BIT_MASK;

    c2 += ((c1 >> 51) as u64) as u128;
    a[1] = (c1 as u64) & LOW_51_BIT_MASK;

    c3 += ((c2 >> 51) as u64) as u128;
    a[2] = (c2 as u64) & LOW_51_BIT_MASK;

    c4 += ((c3 >> 51) as u64) as u128;
    a[3] = (c3 as u64) & LOW_51_BIT_MASK;

    let carry: u64 = (c4 >> 51) as u64;
    a[4] = (c4 as u64) & LOW_51_BIT_MASK;

    // To see that this does not overflow, we need a[0] + carry * 19 < 2^64.
    //
    // c4 < a2^2 + 2*a0*a4 + 2*a1*a3 + (carry from c3)
    //    < 2^(102 + 2*b + lg(5)) + 2^64.
    //
    // When b < 3 we get
    //
    // c4 < 2^110.33  so that carry < 2^59.33
    //
    // so that
    //
    // a[0] + carry * 19 < 2^51 + 19 * 2^59.33 < 2^63.58
    //
    // and there is no overflow.
    a[0] += carry * 19;

    // Now a[1] < 2^51 + 2^(64 -51) = 2^51 + 2^13 < 2^(51 + epsilon).
    a[1] += a[0] >> 51;
    a[0] &= LOW_51_BIT_MASK;

    // Now all a[i] < 2^(51 + epsilon) and a = (input)^2.

    a
}

#[cfg(test)]
mod test {
    use super::*;

    /// Verbatim copy of `pow2k` as it was implemented before `square_limbs` was
    /// factored out of it, kept here as the reference for the differential
    /// tests below. If a future change to the squaring code is not a pure
    /// refactor, these tests are what will say so.
    #[rustfmt::skip]
    fn reference_pow2k(fe: &FieldElement51, mut k: u32) -> FieldElement51 {
        debug_assert!( k > 0 );

        #[inline(always)]
        fn m(x: u64, y: u64) -> u128 {
            (x as u128) * (y as u128)
        }

        let mut a: [u64; 5] = fe.0;

        loop {
            let a3_19 = 19 * a[3];
            let a4_19 = 19 * a[4];

            let     c0: u128 = m(a[0],  a[0]) + 2*( m(a[1], a4_19) + m(a[2], a3_19) );
            let mut c1: u128 = m(a[3], a3_19) + 2*( m(a[0],  a[1]) + m(a[2], a4_19) );
            let mut c2: u128 = m(a[1],  a[1]) + 2*( m(a[0],  a[2]) + m(a[4], a3_19) );
            let mut c3: u128 = m(a[4], a4_19) + 2*( m(a[0],  a[3]) + m(a[1],  a[2]) );
            let mut c4: u128 = m(a[2],  a[2]) + 2*( m(a[0],  a[4]) + m(a[1],  a[3]) );

            const LOW_51_BIT_MASK: u64 = (1u64 << 51) - 1;

            c1 += ((c0 >> 51) as u64) as u128;
            a[0] = (c0 as u64) & LOW_51_BIT_MASK;

            c2 += ((c1 >> 51) as u64) as u128;
            a[1] = (c1 as u64) & LOW_51_BIT_MASK;

            c3 += ((c2 >> 51) as u64) as u128;
            a[2] = (c2 as u64) & LOW_51_BIT_MASK;

            c4 += ((c3 >> 51) as u64) as u128;
            a[3] = (c3 as u64) & LOW_51_BIT_MASK;

            let carry: u64 = (c4 >> 51) as u64;
            a[4] = (c4 as u64) & LOW_51_BIT_MASK;

            a[0] += carry * 19;

            a[1] += a[0] >> 51;
            a[0] &= LOW_51_BIT_MASK;

            k -= 1;
            if k == 0 {
                break;
            }
        }

        FieldElement51(a)
    }

    /// Pre-refactor `square2`, for the same reason.
    fn reference_square2(fe: &FieldElement51) -> FieldElement51 {
        let mut square = reference_pow2k(fe, 1);
        for i in 0..5 {
            square.0[i] *= 2;
        }
        square
    }

    /// Deterministic xorshift64, so the differential tests are reproducible and
    /// need no dependencies.
    struct Rng(u64);

    impl Rng {
        fn next(&mut self) -> u64 {
            self.0 ^= self.0 << 13;
            self.0 ^= self.0 >> 7;
            self.0 ^= self.0 << 17;
            self.0
        }

        /// A field element whose limbs are uniform below `2^bits`.
        ///
        /// The documented precondition of `mul`/`pow2k` is a bit excess
        /// `b < 3`, i.e. limbs below `2^54`, so `bits == 54` is exactly the
        /// upper edge of the supported input range.
        fn field_element(&mut self, bits: u32) -> FieldElement51 {
            let mask = (1u64 << bits) - 1;
            FieldElement51([
                self.next() & mask,
                self.next() & mask,
                self.next() & mask,
                self.next() & mask,
                self.next() & mask,
            ])
        }
    }

    /// The refactored `square` must agree with the old `pow2k(1)` limb for
    /// limb, not merely as a field element: the limbs are an input to the next
    /// operation's bit-excess accounting.
    #[test]
    fn square_matches_reference_bit_for_bit() {
        let mut rng = Rng(0x1234_5678_9abc_def1);

        // 51 bits is the reduced case; 54 bits is the documented upper edge of
        // the bit excess.
        for bits in [51u32, 52, 53, 54] {
            for _ in 0..512 {
                let x = rng.field_element(bits);
                assert_eq!(x.square().0, reference_pow2k(&x, 1).0);
                assert_eq!(x.square2().0, reference_square2(&x).0);
            }
        }

        // The all-ones limbs at the top of the documented range.
        let edge = FieldElement51([(1u64 << 54) - 1; 5]);
        assert_eq!(edge.square().0, reference_pow2k(&edge, 1).0);
        assert_eq!(edge.square2().0, reference_square2(&edge).0);

        // And the degenerate inputs.
        for limbs in [[0u64; 5], [1, 0, 0, 0, 0], [0, 0, 0, 0, 1]] {
            let x = FieldElement51(limbs);
            assert_eq!(x.square().0, reference_pow2k(&x, 1).0);
            assert_eq!(x.square2().0, reference_square2(&x).0);
        }
    }

    /// `pow2k(k)` is now a loop over the same factored-out step, so it must
    /// still agree with the old monolithic loop for every `k`.
    #[test]
    fn pow2k_matches_reference_bit_for_bit() {
        let mut rng = Rng(0xfeed_face_dead_beef);

        for bits in [51u32, 54] {
            for _ in 0..128 {
                let x = rng.field_element(bits);
                for k in 1..=8u32 {
                    assert_eq!(x.pow2k(k).0, reference_pow2k(&x, k).0);
                }
                // The k values actually used on the X25519 inversion path.
                for k in [10u32, 20, 50, 100] {
                    assert_eq!(x.pow2k(k).0, reference_pow2k(&x, k).0);
                }
            }
        }
    }

    /// `k` chained squarings must equal one `pow2k(k)`: this is the property
    /// that lets the ladder use `square()` and the inversion use `pow2k`.
    #[test]
    fn chained_square_matches_pow2k() {
        let mut rng = Rng(0x0bad_c0de_0bad_c0de);

        for _ in 0..128 {
            let x = rng.field_element(54);
            let mut chained = x;
            for k in 1..=16u32 {
                chained = chained.square();
                assert_eq!(chained.0, x.pow2k(k).0);
            }
        }
    }

    /// An independent cross-check that does not go through the reference copy:
    /// squaring and multiplying a value by itself must give the same field
    /// element, including at the top of the documented bit excess.
    #[test]
    fn square_agrees_with_mul_at_bit_excess_bound() {
        let mut rng = Rng(0xa5a5_5a5a_a5a5_5a5a);

        for bits in [51u32, 54] {
            for _ in 0..512 {
                let x = rng.field_element(bits);
                assert_eq!(x.square().to_bytes(), (&x * &x).to_bytes());

                let two_x_sq = {
                    let sq = x.square();
                    &sq + &sq
                };
                assert_eq!(x.square2().to_bytes(), two_x_sq.to_bytes());
            }
        }
    }

    /// `mul121666` must agree with the general multiplication by
    /// `APLUS2_OVER_FOUR` limb for limb: it is a specialization of exactly that
    /// product, and the ladder feeds its output straight into the next
    /// operation's bit-excess accounting.
    #[test]
    fn mul121666_matches_general_mul() {
        use crate::backend::serial::u64::constants::APLUS2_OVER_FOUR;

        let mut rng = Rng(0xc0ff_ee00_c0ff_ee00);

        for bits in [51u32, 52, 53, 54] {
            for _ in 0..512 {
                let x = rng.field_element(bits);
                assert_eq!(x.mul121666().0, (&x * &APLUS2_OVER_FOUR).0);
            }
        }

        // Top of the documented bit excess, and the degenerate inputs.
        let edge = FieldElement51([(1u64 << 54) - 1; 5]);
        assert_eq!(edge.mul121666().0, (&edge * &APLUS2_OVER_FOUR).0);

        for limbs in [
            [0u64; 5],
            [1, 0, 0, 0, 0],
            [0, 0, 0, 0, 1],
            [(1 << 54) - 1, 0, 0, 0, 0],
        ] {
            let x = FieldElement51(limbs);
            assert_eq!(x.mul121666().0, (&x * &APLUS2_OVER_FOUR).0);
        }
    }

    /// Independent of the limb comparison above: `mul121666` must equal
    /// multiplying by 121666 built only out of `add`, which shares no code with
    /// either `mul121666` or `mul`. This is what catches the two of them being
    /// wrong in the same way.
    #[test]
    fn mul121666_is_multiplication_by_121666() {
        let mut rng = Rng(0x1357_9bdf_1357_9bdf);

        for _ in 0..64 {
            let x = rng.field_element(51);

            // Double-and-add over the bits of 121666, reducing after every step
            // so the limbs never leave the documented bit excess. `reduce` is
            // the same carry chain `add`'s callers rely on.
            let mut acc = FieldElement51::ZERO;
            for i in (0..17).rev() {
                acc = FieldElement51::reduce((&acc + &acc).0);
                if (121666u32 >> i) & 1 == 1 {
                    acc = FieldElement51::reduce((&acc + &x).0);
                }
            }

            assert_eq!(x.mul121666().to_bytes(), acc.to_bytes());
        }
    }

    /// `mul` is unchanged by this refactor, but pin its output down anyway so
    /// that a future change to the multiplication has a differential test
    /// waiting for it: `(x + y)^2 == x^2 + 2xy + y^2` exercises `mul` and
    /// `square` against each other on random inputs.
    #[test]
    fn mul_and_square_are_consistent() {
        let mut rng = Rng(0x5eed_1234_5eed_1234);

        for _ in 0..512 {
            // Keep the sum inside the documented bit excess: two 53-bit limbs
            // add to at most 2^54.
            let x = rng.field_element(53);
            let y = rng.field_element(53);

            let sum = &x + &y;
            let lhs = sum.square();

            let xy = &x * &y;
            let rhs = &(&x.square() + &y.square()) + &(&xy + &xy);

            assert_eq!(lhs.to_bytes(), rhs.to_bytes());
        }
    }

    /// `sub_unreduced` is not limb-for-limb equal to `Sub` — it returns a
    /// different representative of the same field element — so it is checked on
    /// the two properties the ladder actually depends on: the value is right,
    /// and the limbs land inside the bit excess the next operation requires.
    #[test]
    fn sub_unreduced_agrees_with_sub() {
        let mut rng = Rng(0x5eed_dead_beef_0051);

        // The precondition is "both operands are `mul`/`square` outputs". Rather
        // than assume what those look like, produce them: multiply random
        // elements and subtract the actual results.
        for _ in 0..1024 {
            let x = &rng.field_element(51) * &rng.field_element(51);
            let y = rng.field_element(51).square();

            let fast = x.sub_unreduced(&y);

            // Same field element as the general subtraction.
            assert_eq!(fast.to_bytes(), (&x - &y).to_bytes());

            // The documented postcondition is the tighter `< 2^53`, which is
            // what is asserted; `b < 3` only requires `< 2^54`.
            assert!(fast.0.iter().all(|&l| l < (1u64 << 53)));
        }

        // The worst case the precondition permits: both operands one below the
        // smallest limb of 2p. This is the input that would underflow first if
        // the offset were ever reduced.
        let hi = FieldElement51([4503599627370457; 5]);
        let lo = FieldElement51([0; 5]);
        assert_eq!(hi.sub_unreduced(&lo).to_bytes(), (&hi - &lo).to_bytes());
        assert_eq!(lo.sub_unreduced(&hi).to_bytes(), (&lo - &hi).to_bytes());
        assert!(hi.sub_unreduced(&lo).0.iter().all(|&l| l < (1u64 << 53)));

        // The ladder's *first* iteration, whose operands are not products:
        // `x0` is the identity `(1, 0)` and `x1` is `(from_bytes(u), 1)`. This
        // is the one case the "both operands are mul/square outputs" phrasing
        // does not literally cover, so it is checked on its own.
        let one = FieldElement51::ONE;
        let zero = FieldElement51::ZERO;
        let u = FieldElement51::from_bytes(&[0xff; 32]);
        for (a, b) in [(one, zero), (u, one), (zero, one), (one, u)] {
            assert_eq!(a.sub_unreduced(&b).to_bytes(), (&a - &b).to_bytes());
            assert!(a.sub_unreduced(&b).0.iter().all(|&l| l < (1u64 << 53)));
        }

        // Degenerate cases, including x - x == 0.
        for limbs in [[0u64; 5], [1, 0, 0, 0, 0], [0, 0, 0, 0, 1]] {
            let x = FieldElement51(limbs);
            assert_eq!(
                x.sub_unreduced(&x).to_bytes(),
                FieldElement51::ZERO.to_bytes()
            );
            assert_eq!(x.sub_unreduced(&x).to_bytes(), (&x - &x).to_bytes());
        }
    }
}
