//! Benchmarks for the field arithmetic that X25519 actually spends its time
//! in, plus the end-to-end X25519 operations a consumer of this crate calls.
//!
//! The field-level benchmarks require the internal hook:
//!
//! ```sh
//! RUSTFLAGS='--cfg curve25519_dalek_bench_internals' \
//!     cargo bench --bench x25519_field
//! ```
//!
//! Without that cfg only the end-to-end group is built, so the bench target
//! still compiles in a plain `cargo bench`.
//!
//! Backends are selected the usual way, e.g.
//!
//! ```sh
//! RUSTFLAGS='--cfg curve25519_dalek_backend="fiat" --cfg curve25519_dalek_bench_internals' \
//!     cargo bench --bench x25519_field
//! ```

use criterion::{Criterion, criterion_main};
use std::hint::black_box;

use curve25519_dalek::constants;
use curve25519_dalek::montgomery::MontgomeryPoint;

/// Deterministic 32-byte string, so every backend and every target is fed the
/// same input and runs are comparable across machines.
fn scalar_bytes(seed: u8) -> [u8; 32] {
    let mut bytes = [0u8; 32];
    // A cheap xorshift-ish fill; the values only need to be fixed and
    // "random-looking", they are not used for anything secret.
    let mut x = 0x9E37_79B9_7F4A_7C15u64 ^ (seed as u64);
    for chunk in bytes.chunks_mut(8) {
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        chunk.copy_from_slice(&x.to_le_bytes());
    }
    bytes
}

/// End-to-end X25519, i.e. what a protocol implementation calls.
mod x25519_benches {
    use super::*;

    pub(crate) fn benches(c: &mut Criterion) {
        let mut g = c.benchmark_group("x25519");

        // Variable-base: the Montgomery ladder, entirely serial field
        // arithmetic. This is the `x25519_agreement`-shaped operation.
        g.bench_function("MontgomeryPoint::mul_clamped", |b| {
            let point = constants::X25519_BASEPOINT;
            let s = scalar_bytes(1);
            b.iter(|| black_box(point).mul_clamped(black_box(s)))
        });

        // Fixed-base: goes through EdwardsPoint::mul_base, so on x86_64 this
        // one *does* reach the vectorized backend, and then converts.
        g.bench_function("MontgomeryPoint::mul_base_clamped", |b| {
            let s = scalar_bytes(2);
            b.iter(|| MontgomeryPoint::mul_base_clamped(black_box(s)))
        });

        g.finish();
    }
}

/// Isolated field operations. Diagnostic for the end-to-end numbers above:
/// a change to `mul`/`pow2k` has to show up here first, and then survive the
/// serial dependency chain of the ladder.
#[cfg(curve25519_dalek_bench_internals)]
mod field_benches {
    use super::*;
    use curve25519_dalek::bench_internals::{
        FieldElement, MULS_PER_MUL_CLAMPED, SQUARES_PER_MUL_CLAMPED, invert,
    };

    fn fe(seed: u8) -> FieldElement {
        let mut bytes = scalar_bytes(seed);
        bytes[31] &= 0x7f;
        FieldElement::from_bytes(&bytes)
    }

    pub(crate) fn benches(c: &mut Criterion) {
        let mut g = c.benchmark_group("field");

        let a = fe(3);
        let b_ = fe(4);

        g.bench_function("FieldElement::mul", |b| {
            b.iter(|| black_box(&a) * black_box(&b_))
        });

        g.bench_function("FieldElement::square", |b| {
            b.iter(|| black_box(&a).square())
        });

        g.bench_function("FieldElement::square2", |b| {
            b.iter(|| black_box(&a).square2())
        });

        g.bench_function("FieldElement::add", |b| {
            b.iter(|| black_box(&a) + black_box(&b_))
        });

        g.bench_function("FieldElement::sub", |b| {
            b.iter(|| black_box(&a) - black_box(&b_))
        });

        // `pow2k(k)` exists to amortize the per-call setup of a squaring over
        // `k` squarings. Measuring a few `k` shows how much of `square()` is
        // amortizable overhead and how much is the squaring itself. Compare
        // `pow2k(k)` against `k` times the `square` number above.
        for k in [1u32, 2, 4, 10, 50, 100] {
            g.bench_function(format!("FieldElement::pow2k({k})"), |b| {
                b.iter(|| black_box(&a).pow2k(black_box(k)))
            });
        }

        // Same total number of squarings, expressed as repeated `square()`
        // calls instead of one `pow2k`. The ratio against `pow2k(50)` is the
        // amortization pow2k buys.
        g.bench_function("FieldElement::square x50", |b| {
            b.iter(|| {
                let mut x = *black_box(&a);
                for _ in 0..50 {
                    x = x.square();
                }
                x
            })
        });

        // `invert` is the tail of every `mul_clamped`: 254 squarings + 12
        // multiplications, and the only place `pow2k(k > 1)` is used on the
        // X25519 path.
        g.bench_function("FieldElement::invert", |b| b.iter(|| invert(black_box(&a))));

        // A synthetic stand-in for the ladder's arithmetic mix: it performs the
        // same number of muls and squares as one `mul_clamped`, but without the
        // conditional swaps and with a shorter dependency chain. The gap
        // between this and `x25519/MontgomeryPoint::mul_clamped` is everything
        // the ladder costs on top of raw field arithmetic.
        g.bench_function("field mix of one mul_clamped", |b| {
            b.iter(|| {
                let mut x = *black_box(&a);
                let y = *black_box(&b_);
                for _ in 0..MULS_PER_MUL_CLAMPED {
                    x = &x * &y;
                }
                for _ in 0..SQUARES_PER_MUL_CLAMPED {
                    x = x.square();
                }
                x
            })
        });

        g.finish();
    }
}

#[cfg(not(curve25519_dalek_bench_internals))]
mod field_benches {
    pub(crate) fn benches(_c: &mut criterion::Criterion) {}
}

fn all() {
    let mut c = Criterion::default().configure_from_args();
    x25519_benches::benches(&mut c);
    field_benches::benches(&mut c);
}

criterion_main!(all);
