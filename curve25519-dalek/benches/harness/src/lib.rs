//! Benchmark kernels for the X25519 field-arithmetic measurement.
//!
//! The same kernels are compiled for x86_64 and for wasm32, and are driven
//! either by `src/main.rs` (native, `std::time::Instant`) or by `run.mjs`
//! (wasm32, Node's `process.hrtime.bigint()` around an exported function).
//! Keeping one set of kernels is the point: the wasm32 and native numbers then
//! measure literally the same code.
//!
//! Every kernel is a *dependent* chain: the output of one operation feeds the
//! next. That is deliberate. The Montgomery ladder is latency bound, not
//! throughput bound, so a chained kernel is the honest model of what
//! `mul_clamped` sees. It also removes the need for `black_box`, which is not
//! uniformly available across targets.

use core::hint::black_box;
use curve25519_dalek::constants;
use std::sync::OnceLock;

use curve25519_dalek::edwards::{
    EdwardsBasepointTableRadix32, EdwardsBasepointTableRadix64, EdwardsPoint,
};
use curve25519_dalek::montgomery::MontgomeryPoint;
use curve25519_dalek::scalar::Scalar;
use curve25519_dalek::traits::BasepointTable;

use ed25519_dalek::{Signature, Signer, SigningKey, Verifier, VerifyingKey};

#[cfg(curve25519_dalek_bench_internals)]
use curve25519_dalek::bench_internals::FieldElement;

/// Kernel selector. Kept as small integers so the wasm export stays a plain
/// `extern "C"` function with no imports.
pub const K_FE_MUL: u32 = 0;
pub const K_FE_SQUARE: u32 = 1;
pub const K_FE_POW2K50: u32 = 2;
pub const K_MUL_CLAMPED: u32 = 3;
pub const K_MUL_BASE_CLAMPED: u32 = 4;
pub const K_FE_MUL121666: u32 = 5;
pub const K_ED_MUL_BASE: u32 = 6;
pub const K_TO_MONTGOMERY: u32 = 7;
pub const K_FE_INVERT: u32 = 8;
pub const K_ED_MUL_BASE_R64: u32 = 9;
pub const K_ED_MUL_BASE_R32: u32 = 10;
pub const K_ED_VARTIME_DOUBLE: u32 = 11;
pub const K_ED25519_VERIFY: u32 = 12;
pub const K_FE_INVERT_X8: u32 = 13;
pub const K_FE_BATCH_INVERT_8: u32 = 14;
pub const K_FE_INVERT_X16: u32 = 15;
pub const K_FE_BATCH_INVERT_16: u32 = 16;
pub const K_ED_TABLE_CREATE: u32 = 17;
pub const K_ED_MUL_BASE_VEC: u32 = 18;

pub const KERNELS: &[(u32, &str, u32)] = &[
    // (selector, name, field operations per iteration)
    (K_FE_MUL, "fe_mul", 1),
    (K_FE_SQUARE, "fe_square", 1),
    (K_FE_POW2K50, "fe_pow2k50", 50),
    (K_FE_MUL121666, "fe_mul121666", 1),
    (K_FE_INVERT, "fe_invert", 1),
    (K_ED_MUL_BASE, "edwards_mul_base", 1),
    (K_ED_MUL_BASE_R32, "edwards_mul_base_radix32", 1),
    (K_ED_MUL_BASE_R64, "edwards_mul_base_radix64", 1),
    (K_TO_MONTGOMERY, "edwards_to_montgomery", 1),
    (K_ED_VARTIME_DOUBLE, "edwards_vartime_double_base", 1),
    (K_MUL_CLAMPED, "x25519_mul_clamped", 1),
    (K_MUL_BASE_CLAMPED, "x25519_mul_base_clamped", 1),
    (K_ED25519_VERIFY, "ed25519_verify", 1),
    (K_FE_INVERT_X8, "fe_invert_x8", 8),
    (K_FE_BATCH_INVERT_8, "fe_batch_invert_8", 8),
    (K_FE_INVERT_X16, "fe_invert_x16", 16),
    (K_FE_BATCH_INVERT_16, "fe_batch_invert_16", 16),
    (K_ED_TABLE_CREATE, "edwards_table_create", 1),
    (K_ED_MUL_BASE_VEC, "vec_edwards_mul_base", 1),
];

/// Whether the field-level kernels were compiled in. They need
/// `--cfg curve25519_dalek_bench_internals`, because `FieldElement` is
/// `pub(crate)`.
pub const HAS_FIELD_KERNELS: bool = cfg!(curve25519_dalek_bench_internals);

/// Whether the prototype vector fixed-base kernel was compiled in. It needs the
/// AVX2 backend, so it is absent on wasm32 by construction. Kernels named
/// `vec_*` are gated on this the way `fe_*` are gated on `HAS_FIELD_KERNELS`.
pub const HAS_VECTOR_KERNELS: bool = cfg!(all(
    curve25519_dalek_bench_internals,
    curve25519_dalek_backend = "simd",
    target_arch = "x86_64"
));

/// The larger basepoint tables, built once. Constructing a radix-64 table costs
/// on the order of a thousand point operations; leaving that inside the timed
/// region would fold it into every repetition and inflate the very numbers the
/// radix comparison exists to produce.
fn table_radix32() -> &'static EdwardsBasepointTableRadix32 {
    static TABLE: OnceLock<EdwardsBasepointTableRadix32> = OnceLock::new();
    TABLE
        .get_or_init(|| EdwardsBasepointTableRadix32::create(&EdwardsPoint::mul_base(&Scalar::ONE)))
}

fn table_radix64() -> &'static EdwardsBasepointTableRadix64 {
    static TABLE: OnceLock<EdwardsBasepointTableRadix64> = OnceLock::new();
    TABLE
        .get_or_init(|| EdwardsBasepointTableRadix64::create(&EdwardsPoint::mul_base(&Scalar::ONE)))
}

/// Deterministic pseudo-random 32 bytes; fixed so that every target and every
/// backend is fed identical input.
pub fn seed_bytes(seed: u8) -> [u8; 32] {
    let mut bytes = [0u8; 32];
    let mut x = 0x9E37_79B9_7F4A_7C15u64 ^ (seed as u64);
    let mut i = 0;
    while i < 32 {
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        bytes[i..i + 8].copy_from_slice(&x.to_le_bytes());
        i += 8;
    }
    bytes
}

/// Run `iters` iterations of kernel `which`, returning a checksum so nothing
/// can be optimized away.
pub fn run_kernel(which: u32, iters: u32) -> u64 {
    match which {
        K_MUL_CLAMPED => {
            let point = constants::X25519_BASEPOINT;
            let mut s = seed_bytes(1);
            let mut acc = 0u64;
            for _ in 0..iters {
                // Feed the result back in as the next scalar: this keeps a
                // real data dependency between iterations, so the compiler
                // cannot hoist the call out of the loop.
                s = point.mul_clamped(s).to_bytes();
                acc = acc.wrapping_add(s[0] as u64);
            }
            acc
        }
        K_MUL_BASE_CLAMPED => {
            let mut s = seed_bytes(2);
            let mut acc = 0u64;
            for _ in 0..iters {
                s = MontgomeryPoint::mul_base_clamped(s).to_bytes();
                acc = acc.wrapping_add(s[0] as u64);
            }
            acc
        }
        K_ED_MUL_BASE => {
            // The fixed-base half of `mul_base_clamped`, without the
            // Edwards-to-Montgomery conversion. The points are accumulated with
            // Edwards addition rather than compressed, because `compress`
            // performs a field inversion and would dominate what is being
            // measured; one addition is a few field multiplications, under 2%
            // here.
            let mut s = Scalar::from_bytes_mod_order(seed_bytes(5));
            let mut sum = EdwardsPoint::default();
            for _ in 0..iters {
                sum += EdwardsPoint::mul_base(&s);
                s += Scalar::ONE;
            }
            sum.compress().to_bytes()[0] as u64
        }
        K_ED_MUL_BASE_R32 => {
            let table = table_radix32();
            let mut s = Scalar::from_bytes_mod_order(seed_bytes(5));
            let mut sum = EdwardsPoint::default();
            for _ in 0..iters {
                sum += table.mul_base(&s);
                s += Scalar::ONE;
            }
            sum.compress().to_bytes()[0] as u64
        }
        K_ED_MUL_BASE_R64 => {
            // Same as K_ED_MUL_BASE, but against a radix-64 table (120 KB, 43
            // additions) instead of the 30 KB radix-16 table `mul_base` uses.
            // The table is built once, outside the timed loop, the way a
            // consumer holding one in a static would.
            let table = table_radix64();
            let mut s = Scalar::from_bytes_mod_order(seed_bytes(5));
            let mut sum = EdwardsPoint::default();
            for _ in 0..iters {
                sum += table.mul_base(&s);
                s += Scalar::ONE;
            }
            sum.compress().to_bytes()[0] as u64
        }
        K_ED_VARTIME_DOUBLE => {
            // The signature-verification shape (aA + bB). This is one of the
            // few operations the vector backend actually covers, so it is the
            // control for "does AVX2 do anything at all in this build".
            let a = Scalar::from_bytes_mod_order(seed_bytes(7));
            let mut b = Scalar::from_bytes_mod_order(seed_bytes(8));
            let point_a = EdwardsPoint::mul_base(&a);
            let mut sum = EdwardsPoint::default();
            for _ in 0..iters {
                sum += EdwardsPoint::vartime_double_scalar_mul_basepoint(&a, &point_a, &b);
                b += Scalar::ONE;
            }
            sum.compress().to_bytes()[0] as u64
        }
        K_ED_TABLE_CREATE => {
            // `EdwardsBasepointTable::create`: 32 `LookupTable<AffineNielsPoint>`,
            // each converting 8 points to affine. It was the only place in the
            // crate where many independent inversions came from a single call —
            // 256 of them — and so the only candidate for Montgomery's trick.
            // It now performs 32 eight-element batch inversions instead; this
            // kernel is what measured the 3.97x that change is worth, and what
            // would catch it regressing.
            // The input point is built once. An earlier version of this kernel
            // called `mul_base` and `compress` *inside* the loop, which put a
            // fixed-base scalar multiplication and an inversion into every
            // repetition — work `create` does not do. `black_box` on the input
            // is what keeps `create` from being hoisted as loop-invariant now
            // that the point no longer changes; the barrier does the job the
            // varying scalar was doing, without paying for a `mul_base`.
            let p = EdwardsPoint::mul_base(&Scalar::from_bytes_mod_order(seed_bytes(11)));
            let mut last = p;
            for _ in 0..iters {
                let table =
                    curve25519_dalek::edwards::EdwardsBasepointTable::create(black_box(&p));
                // `basepoint()` reads only the first of the 32 sub-tables, so
                // without a barrier the optimizer would be entitled to discard
                // the construction of the other 31 — most of what this kernel
                // exists to measure.
                // A *reference* barrier: it forces the whole table to be
                // materialised without copying its 30 KiB, which a by-value
                // `black_box` would do — and that copy costs one instruction per
                // byte under callgrind, inflating this kernel by 1%.
                let table = black_box(&table);
                // A point copy, not a compression: keeping the result alive
                // must not cost an inversion of its own.
                last = table.basepoint();
            }
            // Compressed once, outside the timed loop, only to consume `last`.
            last.compress().to_bytes()[0] as u64
        }
        K_ED_MUL_BASE_VEC if !HAS_VECTOR_KERNELS => 0,
        #[cfg(all(
            curve25519_dalek_bench_internals,
            curve25519_dalek_backend = "simd",
            target_arch = "x86_64"
        ))]
        K_ED_MUL_BASE_VEC => {
            // Prototype vectorised fixed base (§13.11). The table is built once
            // outside the loop, exactly as a shipped static constant would be.
            let table = curve25519_dalek::bench_internals::VectorBasepointTable::create(
                &EdwardsPoint::mul_base(&Scalar::ONE),
            );
            // A wrong ladder would measure the wrong thing, so it is checked
            // against the serial answer before anything is timed.
            for k in 1..8u8 {
                let s = Scalar::from_bytes_mod_order(seed_bytes(k));
                assert_eq!(
                    table.mul_base(&s).compress(),
                    EdwardsPoint::mul_base(&s).compress(),
                    "vector fixed-base disagrees with serial"
                );
            }
            // Accumulated by Edwards addition, matching `edwards_mul_base`
            // exactly so the two kernels are comparable.
            let mut s = Scalar::from_bytes_mod_order(seed_bytes(5));
            let mut sum = EdwardsPoint::default();
            for _ in 0..iters {
                sum += table.mul_base(&s);
                s += Scalar::ONE;
            }
            sum.compress().to_bytes()[0] as u64
        }
        K_ED25519_VERIFY => {
            // A full Ed25519 signature verification: SHA-512 over the message,
            // then `vartime_double_scalar_mul_basepoint` and `compress`. This is
            // the operation the group-messaging profile is dominated by, and the
            // only hot path in this workspace that reaches the vector backend.
            //
            // Key, signature and message are built once, outside the timed loop.
            // Verification is a pure function of them, so re-verifying the same
            // signature measures exactly the work a consumer does per message;
            // what must not leak into the timing is the *signing*, which is a
            // different operation entirely.
            let signing = SigningKey::from_bytes(&seed_bytes(9));
            let verifying: VerifyingKey = signing.verifying_key();
            let message = seed_bytes(10);
            let signature: Signature = signing.sign(&message);

            // The verification is loop-invariant, so under this crate's fat LTO
            // the optimizer would be entitled to perform it once and turn the
            // loop into repeated additions. (It demonstrably does not: the
            // marginal cost of an iteration is ~316k instructions, the right
            // order for a verification, and it moved when the AVX2 backend
            // changed. The barriers make that a guarantee rather than an
            // observation.)
            let mut acc = 0u64;
            for _ in 0..iters {
                let ok = black_box(&verifying)
                    .verify(black_box(&message), black_box(&signature))
                    .is_ok();
                acc = acc.wrapping_add(black_box(ok) as u64);
            }
            // Fail loudly rather than silently timing a rejected signature.
            assert_eq!(acc, iters as u64, "verification failed");
            acc
        }
        K_TO_MONTGOMERY => {
            // The conversion half: one field inversion plus a multiplication.
            // The point is perturbed by a cheap Edwards addition each iteration
            // so the conversion cannot be hoisted out of the loop.
            let mut p = EdwardsPoint::mul_base(&Scalar::from_bytes_mod_order(seed_bytes(6)));
            let basepoint = EdwardsPoint::mul_base(&Scalar::ONE);
            let mut acc = 0u64;
            for _ in 0..iters {
                acc = acc.wrapping_add(p.to_montgomery().to_bytes()[0] as u64);
                p += basepoint;
            }
            acc
        }
        _ => run_field_kernel(which, iters),
    }
}

#[cfg(curve25519_dalek_bench_internals)]
fn run_field_kernel(which: u32, iters: u32) -> u64 {
    let mut bytes_a = seed_bytes(3);
    let mut bytes_b = seed_bytes(4);
    bytes_a[31] &= 0x7f;
    bytes_b[31] &= 0x7f;
    let mut x = FieldElement::from_bytes(&bytes_a);
    let y = FieldElement::from_bytes(&bytes_b);

    match which {
        K_FE_MUL => {
            for _ in 0..iters {
                x = &x * &y;
            }
        }
        K_FE_SQUARE => {
            for _ in 0..iters {
                x = x.square();
            }
        }
        K_FE_POW2K50 => {
            for _ in 0..iters {
                x = x.pow2k(50);
            }
        }
        K_FE_MUL121666 => {
            for _ in 0..iters {
                x = x.mul121666();
            }
        }
        K_FE_INVERT => {
            for _ in 0..iters {
                x = curve25519_dalek::bench_internals::invert(&x);
            }
        }
        // The batch-versus-repeated comparison. `n` independent inversions
        // against one `invert_batch` over the same `n`: Montgomery's trick
        // replaces n inversions with 1 inversion and 3(n-1) multiplications.
        //
        // n = 8 is the size of a lookup table, the one place in this crate where
        // several independent inversions sit next to each other — that batch is
        // now taken, and this pair is what prices it. n = 16 is the count the
        // group profile attributes to `as_affine`; that one is *not* batchable,
        // and pricing it is what makes the refusal checkable rather than
        // asserted.
        //
        // The chain between iterations is preserved by folding limb 0 of the
        // first output back into every input, so the batch cannot be hoisted.
        K_FE_INVERT_X8 | K_FE_INVERT_X16 => {
            let n = if which == K_FE_INVERT_X8 { 8 } else { 16 };
            for _ in 0..iters {
                let mut acc = FieldElement::ZERO;
                for i in 0..n {
                    let t = &x + &FieldElement::from_bytes(&seed_bytes(0x40 + i as u8));
                    acc = &acc + &curve25519_dalek::bench_internals::invert(&t);
                }
                x = acc;
            }
        }
        K_FE_BATCH_INVERT_8 | K_FE_BATCH_INVERT_16 => {
            let n = if which == K_FE_BATCH_INVERT_8 { 8 } else { 16 };
            let mut buf = vec![FieldElement::ZERO; n];
            for _ in 0..iters {
                for (i, slot) in buf.iter_mut().enumerate() {
                    *slot = &x + &FieldElement::from_bytes(&seed_bytes(0x40 + i as u8));
                }
                curve25519_dalek::bench_internals::invert_batch(&mut buf);
                let mut acc = FieldElement::ZERO;
                for slot in buf.iter() {
                    acc = &acc + slot;
                }
                x = acc;
            }
        }
        _ => return 0,
    }

    let out = x.to_bytes();
    u64::from_le_bytes([
        out[0], out[1], out[2], out[3], out[4], out[5], out[6], out[7],
    ])
}

/// Without the internal hook the field kernels are not measurable; return 0 so
/// the driver reports them as unavailable rather than failing to build.
#[cfg(not(curve25519_dalek_bench_internals))]
fn run_field_kernel(_which: u32, _iters: u32) -> u64 {
    0
}

/// What `curve25519-dalek` *actually* compiled, read out of the dependency
/// itself rather than out of the flags we hoped we passed.
///
/// bit 0: 64-bit limbs · bit 1: fiat backend · bit 2: simd backend ·
/// bit 3: field kernels available.
///
/// Bit 3 doubles as the validity flag for bits 0-2. The limb width and backend
/// are read out of `bench_internals`, so without that hook they are not merely
/// unavailable but *unknown*, and a driver must not decode a cleared bit 0 as
/// "32-bit limbs". The crate's configuration cfgs are set by its own build
/// script and are not visible to a dependent, so there is no second source to
/// fall back on — reproducing the selection logic here would attest to what
/// this file believes rather than to what the dependency compiled, which is
/// the whole point of reading it out of the dependency.
pub fn config_code_impl() -> u32 {
    // Nothing assigns to `code` when the hook is absent, which is the whole
    // point of the unknown encoding above; that build mode should still be
    // warning-free.
    #[allow(unused_mut)]
    let mut code = 0;
    #[cfg(curve25519_dalek_bench_internals)]
    {
        use curve25519_dalek::bench_internals as bi;
        if bi::LIMB_BITS == 64 {
            code |= 1;
        }
        if bi::BACKEND == "fiat" {
            code |= 2;
        }
        if bi::BACKEND == "simd" || bi::BACKEND == "avx512" {
            code |= 4;
        }
        code |= 8;
    }
    code
}

/// wasm32 entry point. Node instantiates the module and times this call.
#[cfg(target_arch = "wasm32")]
#[unsafe(no_mangle)]
pub extern "C" fn bench_kernel(which: u32, iters: u32) -> u64 {
    run_kernel(which, iters)
}

#[cfg(target_arch = "wasm32")]
#[unsafe(no_mangle)]
pub extern "C" fn config_code() -> u32 {
    config_code_impl()
}
