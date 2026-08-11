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

use curve25519_dalek::constants;
use std::sync::OnceLock;

use curve25519_dalek::edwards::{
    EdwardsBasepointTableRadix32, EdwardsBasepointTableRadix64, EdwardsPoint,
};
use curve25519_dalek::montgomery::MontgomeryPoint;
use curve25519_dalek::scalar::Scalar;
use curve25519_dalek::traits::BasepointTable;

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
];

/// Whether the field-level kernels were compiled in. They need
/// `--cfg curve25519_dalek_bench_internals`, because `FieldElement` is
/// `pub(crate)`.
pub const HAS_FIELD_KERNELS: bool = cfg!(curve25519_dalek_bench_internals);

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
