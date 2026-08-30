//! Profiling entry point: runs one kernel a fixed number of times and exits.
//!
//! Exists so that the instruction attribution in
//! `docs/perf-x25519-field-arithmetic.md` can be reproduced. Timing on a shared
//! virtual machine is noisy enough that the minimum over repetitions is the only
//! statistic worth reporting; instruction counts have no such problem, so
//! `callgrind` answers "where does the work go" exactly rather than
//! statistically.
//!
//! ```sh
//! cargo build --release --bin cg
//! valgrind --tool=callgrind --callgrind-out-file=cg.out \
//!     target/release/cg 3 20          # kernel 3 = x25519_mul_clamped
//! callgrind_annotate cg.out
//! ```
//!
//! The kernel selectors are the `K_*` constants in `lib.rs`.

fn main() {
    let mut args = std::env::args().skip(1);

    // Defaults apply only when an argument is absent. A typo must not silently
    // profile a different kernel, and `iters = 0` must not silently profile
    // nothing at all.
    let which: u32 = match args.next() {
        Some(a) => a.parse().expect("kernel selector must be a u32"),
        None => 3,
    };
    let iters: u32 = match args.next() {
        Some(a) => a.parse().expect("iteration count must be a u32"),
        None => 20,
    };
    assert!(iters > 0, "iteration count must be greater than zero");

    // `run_kernel` routes any unrecognized selector to the field-kernel
    // dispatcher, which without the internals hook returns zero immediately. A
    // typo, or a field kernel asked for in a build that cannot run one, would
    // therefore produce a *successful* callgrind run over no work at all — the
    // worst failure mode for a measurement tool, since the output looks real.
    let (_, name, _) = x25519_field_harness::KERNELS
        .iter()
        .find(|(selector, _, _)| *selector == which)
        .unwrap_or_else(|| panic!("unknown kernel selector {which}; see K_* in lib.rs"));
    assert!(
        x25519_field_harness::HAS_FIELD_KERNELS || !name.starts_with("fe_"),
        "kernel {name} needs --cfg curve25519_dalek_bench_internals"
    );

    std::hint::black_box(x25519_field_harness::run_kernel(which, iters));
}
