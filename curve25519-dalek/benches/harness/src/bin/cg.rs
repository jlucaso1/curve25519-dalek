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

    std::hint::black_box(x25519_field_harness::run_kernel(which, iters));
}
