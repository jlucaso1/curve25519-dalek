//! Native driver for the X25519 field-arithmetic kernels.
//!
//! Criterion reports a bootstrapped mean, which is the right statistic on a
//! quiet machine and a misleading one on a shared virtual machine, where the
//! mean tracks the hypervisor's mood. This driver instead runs each kernel
//! `--reps` times and reports the **minimum**, which is the closest thing to a
//! noise-free estimate available: interference can only ever make a run slower.
//! The median and the spread are printed alongside so the reader can see how
//! noisy the host was.
//!
//! Output is one TSV line per kernel:
//!
//! ```text
//! kernel  iters  reps  min_ns_per_op  median_ns_per_op  max_ns_per_op  spread_pct
//! ```

use std::time::Instant;

use x25519_field_harness::{HAS_FIELD_KERNELS, KERNELS, config_code_impl, run_kernel};

/// Target wall time for a single repetition, in nanoseconds. Long enough to
/// dwarf the timer's resolution, short enough that a scheduling hiccup lands in
/// only a few repetitions.
const TARGET_NS: u128 = 20_000_000;

fn arg_value(name: &str, default: u32) -> u32 {
    let mut args = std::env::args().skip(1);
    while let Some(a) = args.next() {
        if a == name
            && let Some(v) = args.next()
        {
            return v.parse().unwrap_or(default);
        }
    }
    default
}

/// Find an iteration count whose wall time is close to `TARGET_NS`.
fn calibrate(which: u32) -> u32 {
    let mut iters = 1u32;
    loop {
        let start = Instant::now();
        std::hint::black_box(run_kernel(which, iters));
        let elapsed = start.elapsed().as_nanos().max(1);
        if elapsed >= TARGET_NS || iters >= 1 << 28 {
            return iters;
        }
        // Grow towards the target, but never by more than 16x at a time so a
        // single slow probe cannot overshoot wildly.
        let factor = ((TARGET_NS / elapsed) as u32).clamp(2, 16);
        iters = iters.saturating_mul(factor);
    }
}

fn main() {
    let reps = arg_value("--reps", 15).max(3);
    let code = config_code_impl();

    println!("# config_code={code} field_kernels={HAS_FIELD_KERNELS}");
    println!("kernel\titers\treps\tmin_ns_op\tmed_ns_op\tmax_ns_op\tspread_pct");

    for &(which, name, ops_per_iter) in KERNELS {
        if !HAS_FIELD_KERNELS && name.starts_with("fe_") {
            println!("{name}\t-\t-\tunavailable\t-\t-\t-");
            continue;
        }

        let iters = calibrate(which);
        let mut samples = Vec::with_capacity(reps as usize);
        for _ in 0..reps {
            let start = Instant::now();
            std::hint::black_box(run_kernel(which, iters));
            let elapsed = start.elapsed().as_nanos();
            let ops = (iters as u128) * (ops_per_iter as u128);
            samples.push(elapsed as f64 / ops as f64);
        }
        samples.sort_by(|a, b| a.partial_cmp(b).expect("no NaN timings"));

        let min = samples[0];
        let med = samples[samples.len() / 2];
        let max = samples[samples.len() - 1];
        let spread = 100.0 * (max - min) / min;
        println!("{name}\t{iters}\t{reps}\t{min:.4}\t{med:.4}\t{max:.4}\t{spread:.1}");
    }
}
