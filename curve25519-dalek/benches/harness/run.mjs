// wasm32 driver for the X25519 field-arithmetic kernels.
//
// Criterion does not run on wasm32-unknown-unknown, so the module is
// instantiated directly and timed from the host with process.hrtime.bigint().
// The module has no imports and no wasm-bindgen glue, so nothing but the
// kernels themselves is inside the timed region except one cross-boundary call
// per repetition, which is amortized over a calibrated iteration count.
//
// The reported statistic is the MINIMUM over repetitions, matching the native
// driver in src/main.rs: interference can only make a run slower.
//
// Usage:
//   node run.mjs <module.wasm> [reps]

import { readFile } from "node:fs/promises";

const KERNELS = [
  // [selector, name, field operations per iteration]
  [0, "fe_mul", 1],
  [1, "fe_square", 1],
  [2, "fe_pow2k50", 50],
  [5, "fe_mul121666", 1],
  [8, "fe_invert", 1],
  [6, "edwards_mul_base", 1],
  [10, "edwards_mul_base_radix32", 1],
  [9, "edwards_mul_base_radix64", 1],
  [7, "edwards_to_montgomery", 1],
  [3, "x25519_mul_clamped", 1],
  [4, "x25519_mul_base_clamped", 1],
];

const TARGET_NS = 20_000_000n;

const path = process.argv[2];
if (!path) {
  console.error("usage: node run.mjs <module.wasm> [reps]");
  process.exit(2);
}

// Mirror the native driver's parsing (src/main.rs): a u32, falling back to 15
// on anything unparseable, and never fewer than 3 repetitions. Without the
// floor, `reps < 1` would leave `samples` empty and the reporting below would
// read `undefined`.
const repsArg = process.argv[3] ?? "15";
const reps =
  /^\d+$/.test(repsArg) && Number(repsArg) <= 0xffff_ffff
    ? Math.max(3, Number(repsArg))
    : 15;

const bytes = await readFile(path);
const { instance } = await WebAssembly.instantiate(bytes, {});
const { bench_kernel, config_code } = instance.exports;

const code = config_code ? config_code() : 0;
const hasFieldKernels = (code & 8) !== 0;
console.log(`# config_code=${code} field_kernels=${hasFieldKernels}`);
console.log("kernel\titers\treps\tmin_ns_op\tmed_ns_op\tmax_ns_op\tspread_pct");

function timeOnce(which, iters) {
  const t0 = process.hrtime.bigint();
  const out = bench_kernel(which, iters);
  const t1 = process.hrtime.bigint();
  // Consume the checksum so a future engine cannot elide the call.
  if (out === 0xdeadbeefn) console.error("impossible");
  return t1 - t0;
}

function calibrate(which) {
  let iters = 1;
  for (;;) {
    const ns = timeOnce(which, iters) || 1n;
    if (ns >= TARGET_NS || iters >= 1 << 28) return iters;
    const factor = Math.min(16, Math.max(2, Number(TARGET_NS / ns)));
    iters = Math.min(1 << 28, iters * factor);
  }
}

for (const [which, name, opsPerIter] of KERNELS) {
  if (!hasFieldKernels && name.startsWith("fe_")) {
    console.log(`${name}\t-\t-\tunavailable\t-\t-\t-`);
    continue;
  }
  const iters = calibrate(which);
  const samples = [];
  for (let i = 0; i < reps; i++) {
    const ns = Number(timeOnce(which, iters));
    samples.push(ns / (iters * opsPerIter));
  }
  samples.sort((a, b) => a - b);
  const min = samples[0];
  const med = samples[samples.length >> 1];
  const max = samples[samples.length - 1];
  const spread = (100 * (max - min)) / min;
  console.log(
    `${name}\t${iters}\t${reps}\t${min.toFixed(4)}\t${med.toFixed(4)}\t` +
      `${max.toFixed(4)}\t${spread.toFixed(1)}`,
  );
}
