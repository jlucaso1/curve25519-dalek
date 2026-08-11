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
import { fileURLToPath } from "node:url";
import { dirname, join } from "node:path";

// The kernel list is *derived from* `src/lib.rs`, not copied from it.
//
// It used to be a hand-maintained array here, and it silently fell behind:
// kernels 12-17 (`ed25519_verify`, the batch-inversion comparisons and
// `edwards_table_create`) were added on the Rust side and never mirrored, so
// every wasm run quietly measured twelve of eighteen kernels while the module
// docs claim one set of kernels for both targets. Parsing the Rust table makes
// that drift impossible rather than merely fixed once, and a parse that finds
// nothing is a hard error — silently running a subset is the failure being
// removed.
async function loadKernels() {
  const libPath = join(dirname(fileURLToPath(import.meta.url)), "src", "lib.rs");
  const src = await readFile(libPath, "utf8");

  const consts = new Map();
  for (const m of src.matchAll(/pub const (K_[A-Z0-9_]+): u32 = (\d+);/g)) {
    consts.set(m[1], Number(m[2]));
  }

  const table = src.match(
    /pub const KERNELS: &\[\(u32, &str, u32\)\] = &\[([\s\S]*?)\n\];/,
  );
  if (!table) throw new Error(`could not find KERNELS in ${libPath}`);

  const kernels = [];
  for (const m of table[1].matchAll(/\(\s*(K_[A-Z0-9_]+),\s*"([^"]+)",\s*(\d+)\s*\)/g)) {
    const which = consts.get(m[1]);
    if (which === undefined) throw new Error(`unresolved selector ${m[1]}`);
    kernels.push([which, m[2], Number(m[3])]);
  }
  if (kernels.length === 0) throw new Error(`parsed no kernels from ${libPath}`);
  return kernels;
}

const KERNELS = await loadKernels();

const TARGET_NS = 20_000_000n;

const path = process.argv[2];
if (!path) {
  console.error("usage: node run.mjs <module.wasm> [reps]");
  process.exit(2);
}

// Mirror the native driver's parsing (src/main.rs): a u32, falling back to 15
// on anything unparsable, and never fewer than 3 repetitions. Without the
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

// Bits, as set by `config_code_impl` in lib.rs: 1 = 64-bit limbs,
// 2 = fiat backend, 4 = simd/avx512, 8 = field kernels compiled in.
//
// Bit 3 is also the validity flag for bits 0-2: they are read out of
// `bench_internals`, so a module built without `--cfg
// curve25519_dalek_bench_internals` reports zero for all of them and the limb
// width and backend are unknown, not 32-bit serial. Printing a guess there
// would put a false configuration line above a real measurement.
const code = config_code ? config_code() : 0;
const hasFieldKernels = (code & 8) !== 0;
const limbBits = hasFieldKernels ? (code & 1 ? 64 : 32) : "unknown";
const backend = !hasFieldKernels
  ? "unknown"
  : code & 2
    ? "fiat"
    : code & 4
      ? "simd"
      : "serial";
console.log(
  `# limb_bits=${limbBits} backend=${backend} ` +
    `field_kernels=${hasFieldKernels} (config_code=${code})`,
);
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
