# M2.11.1.5 performance pass: report

Unattended run, 2026-11-10. Benchmark: `bench_layers --mode moist --levels 0.25,0.5,0.85`
(see `docs/BENCHMARKS.md`). "ms/step" below is 2,562 cells, 8 threads, quick
runs (10 days x 3, median) unless marked formal.

Formal baseline (`775c84a`): 2.5k cells 11.19 ms/step (8 threads), 11.63 (15),
35.53 (1); 642 cells 6.37 (4 threads), 10.85 (1).
Climate baseline (642 cells, `FMS_SLAB=10 climate 8 4`, final year):
T_s 291.194 K, planetary albedo 0.315, OLR 233.3, P 3.750 mm/day.

## Changes

| Commit | Change | Before | After | Hash |
|---|---|---|---|---|
| `f74cd1c` | `lift`: reuse the saturation humidity between steps; stop the convective reference once the parcel is below the stratosphere floor | 11.19 | 9.39 (1 thread: 35.5 -> 24.5) | unchanged |
| `fb597f3` | Semi-implicit step: vertical modes solved in parallel, per-layer loops in parallel | 9.39 | 8.72 | unchanged |
| `d2d6cfe` | Tendency, hyperdiffusion and energy rates: layers in parallel (per-layer scratch) | 8.72 | 7.47 | unchanged |
| `eaaf764` | Humidity transport: layers in parallel | 7.47 | 7.05 | unchanged |
| `ab134a2` | Longwave: each slab's transmissivity computed once per term instead of four times | 7.05 | 7.00 (1 thread: 24.5 -> 23.8) | unchanged |
| `906ea95` | Moist physics (the one hash change): `lift` steps pressure by a ratio (no `exp(ln p)`) and carries saturation between calls; in-layer profile factors `(σ/σ̄ₖ)^β` memoized per `Levels` (no `pow` per node); `planck_tail` gets `e^{-nx}` by multiplication | 7.00 | 5.91 (1 thread: 23.8 -> 18.6) | updated to `0xe70e_3500_47bd_d6b0`; climate (642 cells, 4 years, last year): T_s 291.171 vs 291.194 K (-0.02), albedo 0.314 vs 0.315, OLR 233.4 vs 233.3, P 3.765 vs 3.750 mm/day |
| `c482ee5` | Column physics split into chunks of 32 columns, not 256 (uneven column cost); its tendency assembly per layer in parallel | 5.91 | 5.63 | unchanged |
| `c482ee5` | `LayeredModel::step_tendencies`: the forcing (column physics) runs beside the tendency and hyperdiffusion (`rayon::join`) | 5.63 | 5.06 (15 threads: 6.25 -> 5.28) | unchanged |
| `2ad417a` | `LayeredSemiImplicit::step` runs on a pool thread (`rayon::scope`); `dot`'s chunk and the column chunk get their own constants (`DOT_CHUNK`, `MIN_COLUMNS`) | 5.06 | 4.87 (15 threads: 5.28 -> 4.94) | unchanged |
| (next) | Humidity interface ratios computed beside the Helmholtz solves; longwave scratch allocated once per column | 4.87 | 4.63 | unchanged |

## Tried and rejected

- Coarser `lift` steps (0.1, 0.2 in ln p) or RK4: the scheme is first order
  through the condensation level, so the reference profile moves by ~1 K.
  Not done.

- `MIN_CHUNK` 1024 with the step in the pool: 4.75 ms at 8 threads but 5.70 at
  15 (from 4.94). Kept 4096.

- Stopping the convective lift where no layer above can be unstable (exact):
  4.87 -> 4.80 only; nearly every column is tested up to the top layer, whose
  reference is floored at the stratosphere temperature. Reverted.

## Unfinished / uncommitted

## Observations

- `KineticHeating::rate` (`crates/atmos/src/dissipation.rs`) sums with a parallel `.sum()`; above 4,096 edges its float result can depend on how rayon splits. It only feeds the `dissipation` diagnostic. Not touched.

## Architectural suggestions

## Missing tools or permissions

- None so far.
