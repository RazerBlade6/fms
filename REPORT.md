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
| (this) | `lift`: reuse the saturation humidity between steps; stop the convective reference once the parcel is below the stratosphere floor | 11.19 | 9.39 (1 thread: 35.5 -> 24.5) | unchanged |

## Tried and rejected

- Coarser `lift` steps (0.1, 0.2 in ln p) or RK4: the scheme is first order
  through the condensation level, so the reference profile moves by ~1 K.
  Not done.

## Unfinished / uncommitted

## Architectural suggestions

## Missing tools or permissions

- None so far.
