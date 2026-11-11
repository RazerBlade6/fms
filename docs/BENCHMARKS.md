# Benchmarks

Criterion suites for the performance-critical building blocks.

```sh
cargo bench -p fms-grid --bench geo            # grids, geometry, transfers, point location
cargo bench -p fms-sched --bench sched        # scheduler overhead
cargo bench -p fms-planet --bench planet      # insolation, orbit
```

`-- --quick` gives a fast estimate. To check a change for regressions, save a
baseline first:

```sh
cargo bench -p fms-grid --bench geo -- --save-baseline before
# … make the change …
cargo bench -p fms-grid --bench geo -- --baseline before
```

HTML reports go to `target/criterion/`.

## Baseline

Taken on 2026-11-07 on the development Mac (Apple silicon), release profile,
`--quick`. Sizes: f45 = 20,252 cells (default climate grid); f160 = 256,002 cells
(default CLI terrain grid).

| Benchmark                             | Time    | Notes                                                                      |
| ------------------------------------- | ------- | -------------------------------------------------------------------------- |
| `build/geodesic_f160`                 | 11.3 ms | Subdivision only                                                           |
| `build/centroidal_f45`                | 406 ms  | Multilevel Lloyd to 1e-6                                                   |
| `build/geometry_f160`                 | 54.7 ms | Metric terms, TRiSK and Perot weights                                      |
| `build/hilbert_order_f160`            | 89.0 ms | Renumbering and rebuilding adjacency                                       |
| `two_level_f160_to_f45/build_map`     | 3.96 ms |                                                                            |
| `two_level_f160_to_f45/restrict_mean` | 155 µs  |                                                                            |
| `two_level_f160_to_f45/restrict_sum`  | 154 µs  |                                                                            |
| `two_level_f160_to_f45/prolong`       | 214 µs  |                                                                            |
| `operators_f45/perot_reconstruction`  | 149 µs  |                                                                            |
| `locate_f160/cells_10k_points`        | 636 µs  | ~64 ns per point                                                           |
| `locate_f160/triangles_10k_points`    | 1.39 ms |                                                                            |
| `planet/insolation_20k`               | 128 µs  | Parallel                                                                   |
| `planet/daily_mean_20k`               | 3.07 ms | Serial loop in the bench; the system runs it in parallel                   |
| `planet/orbit_state`                  | 93 ns   | Kepler solve                                                               |
| `scheduler/step_20_systems`           | 119 µs  | 20 trivial systems in 2 stages; see [OBSERVATIONS.md](OBSERVATIONS.md) O-8 |

## Conditions

Every benchmark is run one configuration at a time, never alongside tests or
builds, with repeats and the median reported, and logged here with:
the machine, power source and mode, what else is running, CPU idle and free
memory at the start (`top -l 2`) and end, macOS and rustc versions, the
commits compared and the date.

## Model cost by number of layers (M2.9.5)

`cargo run --release -p fms-atmos --example bench_layers -- --mode moist|hs --levels … --days 30`
(step 9's code: the same program against commit `6692284`, two layers only).
2,562 cells (k = 16), 30 model days per run, 3 runs each, median shown; the
runs agreed to within 1.5%.

**Conditions** (2026-11-08):

- MacBook Pro, Apple M5 Pro: 15 cores (5 "Super", 10 "Performance"), 24 GB;
  all 15 used (rayon).
- macOS 26.6.2 (25G83); rustc 1.99.1; release profile
  (`debug = "line-tables-only"`).
- On AC power, High Power mode (`pmset` powermode 2).
- Running: only Activity Monitor and the Ghostty terminal; memory-heavy apps
  closed.
- CPU idle 98.2% before, 99.3% after (`top`). Memory: 9.4 GB unused
- Commits: step 9 `6692284`; step 9.5 `6dd3579` plus the benchmark example.

| Configuration                        | Moist aquaplanet (dt 2057 s) | Held–Suarez, dry (dt 3086 s) |
| ------------------------------------ | ---------------------------- | ---------------------------- |
| Step 9, 2 layers                     | 6.60 ms/step (3.6 days/s)    | 2.75 ms/step (13.0 days/s)   |
| 9.5, 2 layers, step 9's reference    | 7.31 (+11%)                  | 3.09 (+12%)                  |
| 9.5, 2 layers (isothermal reference) | 7.41 (+12%)                  | 3.14 (+14%)                  |
| 9.5, 3 layers (σ 0.4, 0.85)          | 9.46 (×1.28 of 9.5's two)    | 4.33 (×1.38)                 |
| 9.5, 5 equal layers                  | 13.66 (×1.84)                | 6.85 (×2.18)                 |

- The generalization costs ~12% at two layers with identical arithmetic
  (same results bit for bit): per-layer vectors, per-cell allocations in the
  physics, and scratch copies ([OBSERVATIONS.md](OBSERVATIONS.md) O-29).
- The isothermal reference adds ~1%: one more CG iteration per solve
  (12.1 against 11.0, moist).
- Extra layers cost less than in proportion: the moist model's third layer
  adds 28% and five layers 84%, since radiation's spectral work, the
  Helmholtz solves of the slow internal modes (9–12 iterations) and
  per-cell overheads don't triple. A year on three layers at this resolution:
  about 2.5 minutes of stepping, against 1.7 for step 9.

## Early performance pass (M2.11.1.5)

`RAYON_NUM_THREADS=T cargo run --release -p fms-atmos --example bench_layers -- --mode moist --levels 0.25,0.5,0.85 --k K --days D --repeats 3`:
the moist aquaplanet on the climate's four layers, median of 3 runs
(2,562 cells: dt 2057 s, 15 days; 642 cells: dt 4114 s, 60 days).

**Conditions** (2026-11-10): the
machine, macOS and rustc as above; on AC power,
High Power mode (`pmset` powermode 2), `caffeinate` active; nothing else
running but the terminal and the agent driving the run; CPU idle 98.6%
before and 98.4% after the baseline (`top -l 2`).

| Commit                                                                  | 2,562 cells, 8 threads | 15 threads | 1 thread | 642 cells, 4 threads | 1 thread |
| ----------------------------------------------------------------------- | ---------------------- | ---------- | -------- | -------------------- | -------- |
| `775c84a` (baseline)                                                    | 11.19 ms/step          | 11.63      | 35.53    | 6.37                 | 10.85    |
| `90259f9` (05:15 IST, same conditions, CPU idle 98.6% before and after) | 4.30                   | 4.55       | 16.17    | 1.66                 | 5.02     |

## Everyday conditions (M2.11.6)

Not a formal benchmark: this is what the model costs on a laptop in ordinary
use, to set against the quiet, plugged-in figures above. Same 8 threads and
`bench_layers` command as the performance pass, then the full climate run:
`RAYON_NUM_THREADS=8 fms-climate climate --climate-frequency 16 --heightmap data/ETOPO_2022_v1_60s_N90W180_surface.tif`
(Earth, 2,562 cells, six accelerated years, dt 2057 s, 92,030 steps).

**Conditions** (2026-11-10, 18:52 to 19:48 IST): Apple M5 Pro, 15 cores, 24 GB,
macOS 26.6.2, rustc 1.99.1. On battery, normal power mode (`powermode 0`, not
High Power), battery falling from 82% to 31%. The user was working: Steam
(about 60% of a core throughout), Firefox with its GPU and content processes,
the window server at about 45%. CPU 87% idle before the jobs and 52–55% idle
during them (`top`); load average 8–9; 0.5–1 GB of memory unused by `top`'s
count, with 22–23 GB in use. Commit `e2cc65d` plus the ice-sheet work of
this step, which has no effect on the aquaplanet benchmark. The fourth Earth run
has ice sheets; the first ran before the ice-sheet code was built in.

| Job                                                        | Here                                                                                                               | Quiet, plugged in (`90259f9`) |
| ---------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------ | ----------------------------- |
| `bench_layers`, moist aquaplanet, 2,562 cells, median of 3 | 4.52, 4.52 and 4.93 ms/step (first two at the start of each pass, last after 35 minutes as the browser got busier) | 4.30 ms/step                  |
| Earth, six years, inside the time loop                     | 596, 616, 621 and 633 s (155, 149, 148 and 145 steps/s; 6.5–6.9 ms/step)                                           | not measured                  |

- Everyday use cost 5–15% on the like-for-like benchmark. Only 8 of the 15
  cores are in use, which leaves the background work room.
- A year on Earth takes about 1.7 minutes, so a six-year run about 10. Land,
  sea ice and sampling put Earth about 45% above the aquaplanet's ms/step at
  the same grid.
- Whole-process time was 8 s beyond the time loop in one run and 135 s in
  another, with the same outputs. The cause wasn't identified; memory was tight
  and the browser busy at the time.
- Battery fell about 55–60% per hour during these runs, so a laptop on battery
  lasts a little under two hours of continuous simulation with this workload.
- `cargo t2`, the workspace tests, took 4 min 9 s under the same conditions.
