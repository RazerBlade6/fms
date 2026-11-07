# Benchmarks

Criterion suites for the performance-critical building blocks. Run them locally:
CI machines are too noisy for timing.

```sh
cargo bench -p fms-geo --bench geo            # grids, geometry, transfers, point location
cargo bench -p fms-sched --bench sched        # scheduler overhead
cargo bench -p fms-planet --bench planet      # insolation, orbit
```

`-- --quick` gives a fast estimate. To check a change for regressions, save a
baseline first:

```sh
cargo bench -p fms-geo --bench geo -- --save-baseline before
# … make the change …
cargo bench -p fms-geo --bench geo -- --baseline before
```

HTML reports go to `target/criterion/`.

## Baseline

Taken on 2026-11-07 on the development Mac (Apple silicon), release profile,
`--quick`. Sizes: f45 = 20,252 cells (default climate grid); f160 = 256,002 cells
(default CLI terrain grid).

| Benchmark | Time | Notes |
|---|---|---|
| `build/geodesic_f160` | 11.3 ms | Subdivision only |
| `build/centroidal_f45` | 406 ms | Multilevel Lloyd to 1e-6 |
| `build/geometry_f160` | 54.7 ms | Metric terms, TRiSK and Perot weights |
| `build/hilbert_order_f160` | 89.0 ms | Renumbering and rebuilding adjacency |
| `two_level_f160_to_f45/build_map` | 3.96 ms | |
| `two_level_f160_to_f45/restrict_mean` | 155 µs | |
| `two_level_f160_to_f45/restrict_sum` | 154 µs | |
| `two_level_f160_to_f45/prolong` | 214 µs | |
| `operators_f45/perot_reconstruction` | 149 µs | |
| `locate_f160/cells_10k_points` | 636 µs | ~64 ns per point |
| `locate_f160/triangles_10k_points` | 1.39 ms | |
| `planet/insolation_20k` | 128 µs | Parallel |
| `planet/daily_mean_20k` | 3.07 ms | Serial loop in the bench; the system runs it in parallel |
| `planet/orbit_state` | 93 ns | Kepler solve |
| `scheduler/step_20_systems` | 119 µs | 20 trivial systems in 2 stages. This is mostly the cost of starting parallel work per stage; at speed 2 (~340 steps/s) it is ~4% overhead. Candidate for optimization (cache stage plans, run tiny stages serially) if it shows up in M2 profiles |
