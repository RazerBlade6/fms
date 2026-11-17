# Benchmarks

Timings of the climate run, and how to take them. The performance pass of
M2.15 is the reference: every later timing is compared against its tables.
Earlier benchmarks are summarized at the end.

## Conditions

Every benchmark is run one configuration at a time, never alongside tests or
builds, with repeats and the median reported, and logged here with:
the machine, power source and mode, what else is running, CPU idle and free
memory at the start (`top -l 2`) and end, macOS and rustc versions, the
commits compared and the date.

## Timing and profiling the climate run

`scripts/perf.sh` runs `fms-climate run` under `/usr/bin/time -l` or the
[samply](https://github.com/mstange/samply) profiler
(`cargo install --locked samply`), from a restart file it makes once per
resolution:

```sh
scripts/perf.sh spinup 2562       # a full Earth run from rest; saves the restart file
scripts/perf.sh bench 2562        # one timed year from it
scripts/perf.sh sweep 2562        # the same at 1, 4, 8 and 15 threads
scripts/perf.sh profile 2562      # the same under samply
scripts/perf.sh view profile-earth-2562
scripts/perf.sh phases bench-earth-2562           # a run's time loop by phase
scripts/perf.sh compare <before> <after>          # two runs' phases, in ms a step
scripts/perf.sh conditions        # the machine's state, as logged below
```

Each run gets a folder in `output/m2.15/` with its log and the machine's
state before and after; `RUNS.txt` there has every command and `results.tsv`
a row per run: steps a second in the time loop, CG iterations, and the whole
process's user and kernel time, context switches, page faults and peak
memory. `phases` and `compare` read the timing report every run writes
([CLI.md](CLI.md)): wall time by phase of a step, and CG iterations by
vertical mode. The script's header lists the options (`WORLD=aqua`, `THREADS`,
`COLD=1` for a run from rest, `THROTTLE=1`).

## Earlier benchmarks

Superseded by the performance pass; kept for comparison. Same machine (Apple
M5 Pro, 15 cores, 24 GB), quiet and on AC power in High Power mode unless
said otherwise.

**Building blocks** (criterion, 2026-11-07, `--quick`).
`cargo bench -p fms-grid --bench geo`, `-p fms-sched --bench sched` and
`-p fms-planet --bench planet`; `-- --save-baseline before` and
`-- --baseline before` compare a change. f45 is 20,252 cells, f160 256,002.

| Benchmark                                         | Time                      |
| ------------------------------------------------- | ------------------------- |
| `build/geodesic_f160`, `geometry_f160`            | 11.3 ms, 54.7 ms          |
| `build/centroidal_f45`, `hilbert_order_f160`      | 406 ms, 89.0 ms           |
| `two_level_f160_to_f45/build_map`                 | 3.96 ms                   |
| `two_level_f160_to_f45` restrict (mean, sum), prolong | 155 µs, 154 µs, 214 µs |
| `operators_f45/perot_reconstruction`              | 149 µs                    |
| `locate_f160` cells, triangles (10k points)       | 636 µs, 1.39 ms           |
| `planet/insolation_20k`, `daily_mean_20k`, `orbit_state` | 128 µs, 3.07 ms, 93 ns |
| `scheduler/step_20_systems`                       | 119 µs (O-8)              |

**Layers** (M2.9.5, 2026-11-08, `bench_layers` example, 2,562 cells, 15
threads, median of 3). Moist aquaplanet: 7.41 ms a step on two layers, 9.46
on three, 13.66 on five; dry Held–Suarez: 3.14, 4.33 and 6.85. The layered
core cost 12% over the two-layer one it replaced (O-29).

**Early performance pass** (M2.11.1.5, 2026-11-10, `bench_layers`, moist
aquaplanet on four layers, median of 3), in ms a step:

| Commit               | 2,562 cells, 8 threads | 15 threads | 1 thread | 642 cells, 4 threads | 1 thread |
| -------------------- | ---------------------- | ---------- | -------- | -------------------- | -------- |
| `775c84a` (baseline) | 11.19                  | 11.63      | 35.53    | 6.37                 | 10.85    |
| `90259f9`            | 4.30                   | 4.55       | 16.17    | 1.66                 | 5.02     |

**Everyday conditions** (M2.11.6, 2026-11-10; on battery, normal power mode,
the laptop in use at 52–55% idle). `bench_layers` as above took 4.52–4.93 ms
a step against 4.30 quiet. The full Earth run at 2,562 cells on 8 threads
(six accelerated years, 92,030 steps) ran at 145–155 steps/s. Whole-process
time was 8 s beyond the time loop in one run and 135 s in another, cause
unknown. Battery fell 55–60% an hour; `cargo t2` took 4 min 9 s.
