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
scripts/perf.sh same <run> <reference>            # are two runs' results identical?
scripts/perf.sh top profile-earth-2562            # the profile as text: CPU time by function
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

## Performance pass (M2.15)

The target is the full climate run, `fms-climate run`, on Earth at 2,562
cells: speed 2, one model week a second, is 294 steps a second at its 2,057 s
step. That is the floor, not the ceiling. The default 20,252 cells (697 s
step) is timed at the start and the end; speed 1, a model day a second, is
124 steps a second there and speed 2 is 868.

**The benchmark** is one model year from a restart file saved after the
default nine-year spin-up (`scripts/perf.sh bench 2562`: 15,340 steps, with
the year's statistics, monthly means and outputs as in any run). At 20,252
cells it is 1.1 years from rest without ocean jumps
(`COLD=1 YEARS=1.1 scripts/perf.sh bench 20252`: 49,818 steps). Steps a
second are for the time loop; the median of three runs where three were
made, and the spread between repeats is about 0.5%.

**Conditions** (2026-11-17): MacBook Pro, Apple M5 Pro, 15 cores (5 "Super",
10 "Performance"), 24 GB; macOS 26.6.2 (25G83); rustc 1.99.1; release
profile. On AC power in High Power mode (`pmset` powermode 2), battery 95%.
Only the Ghostty terminal, Activity Monitor and the agent driving the runs
were open; CPU 96–98% idle before and after every run (`top`), 8–10 GB of
memory unused. Runs went one at a time.

### Baseline (15.2, commit `22b0240`)

| Threads | 1 | 2 | 4 | 6 | 8 | 10 | 12 | 15 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Steps a second, 2,562 cells | 53 | 94 | 162 | 179 | 189 | 189 | 186 | 176 |
| User time (s) | 306 | 322 | 331 | 361 | 391 | 414 | 436 | 470 |
| Kernel time (s) | 1 | 2 | 8 | 42 | 90 | 143 | 200 | 302 |

Results are identical at every thread count, to the last byte of every
output. The run is 0.60 of speed 2 on all cores and fastest on 8 to 10. A
step is 18.9 ms of work on one thread and takes 5.7 ms on 15: 3.3 times
faster for 2.5 times the CPU time, most of that in the kernel.

Wall time by phase, in ms a step (the physics and the dynamics tendencies run
side by side, so the tendencies take the longer of the two):

| Phase | 2,562 cells, 1 thread | 8 threads | 15 threads | 20,252 cells, 15 threads |
| --- | --- | --- | --- | --- |
| Physics tendency: moisture | 5.02 | 1.03 | 0.84 | 4.88 |
| Physics tendency: column tendencies | 0.57 | 0.40 | 0.53 | 1.97 |
| Dynamics tendency | 1.08 | 1.08 | 0.95 | 5.60 |
| Hyperdiffusion | 1.09 | 0.39 | 0.55 | 2.38 |
| Implicit right-hand side | 0.57 | 0.35 | 0.53 | 1.70 |
| Helmholtz solves (with the interface ratios) | 2.03 | 0.72 | 0.75 | 6.41 |
| Back-substitution | 0.19 | 0.14 | 0.22 | 0.53 |
| Humidity transport | 1.26 | 0.79 | 0.98 | 3.23 |
| Time filter and fixers | 0.10 | 0.13 | 0.14 | 1.48 |
| Surface fluxes | 1.99 | 0.39 | 0.42 | 1.66 |
| Radiation (every 5th step; every 15th at 20,252) | 4.32 | 0.74 | 0.54 | 1.23 |
| Ocean, sea ice, land | 0.16 | 0.18 | 0.18 | 1.52 |
| Health check, sampling, state copy | 0.20 | 0.23 | 0.28 | 0.91 |
| **A step** | **18.87** | **5.29** | **5.69** | **27.34** |
| **Steps a second** | **53** | **189** | **176** | **37** |

- The column physics is 63% of the work (11.9 ms on one thread) and the part
  that gains from threads: 5 to 8 times on 15.
- The dynamics gains almost nothing: 6.3 ms on one thread, 4.1 on 15. At this
  grid every array is below the smallest piece handed to a thread, so the
  only parallelism is across the four layers, in dozens of regions a step
  with tens of microseconds of work each.
- The profile at 15 threads (`scripts/perf.sh top`) has 72% of the CPU time
  in the kernel: the thread pool's idle workers yielding (33%), sleeping and
  being woken (36%), and waiting on the pool's own lock (3%). That is where
  the kernel time of O-54 goes. Allocation is 1.5%, and page faults don't
  grow with the threads.
- Of the time in the model itself, radiation's `column_fluxes` is 28%
  (two calls a column: all-sky and clear-sky), `pow` 10% and `exp` 6%.
- The solver takes 22.0 iterations a solve in the external mode, then 11.0,
  5.0 and 3.8 (10.5 on average): the external mode alone, 0.69 ms on one
  thread, sets the solves' wall time.
- Setup is 7 s (building the world from the heightmap) and output 0.6 s.
- At 20,252 cells the run makes 37 steps a second, 0.30 of speed 1, with 35%
  of its CPU time in the kernel.

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
