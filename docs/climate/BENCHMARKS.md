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
vertical mode. The script's header lists the options (`WORLD=aqua`, `THREADS`
and `IDLE=spin|sleep` to override the run's own choice of threads, `COLD=1`
for a run from rest, `STEPS=N` to stop after N steps with the timing report
alone, `THROTTLE=1`).

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
made. The spread between repeats was about 0.5% in the morning's runs and
several percent later in the day ("Drift on the day" below).

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

### Thread scheduling and allocation in the dynamics (15.3, 15.4)

Bit-identical to the baseline. Steps a second at 2,562 cells:

| Commit | 15 threads | 8 threads | 4 threads |
| --- | --- | --- | --- |
| Baseline `22b0240` | 176 | 189 | 162 |
| 15.3, thread scheduling `83fbd5b` | 209 | 214 | |
| 15.4, no allocation in a step `2d60891` | 208 | 216 | 176 |

- 15.3: each phase of the dynamics is one parallel region with the layers
  (or the vertical modes) as its tasks, each running its layer's whole chain
  of operators; the couplings between layers are plain loops. A loop over one
  layer's cells or edges goes to the pool only when it is long enough to be
  split (8,192 elements), and otherwise runs where it is, waking no worker.
  The implicit right-hand side fell from 0.53 to 0.24 ms a step, the solves
  from 0.75 to 0.54, back-substitution from 0.22 to 0.14, humidity transport
  from 0.98 to 0.75 and the filter from 0.14 to 0.10.
- A parallel region still costs 50 to 70 µs on the critical path at 15
  threads, mostly waking sleeping workers, and the kernel's share of the CPU
  time is still 70%: the regions that remain are the floor of this design.
- 15.4: every array the dynamics used to allocate each step is owned by the
  stepper, the model, the transport or the solver. Worth 1 to 2%: allocation
  was never the cost at this grid (O-29).

### Fusion, solver, column physics and sampling (15.5 to 15.8)

Steps a second at 2,562 cells, each row on top of the one above. Every row
but the stencil gives results identical to the baseline, to the last byte of
every output, at every thread count tried.

| Commit | 15 threads | 8 threads | 10 threads |
| --- | --- | --- | --- |
| 15.3 and 15.4 `2d60891` | 208 | 216 | |
| 15.7 column physics and radiation, 15.8 sampling `2a475dd` | 285 | 306 | |
| 15.7 kinetic heating reads the moist physics' Exner functions `f2a56e5` | 293 | | 311 |
| 15.5 fused loops and merged regions `9b9ecb2` | 306 | | |
| 15.6 Helmholtz operator from a stencil `0a6ae45` | 302 to 328 | 342 | 345 to 347 |

One thread, with 15.7 and 15.8 alone: 53 to 79 steps a second.

- **15.7, column physics and radiation.** The all-sky and clear-sky fluxes
  are computed in one pass that shares the profile, the Planck emission of
  every level, the ozone transmissivities, and each slab whose optics are the
  same in both skies. Nothing is allocated per column; each of moisture,
  surface fluxes, radiation and the column tendencies is one parallel region
  (chunks of 8 columns, 128 for the column tendencies, measured against 16 to
  4,096); convective parcels are lifted four columns side by side; the
  surface exchange shares its logarithms. On one thread moisture fell from
  5.02 to 2.12 ms a step, a radiation call from 21.6 to 9.7 ms, surface
  fluxes from 1.99 to 1.66 ms. The kinetic heating reused the Exner functions
  the moist physics had just computed (column tendencies 0.36 to 0.23 ms).
- **15.8, sampling.** The zonal climatology computed each cell's Exner
  functions 16 times a sample, and the run built its per-cell fields twice:
  0.20 to 0.05 ms a step. The ocean's per-step code was left: 0.13 ms a step
  at 2,562 cells, 60% of it the gyres' solve.
- **15.5, fusion (O-12).** Operators have per-element forms, so one pass
  computes what several did: the dynamics' second layer region is 6 passes
  instead of about 16, a conjugate-gradient iteration is 4 passes instead of
  9, the back-substitution's region also does the mass fluxes, and the
  transport's last region also does the time filter. Mass fluxes and filter
  together fell from 0.23 to 0.01 ms a step.
- **15.6, solver (O-13).** The Helmholtz operator is applied from a per-cell
  stencil (a diagonal and one coupling per edge, the same number in both
  rows, so the matrix is exactly symmetric) instead of a gradient and a
  divergence: the solves fell from 0.51 to 0.29 ms a step on 8 threads.
  Iterations are unchanged (22.0, 11.0, 5.0 and 3.8 a solve by mode, 10.5 on
  average). This is the one change of the pass that alters results, at
  rounding level: the tolerance (1e-8), the preconditioner and the first
  guess are as before, and the golden hash was updated.
- The sum of the kinetic-energy dissipation used the pool's own order of
  additions, which depends on the thread count from 8,192 edges; it now adds
  in a fixed order. It feeds only the energy budget's printed line.

The stencil's year against the baseline's, and against a third year whose
solver differed within the same tolerance (an extrapolated first guess, not
kept): the differences are those of any two years.

| Final-year value | Baseline | Stencil | Third trajectory |
| --- | --- | --- | --- |
| Surface temperature (K) | 288.06 | 288.02 | 288.04 |
| Outgoing longwave (W m⁻²) | 228.2 | 228.1 | 228.1 |
| Precipitation (mm/day) | 3.44 | 3.45 | 3.45 |
| Northern Hadley cell (10⁹ kg/s) | 40.8 | 47.3 | 44.7 |
| Southern Hadley cell (10⁹ kg/s) | −59.8 | −56.7 | −52.1 |
| Northern Hadley edge | 26.3° | 26.1° | 26.2° |
| Southern jet (m/s) | 28.7 at −35.3° | 29.3 at −35.8° | 29.2 at −35.3° |
| Equator − pole (K) | 44.3 | 44.4 | 44.6 |
| Rain belt's seasonal swing | 7.6° | 8.0° | 7.7° |

Wall time by phase after 15.8, 2,562 cells on 10 threads (run `auto-2562-1`,
348 steps a second), against the baseline's best (8 threads), ms a step:

| Phase | Baseline, 8 threads | After, 10 threads |
| --- | --- | --- |
| Tendencies (physics, dynamics and hyperdiffusion side by side) | 1.58 | 0.83 |
| Implicit right-hand side | 0.35 | 0.19 |
| Helmholtz solves | 0.72 | 0.29 |
| Back-substitution, mass fluxes | 0.23 | 0.16 |
| Flux-corrected transport, time filter | 0.81 | 0.53 |
| Surface fluxes | 0.39 | 0.26 |
| Radiation | 0.74 | 0.25 |
| Ocean, sea ice, land | 0.18 | 0.18 |
| Health check, sampling, state copy | 0.23 | 0.13 |
| **A step** | **5.29** | **2.87** |

### Threads and idle workers (O-54)

The pool's workers yield and then sleep when they have no work, and a step is
a sequence of short parallel regions: every region wakes them and every gap
puts them back. That is the kernel time and the context switches of O-54.
Keeping workers spinning instead (looking for work in a loop, never
sleeping) was measured at both grids, results identical in every case.

2,562 cells, steps a second, after 15.6:

| Threads | 6 | 8 | 10 | 12 | 15 |
| --- | --- | --- | --- | --- | --- |
| Workers sleep | 328 | 342 | 345 | 342 | 326 |
| Kernel time (s) | 24 | 50 | 79 | 106 | 160 |
| Spinning workers (how many) | | 350 (4) | 348 (5) | | 355 (5), 349 (14) |

20,252 cells, the first 3,000 steps from rest, ms a step (repeats in the same
state of the machine; see the note on drift below):

| Threads | 1 | 5 | 8 | 10 | 12 | 15 |
| --- | --- | --- | --- | --- | --- | --- |
| Workers sleep | 75.4 | 20.5 | 17.0 | 16.5 | 16.1 | 16.2 to 17.4 |
| Context switches (millions) | 0 | 1.2 | 2.3 | 4.4 | 6.8 | 10 to 14 |
| All but one spinning | | 21.0 | | 15.3 | | 13.1 to 14.9 |
| Context switches (millions) | | 0.9 | | 1.4 | | 2.6 to 4.1 |

- At 20,252 cells spinning is 15 to 24% faster on 15 threads, in runs
  alternated with plain ones, for 5 to 14% more CPU time in all: the sleeping
  pool was already spending a third of its CPU time in the kernel.
- At 2,562 cells it gains 1 to 3% over the best plain thread count for twice
  the CPU time. The dynamics there has only its four layers to share out, so
  most spinners find nothing.
- Threads beyond 10 lose time at 2,562 cells; at 20,252 cells 10 and 15
  sleeping threads are equally fast, and 5 threads (the Super cores alone)
  reach 80% of 15.
- A spinner on every thread is slower than all but one (15.0 against 13.6 ms):
  the thread that runs the loop is one more than the cores. Twelve of 15 is
  the same as 14 within the noise. A pause of 200 spins between looks for
  work is no better than none, and 2,000 is worse.
- With every core spinning, a step's serial stretches run 10 to 25% slower
  (ocean transport 1.16 to 1.42 ms at 20,252 cells): they land on whichever
  thread took the step, more often a Performance core than a Super one.
  Threads can't be pinned without `unsafe` code or a new dependency.
- Running the whole time loop inside the pool changes nothing (16.2 and
  17.0 ms on 15 and 10 threads).
- Splitting the dynamics' loops finer than 4,096 elements, with spinning, is
  worth about 2% at 20,252 cells and was left.

**Adopted** (`a6a04ba`): the run sizes its own pool. From 8,192 cells, where
a layer's loops are split between threads, it takes every core and keeps all
but one spinning; below, it takes a thread for every 256 cells (10 at 2,562)
and lets idle ones sleep. `--threads` and `--idle spin|sleep` override it,
and fewer threads than cores turn the spinning off unless it is asked for.

### The default grid (20,252 cells)

The same run as the baseline's: 1.1 years from rest, 49,818 steps.

| | Steps a second | ms a step | Share of speed 1 | User time (s) | Kernel time (s) | Context switches (millions) |
| --- | --- | --- | --- | --- | --- | --- |
| Baseline `22b0240`, 15 threads | 36.6 | 27.34 | 0.30 | 9,965 | 5,286 | 323 |
| After 15.6 `47981a9`, 15 threads | 55.4 | 18.05 | 0.45 | 6,470 | 3,170 | 199 |
| After 15.6, 10 threads | 56.2 | 17.81 | 0.45 | 5,772 | 1,350 | 76 |
| 14 of 15 threads spinning | 65.3 | 15.33 | 0.53 | 10,521 | 432 | 34 |

The spinning run's results are identical to the plain run's.

Wall time by phase, ms a step, in those runs:

| Phase | Baseline | After 15.6 | Spinning |
| --- | --- | --- | --- |
| Tendencies | 8.10 | 5.39 | 4.92 |
| Implicit right-hand side | 1.70 | 1.16 | 0.85 |
| Helmholtz solves | 6.41 | 2.78 | 1.69 |
| Back-substitution, mass fluxes | 0.96 | 0.99 | 0.67 |
| Flux-corrected transport, time filter | 4.03 | 2.82 | 2.15 |
| Surface fluxes | 1.66 | 1.21 | 1.28 |
| Radiation (every 15th step) | 1.23 | 0.47 | 0.53 |
| Ocean transport (every 15th step, one thread) | 1.03 | 1.16 | 1.42 |
| Ocean surface, land, budgets, humidity fill | 1.30 | 1.32 | 0.91 |
| Health check, sampling, state copy | 0.91 | 0.74 | 0.87 |
| **A step** | **27.34** | **18.05** | **15.33** |

On one thread a step is 75 ms (first 3,000 steps), so 15 threads make it 4.6
times faster asleep and 5.4 spinning: tendencies 7 times, surface fluxes and
radiation 10 times, the implicit right-hand side, the solves, the
back-substitution and the transport 3 to 5 times.

### Drift on the day

From the afternoon of 2026-11-17, after about three hours of all-core load,
the same run was about 7% slower than in the morning and repeats differed by
several percent (15 threads at 20,252 cells: 16.2, then 16.7 to 17.4 ms; the
2,562-cell run on 15 threads: 302 straight after the test suite, 321 to 328
otherwise). `pmset -g therm` recorded no warning. Comparisons above are
between runs alternated in the same state. **The close-out timings (15.9)
wait for a cold machine.**

### Close-out, still to run (15.9)

On a cold, idle machine, with the conditions logged as above:

```sh
scripts/perf.sh conditions
scripts/perf.sh bench 2562 final-2562-1      # three times: -1, -2, -3; the median
scripts/perf.sh sweep 2562 1 4 8 10 15       # the thread counts, each asked for
THREADS=15 IDLE=spin scripts/perf.sh bench 2562 final-2562-spin
COLD=1 YEARS=1.1 scripts/perf.sh bench 20252 final-20252      # like the baseline's run
COLD=1 YEARS=1.1 IDLE=sleep scripts/perf.sh bench 20252 final-20252-sleep
STEPS=3000 COLD=1 YEARS=1.1 THREADS=1 scripts/perf.sh bench 20252 final-20252-t1
scripts/perf.sh profile 2562 final-profile-2562 && scripts/perf.sh top final-profile-2562
scripts/perf.sh same final-2562-1 auto-2562-1     # identical results to today's
```

Then: the final table here (2,562 and 20,252 cells against the baseline), the
target at 20,252 cells decided, the roadmap's 15.9 ticked. The runs of
2026-11-17 kept in `output/m2.15/` for comparison: `base-2562-1` (baseline
results), `auto-2562-1` (results after 15.6), `base-20252-cold`,
`final-20252-cold` and `spin-20252-cold`.

### Measured and not adopted

- An extrapolated first guess for the solver (twice the last solution minus
  the one before): iterations rose in every mode (642,056 to 648,689 a year
  in all). The last two solutions are of alternating leapfrog levels.
- Pipelining the physics a step behind the dynamics, a longer radiation
  interval, a thinner clear-sky calculation, coarser parcel steps, tabulated
  saturation, Planck and `pow`, a looser solver tolerance: each changes the
  model's results for a few percent, and speed 2 didn't need them.
- Computing the next step's dynamics tendencies beside the physics' advance:
  the advance is shorter than those tendencies on four steps in five, so
  nothing is saved.

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
