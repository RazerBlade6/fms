# Development

How we work on this project: setup, conventions, and the rules that keep the
simulation correct and reproducible.

## Setup

The toolchain is pinned in `rust-toolchain.toml` (Rust 1.99.1 with rustfmt and
clippy); `rustup` installs it automatically on first build.

```sh
cargo test --workspace                                 # all tests (~4 minutes once built; the climate and Held–Suarez runs are most of it)
cargo test -p fms-atmos --release -- --ignored         # the slow tests (multi-year climate runs, minutes)
cargo clippy --workspace --all-targets -- -D warnings  # lints, as CI runs them
cargo fmt --all                                        # format
RUSTDOCFLAGS="-D warnings" cargo doc --workspace --no-deps   # API docs, as CI builds them
cargo run --release -p fms-climate -- --help               # the headless runner
```

Use `--release` for anything that builds full-size grids. Debug builds use
`opt-level = 1` so tests stay reasonably fast, but grid construction is still
~10× slower than release. The numerical crates `fms-grid` and `fms-atmos` build
at `opt-level = 3` even in debug, because their tests run multi-day simulations
(debug info is kept; drop the override in `Cargo.toml` if you need to step
through optimized-out code).

## Repository layout

| Path | Contents |
|---|---|
| `crates/` | The Cargo workspace. Each crate is described in [PLAN.md](PLAN.md) §3.3 |
| `docs/PLAN.md` | Design: architecture, decisions, subsystems, equations |
| `docs/ROADMAP.md` | Milestones and steps, ticked as they're done |
| `docs/FMSMAP.md` | The `.fmsmap` file every `fms` module exchanges |
| `docs/OBSERVATIONS.md` | Known limitations and ideas we've chosen not to act on yet |
| `docs/CLI.md` | The `fms-climate` command-line runner |
| `docs/BENCHMARKS.md` | How to run benchmarks; baseline numbers |
| `output/` | Generated images and other artifacts for review (gitignored) |

The README stays short and points here.

## Git workflow

- **`dev`** is the working branch; all work happens there. **`main`** is the
  release branch and only ever moves by fast-forwarding to `dev`.
- Commit to `dev` after each completed roadmap step, with a message naming the
  step (e.g. `M1.4: point location and fine/coarse grid mapping`).
- Pushing needs an explicit go-ahead.
- "Ship it" / "push to prod" means: commit, push `dev`, fast-forward `main` to
  `dev`, push `main`. `main` is never force-pushed.
- CI runs on every push to `dev`: format, clippy, API docs, and the golden-hash
  tests on macOS, Windows and Linux. The full test suite runs locally before a
  push, not in CI. `main` only ever receives what has already passed on `dev`.
  From M5, pushes to `main` build a GitHub Release (CD).

## Milestones

- Every roadmap step ends with an automated test or a visual check (an image in
  `output/`). Tick it in `ROADMAP.md` (`- [ ]` → `- [x]`) when it's done.
- After the last step of a milestone, **stop and review**: what was built, what
  is complete and verified, and what's missing that later milestones don't
  already cover. The review's findings either become roadmap steps or
  entries in `OBSERVATIONS.md`.

## Determinism rules

The same save on the same build must give **bit-identical** results on any
thread count and on macOS, Windows and Linux. Three rules make that true:

1. **No platform math.** `sin`, `cos`, `atan2`, `exp`, `ln`, `pow`, `hypot` and
   friends come from `fms_math` (a pure-Rust `libm`), never from `f64` methods.
   `clippy.toml` forbids the std versions, so CI catches mistakes. Basic
   arithmetic, `sqrt`, `abs`, `floor`, `round`, `min`/`max` and `powi` are exactly
   specified and fine to use. If you need a function `fms_math` lacks, add a
   wrapper there.
2. **Fixed-order reductions.** Parallel sums must not depend on how work is split
   across threads. Reduce per element in a fixed order (as the two-level
   transfers do), or use associative operations such as `min` and `max`.
3. **Declared access only.** Systems touch fields only through their
   `StepContext`, which checks every access against what they declared. That's
   what lets the scheduler run systems in parallel without changing results.

### Golden-hash tests

Two tests pin exact results and run on all three CI platforms:

| Test | Pins |
|---|---|
| `fms-math` `bit_exact_golden_values` | 7,000 `libm` results |
| `fms-sim` `results_are_bit_identical_across_platforms` | The optimized climate grid and both insolation fields after 100 steps |

If one fails, first ask whether the change was **meant** to alter numerical
results:

- **No:** something became non-deterministic or platform-dependent. Find it;
  don't update the hash.
- **Yes** (e.g. a deliberate change to the grid optimizer): the failure message
  prints the new values. Paste them into the test's `GOLDEN` constant, say so in
  the commit message, and let CI confirm all three platforms agree.

## Conventions

- **Names.** Fields and systems are `namespace.name` in snake_case, with the
  owning crate or system as the namespace (`terrain.elevation`,
  `atmos.dynamics`). The registry and scheduler reject anything else.
- **Units.** SI throughout (m, s, kg, K, Pa, W m⁻²). Geometry in `fms-grid` is on
  the unit sphere; scale lengths by the planet radius and areas by its square.
  Angles are in radians in code; degrees only in user-facing parameters, with a
  `_deg` suffix.
- **Coordinates.** +z is the north pole and +x is (0°, 0°). `t = 0` is the
  northern vernal equinox, at midnight at longitude 0.
- **Errors.** Library crates return `Result` for anything a user can cause (bad
  parameters, unknown field names). They panic only on internal invariants and
  programming errors, with a message that names the problem.
- **Lints.** `unsafe_code` is denied; missing docs warn (and fail CI). Fix
  warnings rather than silencing them; if an `allow` is truly needed, add a
  comment saying why.

## Testing

- **Test against known answers.** Prefer analytic solutions and standard
  benchmarks (Williamson test cases, Held–Suarez, Earth's real values) over
  "looks plausible".
- **Check convergence, not just accuracy.** For numerical operators, test that
  the error shrinks with resolution at the expected rate. If a norm doesn't
  converge for a known reason, assert what *does* hold (e.g. RMS converges, and
  the maximum error stays bounded at the pentagons), document why, and record it
  in `OBSERVATIONS.md`.
- **Look at the output.** Anything spatial gets rendered with `fms-climate` and
  inspected; several real bugs (the coastline colour, an off-by-one date stamp)
  were caught this way, not by tests.
- **Tolerances need headroom.** Set bounds with margin over the measured value,
  and record the measured value in a comment or in the test's output.

## Performance work

Measure before and after with the criterion suites (see
[BENCHMARKS.md](BENCHMARKS.md)):

```sh
cargo bench -p fms-grid --bench geo -- --save-baseline before
# … change …
cargo bench -p fms-grid --bench geo -- --baseline before
```

For whole-simulation costs, `fms-climate dump` and `fms-climate animate` print a per-system
timing breakdown after each run. Update the baseline table in `BENCHMARKS.md`
when a change moves the numbers on purpose.
