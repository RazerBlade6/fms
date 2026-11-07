# Observations

Known limitations, measured quirks and ideas for later: things we noticed while
building and decided **not** to act on yet. Each entry says what we saw, why it
is acceptable for now, and what would make us revisit it. Planned work lives in
[ROADMAP.md](ROADMAP.md); this file is for everything else.

When something here gets scheduled, link the roadmap step. When it is resolved,
set its status to *Resolved* and say how, rather than deleting it.

| # | Observation | Area | Status |
|---|---|---|---|
| O-1 | [Operator errors peak at the 12 pentagons](#o-1-operator-errors-peak-at-the-12-pentagons) | Numerics | Scheduled: M2 step 14 |
| O-2 | [Cell sizes vary by up to 1.76×](#o-2-cell-sizes-vary-by-up-to-176) | Grid | Watching |
| O-3 | [Grid frequencies need small prime factors](#o-3-grid-frequencies-need-small-prime-factors) | Grid | Accepted |
| O-4 | [Terrain cells are assigned whole to climate cells](#o-4-terrain-cells-are-assigned-whole-to-climate-cells) | Grid | Accepted |
| O-5 | [Kepler's third law is not enforced](#o-5-keplers-third-law-is-not-enforced) | Planet | Accepted |
| O-6 | [The equation of time is ignored](#o-6-the-equation-of-time-is-ignored) | Planet | Accepted |
| O-7 | [Daily-mean insolation uses one orbital state per day](#o-7-daily-mean-insolation-uses-one-orbital-state-per-day) | Planet | Accepted |
| O-8 | [Scheduler overhead is ~119 µs per step](#o-8-scheduler-overhead-is-119-µs-per-step) | Performance | Scheduled: M2 step 16 |
| O-9 | [Vector fields render as arrows only](#o-9-vector-fields-render-as-arrows-only) | Tooling | Watching |
| O-10 | [Field names must be `'static` strings](#o-10-field-names-must-be-static-strings) | Fields | Scheduled: M4 step 2 |
| O-11 | [Land below sea level counts as ocean](#o-11-land-below-sea-level-counts-as-ocean) | Terrain | Watching |

---

### O-1: Operator errors peak at the 12 pentagons

**Seen (M1.3, M1.5):** three operators have errors that converge in RMS
(halving when resolution doubles) but whose *maximum* stays fixed, located at or
next to the 12 pentagons:

| Operator | Max error (k=16 → k=32) | RMS error (k=16 → k=32) |
|---|---|---|
| TRiSK tangential reconstruction | 0.62% → 0.71% | 0.094% → 0.039% |
| Perot cell-vector reconstruction | 2.0% → 2.0% | 0.30% → 0.15% |
| Vertex vorticity (Rossby–Haurwitz) | 6.3% → 5.8% of peak | 1.4% → 0.7% of peak |

This matches the literature (Peixoto 2016, *Accuracy analysis of mimetic finite
volume operators on geodesic grids*).

**Why it matters:** fixed local errors at the pentagons can leave a spurious
wavenumber-5 / 12-point pattern in long GCM runs ("grid imprinting"), e.g. in
jets and storm tracks.

**Plan:** M2 step 14 measures imprinting in the Held–Suarez run. If it is
visible: Heikes–Randall grid optimization, or Peixoto's corrected
reconstruction weights.

### O-2: Cell sizes vary by up to 1.76×

**Seen (M1.2):** after Lloyd optimization the largest cell is 1.46× the
smallest at 20k cells, and 1.76× at 1.6M. The smallest cells cluster around the
pentagons.

**Why it matters:** the smallest cells set the CFL timestep, costing ~17% of the
maximum `dt` on the climate grid.

**Revisit if:** the M2 performance pass shows `dt` is the bottleneck. Options:
density-weighted Lloyd, or spring-dynamics grid optimization.

### O-3: Grid frequencies need small prime factors

**Seen (M1.2):** multilevel Lloyd refines by prime factors ≤ 7, so frequencies
like 45 or 400 build in seconds. A frequency with a large prime factor (e.g.
401) falls back to plain Lloyd, which needs O(frequency²) iterations.

**Accepted:** we choose the frequencies. The UI should only offer good ones.

### O-4: Terrain cells are assigned whole to climate cells

**Seen (M1.4):** each terrain cell belongs entirely to the climate cell
containing its centre, rather than being split by exact polygon overlap. A climate
cell's child area is within ±25% of its own area (72–83 children at 1.6M → 20k).

**Accepted:** sums are still exactly conservative globally, and area-weighted
means are accurate for smooth fields. Exact overlap remapping would only matter
if per-cell conservation at the climate scale is ever needed.

### O-5: Kepler's third law is not enforced

**Seen (M1.7):** year length and semi-major axis are set independently, so
physically impossible combinations are allowed (`P² ∝ a³ / M★` is not checked).

**Accepted** for fantasy worlds. If realism is wanted: an optional star-mass
parameter that derives one from the other, or flags inconsistent combinations.

### O-6: The equation of time is ignored

**Seen (M1.7):** the Sun's position uses mean solar time. On an eccentric or
tilted orbit, true local noon drifts from mean noon by minutes (about ±16 min on
Earth).

**Accepted:** irrelevant for climate. Revisit only if exact sundials or sunrise
times are ever displayed.

### O-7: Daily-mean insolation uses one orbital state per day

**Seen (M1.8):** `planet.daily_insolation` uses the orbital state at mid-day
for the whole day.

**Accepted:** declination and distance change negligibly within a day for any
year longer than ~20 days (the minimum allowed).

### O-8: Scheduler overhead is ~119 µs per step

**Seen (M1.5.4):** 20 trivial systems in 2 stages cost ~119 µs per step, mostly
from starting parallel work for each stage. At speed 2 (~340 steps/s) that's ~4%
overhead.

**Plan:** look at it in the M2 performance pass (step 16) if profiles show it.
Options: cache stage plans; run stages of tiny systems serially.

### O-9: Vector fields render as arrows only

**Seen (M1.5.2):** the CLI draws vector fields as speed plus arrows. Streamlines
read better for large-scale circulation (Hadley cells, gyres).

**Revisit if:** M2 images of the general circulation are hard to read.

### O-10: Field names must be `'static` strings

**Seen (M1.5):** fine for every system compiled in, but names loaded at runtime
(unknown fields in a save file, scripting) can't be registered.

**Plan:** switch `FieldSpec::name` to `Cow<'static, str>` in M4 step 2 (save
format), where keeping unknown fields first needs it.

### O-11: Land below sea level counts as ocean

**Seen (M2.1):** land and sea are split at elevation 0, so inland basins below
sea level count as ocean in `terrain.land_fraction` and get a flat sea-level
surface. On ETOPO the Caspian Sea (−28 m) is the only one visible at climate
resolution. The Dead Sea, Qattara Depression and Death Valley are a few terrain
cells at most.

**Why it's acceptable for now:** the Caspian really is water, so for climate it
mostly behaves right, and nothing reads the land/sea split yet.

**Revisit when:** the slab ocean (M2 step 11) needs a real ocean mask: take the
ocean to be the below-sea-level region connected to the open sea, and leave
closed basins to hydrology's lakes and depressions.
