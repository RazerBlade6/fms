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
| O-12 | [Shallow-water steps are overhead-bound](#o-12-shallow-water-steps-are-overhead-bound) | Performance | Scheduled: M2 step 16 |
| O-13 | [Semi-implicit steps are long but expensive](#o-13-semi-implicit-steps-are-long-but-expensive) | Performance | Scheduled: M2 step 16 |
| O-14 | [The undamped core collapses to the grid scale](#o-14-the-undamped-core-collapses-to-the-grid-scale) | Dynamics | Resolved: M2 step 6 |
| O-15 | [Baroclinic growth is ~75% of local theory](#o-15-baroclinic-growth-is-75-of-local-theory) | Dynamics | Accepted |
| O-16 | [No sponge layer](#o-16-no-sponge-layer) | Dynamics | Accepted |
| O-17 | [Zig-zag vertex modes limit explicit hyperdiffusion](#o-17-zig-zag-vertex-modes-limit-explicit-hyperdiffusion) | Numerics | Accepted |
| O-18 | [Two layers have no surface wind](#o-18-two-layers-have-no-surface-wind) | Dynamics | Resolved: M2.10.1 |
| O-19 | [Radiation is calibrated, not resolved](#o-19-radiation-is-calibrated-not-resolved) | Radiation | Accepted |
| O-20 | [No cloud or ice feedbacks yet](#o-20-no-cloud-or-ice-feedbacks-yet) | Radiation | Scheduled: M2 steps 10–11 (clouds done in step 9) |
| O-21 | [Dissipated kinetic energy leaks out of the planet](#o-21-dissipated-kinetic-energy-leaks-out-of-the-planet) | Energy | Resolved: M2.10.1 |
| O-22 | [The subtropics rain out their own evaporation](#o-22-the-subtropics-rain-out-their-own-evaporation) | Moisture | Improved: M2.10.1; rest in M2 step 15 |
| O-23 | [An ice-free aquaplanet runs warmer than Earth](#o-23-an-ice-free-aquaplanet-runs-warmer-than-earth) | Climate | Scheduled: M2 steps 11–12, 15 |
| O-24 | [The water cycle is too strong](#o-24-the-water-cycle-is-too-strong) | Moisture | Scheduled: M2 step 15 |
| O-25 | [Clouds are calibrated from layer-mean humidity](#o-25-clouds-are-calibrated-from-layer-mean-humidity) | Clouds | Accepted |
| O-26 | [Strong easterlies aloft at the equator](#o-26-strong-easterlies-aloft-at-the-equator) | Dynamics | Watching |
| O-27 | [Shortwave water absorption depends on the layering](#o-27-shortwave-water-absorption-depends-on-the-layering) | Radiation | Resolved: M2.10.1 |
| O-28 | [More than two layers make the vertical modes non-normal](#o-28-more-than-two-layers-make-the-vertical-modes-non-normal) | Numerics | Accepted |
| O-29 | [The layered core is ~12% slower than the two-layer one](#o-29-the-layered-core-is-12-slower-than-the-two-layer-one) | Performance | Scheduled: M2 step 16 |
| O-30 | [Terrain is smoothed, and spreads over coastal ocean](#o-30-terrain-is-smoothed-and-spreads-over-coastal-ocean) | Terrain | Watching |

---

### O-1: Operator errors peak at the 12 pentagons

**Seen (M1.3, M1.5, M2.2):** several operators have errors that converge in RMS
(halving when resolution doubles) but whose *maximum* stays fixed, located at or
next to the 12 pentagons:

| Operator | Max error (k = 16 → 32 [→ 64]) | RMS error (k = 16 → 32 [→ 64]) |
|---|---|---|
| TRiSK tangential reconstruction | 0.62% → 0.71% | 0.094% → 0.039% |
| Perot cell-vector reconstruction | 2.0% → 2.0% | 0.30% → 0.15% |
| Vertex vorticity (Rossby–Haurwitz) | 6.3% → 5.8% of peak | 1.4% → 0.7% of peak |
| Divergence (M2.2) | 1.24% → 1.15% → 1.13% | 0.13% → 0.058% → 0.028% |
| Laplacian, div(grad) (M2.2) | 1.11% → 1.15% → 1.16% | 0.22% → 0.075% → 0.032% |
| Kinetic energy (M2.2) | 1.74% → 1.53% → 1.50% | 0.25% → 0.081% → 0.032% |

This matches the literature (Peixoto 2016, *Accuracy analysis of mimetic finite
volume operators on geodesic grids*). The divergence and Laplacian maxima sit in the
hexagons right next to the pentagons. Three or more cells away the error is
about 5× smaller (`cargo run --release -p fms-geo --example divergence_error`).

The first flow to show it is Williamson case 2 (M2.3): the l2 height error
converges at second order, but l∞ stalls at ~5e-4 because of sharp error
dipoles at exactly the 12 pentagons (`fms williamson --case 2`, then look at
`output/williamson2/sw.height_error.png`). They are ±1.5 m against a 3 km
fluid depth after 5 days.

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

### O-12: Shallow-water steps are overhead-bound

**Seen (M2.3):** the RK4 shallow-water model runs about 226 steps/s on the 20k
climate grid with 15 threads, about 4.4 ms per step. That's 4 stages of about
15 operator calls each, so each call costs ~70 µs, mostly in thread
scheduling: a single thread is 7× faster on a 2.5k mesh. Setting a minimum
chunk of 4096 elements per thread (`fms_geo::operators::MIN_CHUNK`) won 11×
on small meshes, but large ones are still dominated by per-call overhead.

**Why it's acceptable for now:** correctness first. The roadmap's target is
speed 2 (~340 steps/s) for the full two-layer model, so this needs work. The
semi-implicit step (M2 step 4) makes it worse per step; see O-13.

**Plan:** M2 step 16. Options: fuse operators into fewer passes (one loop over
edges per stage instead of several), run a whole stage per parallel task over
mesh partitions, and precompute per-edge factors like `lₑ/(dₑ r)`.

### O-13: Semi-implicit steps are long but expensive

**Seen (M2.4):** semi-implicit leapfrog runs the Williamson cases stably at
4.7–7.3× the explicit time step, with the same accuracy. But on the 20k climate
grid it is only 1.35–1.9× faster in wall time (case 5: 15.3 s against 21.3 s
for RK4). Each step needs a Helmholtz solve, and conjugate gradient with a
diagonal preconditioner takes ~36 iterations at 7× the explicit step. The
iteration count grows with the gravity-wave Courant number, roughly in
proportion to Δt. Every iteration is several small parallel calls, which pay
the overhead of O-12.

Also: the step must be chosen from the fastest wind *expected* over the run,
not the current wind (`SemiImplicit::stable_dt` takes a design maximum).
Williamson case 5 starts at 20 m/s and accelerates to 40 m/s; a step sized for
the initial wind blew up.

**Why it's acceptable for now:** the scheme is correct, conserves mass to
round-off, and is already a net win. Per-step cost is a performance problem.

**Also seen (M2.7):** Held–Suarez at 10k cells runs at 42 steps/s. Of the
24 ms per step, 18 ms is the adiabatic core (15 CG iterations per solve),
3 ms hyperdiffusion, 1 ms the forcing and 0.7 ms the energy budget. The
1200-day benchmark takes 28 minutes.

**Plan:** M2 step 16, together with O-12. Options: a better preconditioner
(geometric multigrid on the icosahedral refinement hierarchy, or a Chebyshev
polynomial preconditioner, which needs no dot products); fusing each CG
iteration into fewer passes; a looser tolerance (1e-6 halves the iterations
again without measurable error, since mass is restored exactly). M2 step 13
(deriving `dt` from the CFL limit) should use a design maximum wind for the
same reason as above.

### O-14: The undamped core collapses to the grid scale

**Seen (M2.5):** in `fms baroclinic` the storms mature around day 7. Then
fronts and vortices sharpen to single cells (lower-layer vorticity reaches
7×10⁻⁴ s⁻¹, about 7f), winds pass 130 m/s, and the run blows up on day 8.
RK4 at the explicit step does the same through day 7 (132 vs 134 m/s), so this
is the inviscid equations cascading enstrophy to the grid scale with nothing
to remove it, not a time-stepping fault. The semi-implicit run dies first only
because its step was sized for 80 m/s winds.

**Why it's acceptable for now:** every inviscid grid-point model does this;
dissipation is a planned step.

**Plan:** M2 step 6 adds hyperdiffusion. Held–Suarez runs ~1000 days and
cannot work without it (Held & Suarez 1994 used ∇⁸ diffusion), so the
dissipation step was moved ahead of it: Held–Suarez is now step 7.

**Resolved (M2.6):** the diagnosis above was half right. The day-8 blow-up
was winds outrunning a step sized for 80 m/s. With a step sized for faster
winds the inviscid core survives a month, but it does collapse. On the
default 30 m/s-shear jet at 10k cells, vorticity sharpens to
1.4×10⁻³ s⁻¹ (14f) by day 30, winds reach 129 m/s, and the energy error grows
to ~10⁻⁵. With ∇⁴ hyperdiffusion (`crates/atmos/src/diffusion.rs`), the same
storms keep |ζ| below 3×10⁻⁴ s⁻¹ and peak at 105 m/s, and the energy budget
closes to within 1% of what diffusion removes. Blow-ups now stop the run at
the first unhealthy step and write a diagnostic dump (`fms baroclinic`;
`crates/atmos/src/safety.rs`). The default storm run lasts 30 days.

### O-15: Baroclinic growth is ~75% of local theory

**Seen (M2.5):** on a balanced sheared jet, wavenumber-6 perturbations grow at
0.83 /day (k = 16) and 0.78 /day (k = 32), against 1.10 /day from two-level
quasi-geostrophic theory on the local β-plane at 45° N (Phillips 1954; Holton
§8.2), using the model's own internal gravity-wave speed for the deformation
radius. RK4 and semi-implicit agree to three digits. Below the critical shear
(11 m/s) nothing grows, as theory predicts, and the fastest-growing zonal
wavenumber is 6 in both model and theory.

**Why it's acceptable:** local theory assumes waves with no north–south
structure on a flat β-plane. On the sphere the unstable zone is about 30° of
latitude wide; adding that meridional wavenumber to the theory gives a factor
of about 0.73. So local theory is an upper bound, and the tests check 0.6–1.0
of it.

**Revisit if:** a sharper check is wanted. That needs the linear stability
problem of the spherical jet solved numerically (an eigenvalue problem in
latitude) as the reference.

### O-16: No sponge layer

**Seen (M2.6):** the roadmap step listed a sponge layer: extra damping near the
model top that absorbs upward-propagating gravity waves before they reflect
off the lid.

**Why it's acceptable:** with two layers and the top at p = 0 there is no
vertical room for one. A sponge works over several upper levels; here it
would be Rayleigh friction on half the atmosphere. The two-layer model also
has only one internal gravity-wave mode, which the semi-implicit scheme
already handles, and Held–Suarez's Rayleigh friction (M2 step 7) acts at the
surface, not the top.

**Revisit if:** the model gets more vertical levels.

### O-17: Zig-zag vertex modes limit explicit hyperdiffusion

**Seen (M2.6):** a hexagonal C-grid has about twice as many vertices as cells.
The extra rotational modes zig-zag between neighbouring vertices, and the
vertex Laplacian reaches 3.7× higher than the cell Laplacian
(λ·d² = 28.8 vs 7.8 at 10k cells). Hyperdiffusion set to damp the shortest
*cell* wave in `τ` damps these modes `3.7ⁿ` times faster: 14× for ∇⁴, 190×
for ∇⁸. Explicit diffusion, lagged over a leapfrog step, is only stable for
`Δt` below that fastest e-folding time. So the default `τ` is the shortest the
step allows with a 2× margin: 5 h at 20k cells, 7.5 h at 10k. ∇⁸, as in
Held & Suarez (1994), would need `τ ≈ 370 Δt` (about 3 days at 20k cells),
which leaves the 2–4 cell scales almost undamped. Normalizing the rotational
part to the vertex grid instead was tried and fails: then vorticity at the
cell scale is barely damped, and the storms still collapse.

Separately, the TRiSK vector Laplacian is inconsistent (it doesn't converge
pointwise): ~1% RMS error in its divergent part and ~12% in its rotational
part, at every resolution.

**Why it's acceptable:** ∇⁴ at the default strength does the job. Storms live
30 days with smooth fields, and large-scale waves are untouched (wave-6
baroclinic growth changes by 0.5%). The Laplacian errors only matter at
scales where the damping is negligible anyway.

**Revisit if:** the eddies look too damped at mid scales, or Held–Suarez
(M2 step 7) needs a sharper filter. Options: implicit diffusion (one CG
solve per step; its matrix is well conditioned at these strengths), or
filtering the vertex modes separately.

### O-18: Two layers have no surface wind

**Seen (M2.7):** in the Held–Suarez benchmark (1200 days, 10k cells) the
upper layer's jets are close to HS94: 26 m/s at 49° in both hemispheres. The
lower layer is the mean over σ 0.5–1, and its westerlies peak at only
5.7 m/s, against a 20 m/s difference between the layers. Extrapolated
linearly to the surface (`1.5 u_lower − 0.5 u_upper`), that gives about
−4 m/s, easterly under the jets, where HS94's surface winds are westerly.
The zonal-mean surface pressure agrees: it rises steadily from 987 hPa at the
equator to 1019 hPa at the poles. HS94 has subtropical highs and subpolar
lows. The upper-level climate, the temperature contrast (43 K in the lower
layer) and the energy budget all look right.

**Why it's acceptable for now:** it's a limit of two layers. The model has
one wind for the whole lower half of the atmosphere, and friction acts on
that layer mean. Held–Suarez only tests the dynamics, and those behave.

**Plan:** M2 step 10 (land surface and boundary-layer friction) must
diagnose a surface wind rather than use the lower-layer wind: bulk fluxes,
orographic lift and step 12's wind-stress-driven currents all depend on it.
Options: a boundary-layer profile below the lower layer, or a thin third
(boundary) layer. Check the result against the sign of Held–Suarez's
surface winds.

**Update (M2.9.5):** the core now takes any layers, and step 10 will use a
thin third layer. Held–Suarez (k = 16, days 100–200) with five equal-mass
layers has westerlies of +7.3 m/s in the lowest layer under jets of 26 m/s
at 42°: the right sign. With three layers split at σ 0.4 and 0.85, the
jets sit at 32° (in the 0–0.4 layer) with −2 m/s under them, so a thin
bottom layer alone isn't the whole answer; step 10's surface drag and
boundary-layer mixing replace Held–Suarez's friction there.

**Resolved (M2.10.1):** the climate's four layers diagnose a 10 m wind from
the thin bottom layer (log profile, ocean roughness, stability). The bottom
layer has easterly trades (about −6 m/s at 15°) and westerlies (about +6 at
35°) at fixed APE sea temperatures, and on these layers Held–Suarez has
westerlies of +3 to +4 m/s under its jets, where three layers had −2.

### O-19: Radiation is calibrated, not resolved

**Seen (M2.8):** a few of the radiation scheme's parameters are effective
values, set so a single column reproduces Earth's global-mean column, not
values from spectroscopy:

- **CO₂ band growth** `l = 30` cm⁻¹ per e-fold, against the spectroscopic
  ~10 cm⁻¹. Two layers can't resolve the stratosphere cooling and
  readjusting when CO₂ rises, which adds about a third to the instantaneous
  forcing. With l = 30 the model's forcing per doubling is 3.6 W/m², close to
  the adjusted value (Myhre 3.7). Its instantaneous forcing at l = 20 was
  2.5 W/m², close to Earth's real instantaneous value.
- **Convective lapse rate** 4.4 K/km between the layer means, not 6.5: the
  upper layer (0–500 hPa) includes the stratosphere, so Earth's observed
  layer means (272 and 226 K, 10.5 km apart) sit at 4.4. It is fixed, so
  the single column has no lapse-rate feedback (the moist adiabat
  flattening as the surface warms, about −0.5 W/m²/K on Earth). Since M2
  step 9 the 3-D model convects towards the moist adiabat instead, and has
  it; 4.4 K/km still extrapolates the lower layer to the surface air.
- **In-layer lapse rate** 6.5 K/km for the radiation's temperature and
  water-vapour profiles. Through the lower layer's mean this puts its
  radiating bottom ~5 K above the surface air (extrapolated at 4.4 K/km).
  Making the two consistent (M2 step 9 check) cuts downward longwave at the
  surface by 26 W/m², well below Earth's; restoring it needs more
  continuum, which makes the surface's longwave even more sensitive to
  warming (O-24). So 6.5 K/km is effectively a calibrated value too.
- **Water vapour's spectrum** is two regions (band cores, and weak edges
  where CO₂'s wings lie), each a four-term k-distribution, plus a window
  continuum. Calibrating this mattered more than anything else: one grey
  water band gave a runaway (outgoing longwave collapses once it turns
  opaque); a too-strong continuum cut the clear-sky feedback from ~2 to
  ~0.5 W/m²/K.
- **Cloud reflection** in the single column is a fixed deck at the layer
  interface (albedo 0.27); the 3-D model diagnoses its clouds (O-25).
- **Rayleigh scattering** (added at M2 step 9) is one non-absorbing sheet
  at the atmosphere's mass centre reflecting 6% of sunlight: an aquaplanet's
  clear sky then has an albedo of ~0.08, as Earth's oceans do. It doesn't
  yet scale with the column's mass.

Also simplified: vapour pressure capped at a quarter of the air pressure
(very hot planets).

**Why it's acceptable:** with them, the single column reproduces Earth's
global mean: surface 287 K, layers 232/274 K, outgoing longwave 247 W/m²,
21 kg/m² of water, sensible and latent heat 22 and 77 W/m². It also gives
~2 W/m²/K clear-sky feedback, 3.6 W/m² and 2.4 K per CO₂ doubling, and
monotonic, finite climates across the planet parameters' full ranges
(`crates/atmos/tests/radiation.rs`).

**Revisit if:** fantasy settings far from Earth (dense CO₂ atmospheres, hot
or very thin ones) need to be quantitatively right. Then a correlated-k
scheme fitted to line-by-line spectra would replace the effective
parameters, and Rayleigh scattering would scale with the air's mass (thick
atmospheres are now too dark, and so too warm).

### O-20: No cloud or ice feedbacks yet

**Seen (M2.8):** clouds are a fixed deck at the layer interface reflecting
33% of the sunlight that reaches it (30% in the single column, which needs a
lower value to reach the same global mean as the 3-D model). There is no sea ice or snow. So the
model has the Planck, water-vapour and shortwave-absorption feedbacks only,
and warms 2.4 K per CO₂ doubling in a single column. Earth's ~3 K includes
ice-albedo and cloud feedbacks.

**Why it's acceptable for now:** those need moisture (clouds) and surfaces
(ice and snow), which are later steps; 2.4 K is in the plausible range.

**Plan:** M2 step 9 diagnoses clouds from humidity and condensation and
replaces the fixed deck; steps 10–11 add snow and sea ice with their albedo.

**Update (M2.9):** done for clouds: the 3-D model's clouds now come from
humidity and convection and act on both sunlight and longwave (O-25). Snow
falls (and is reported) but doesn't yet lie anywhere; ice and snow cover
remain for steps 10–11.

### O-21: Dissipated kinetic energy leaks out of the planet

**Seen (M2.8):** at M2 step 6 we chose to keep the kinetic energy that
hyperdiffusion removes in the budget only, not return it as heat (and the
Held–Suarez friction, now a placeholder for surface drag, does the same).
With radiation running this is a permanent leak: in equilibrium, absorbed
sunlight exceeds outgoing longwave by the energy friction and diffusion
destroy: 1.3 W/m² (friction 1.1, diffusion 0.2) in the 3-D Earth run. Earth's
frictional dissipation is ~2 W/m², all of it returned as heat. The budget still closes exactly; it just has this sink.

**Why it's acceptable for now:** the leak is small next to the 240 W/m²
flows and steady, so the climate is calibrated with it in place. It is,
though, a third of the forcing from doubling CO₂, so fix it before
interpreting small climate changes.

**Plan:** M2 step 10 replaces the Rayleigh friction with surface drag; return
the dissipated kinetic energy as heat there (frictional heating in the lower
layer, and the diffusion's loss spread where it occurs), so the top of the
atmosphere balances in equilibrium.

**Resolved (M2.10.1):** surface drag, boundary-layer mixing and
hyperdiffusion return their kinetic energy as heat. After five slab years the
top of the atmosphere balances: ASR − OLR is −0.03 W/m², all of it the slab's
uptake, with the energy budget closing to the time-stepping error.

### O-22: The subtropics rain out their own evaporation

**Seen (M2.9):** in the aquaplanet the Hadley circulation is about Earth's
strength (lower-layer flow toward the equator of ~0.35 m/s at 15°, an
overturning of ~70×10⁹ kg/s per hemisphere), and the equator rains
~6.5 mm/day, near Earth's zonal mean. But the subtropics rain 4.4 mm/day
(Earth: ~2) and so do the midlatitudes: rain is nearly flat from 45° S to
45° N, and the subtropics rain out ~75% of their own evaporation instead of
exporting it to the rain belt. Neither the convection scheme nor its
humidity settings change this (four variants tried: faster convection, a
diluted (entraining) parcel, higher condensation thresholds, a moister
convective reference).

**Why:** on Earth the subtropics' water stays in a ~1.5 km boundary layer
under the trade-wind inversion, where it makes shallow, barely raining
cumulus, and the trades carry it to the equator. The model's lower layer is
the whole 500–1000 hPa slab, so evaporation is mixed through half the
atmosphere at once and condenses where it rose; the same missing boundary
layer gives the weak surface winds of O-18.

**Why it's acceptable for now:** the water budget and the large-scale
circulation are right; this is where rain falls, and land (M2 step 10)
reshapes it anyway. But without a fix the world has no subtropical deserts.

**Plan:** M2 step 10, together with O-18: a boundary layer (a profile
under the lower layer, or a thin third layer) that holds the evaporated
water near the surface. Check that the subtropics then rain ≲2 mm/day.

**Update (M2.10.1):** the boundary layer is in, and at fixed APE CONTROL sea
temperatures the subtropics rain 2.5–2.9 mm/day at 25° (down from 4.4), still
above the target. The cause is not the boundary layer. Raising the
condensation thresholds to saturation only reaches 2.5, and venting the
boundary layer harder rains the same water out a layer higher: evaporation
minus rain at 25° stays at 1.2–1.4 mm/day in every experiment. It is set by
the water the circulation exports. The Hadley cell ends near 25°, so it barely
dries the 20–30° band, and transient eddies carry water poleward across 15°
at about three-quarters of the rate the cell brings it equatorward, which also
weakens the ITCZ (6 mm/day). This doesn't change with resolution (k = 16, 24
and 32 agree: 5.8, 5.1 and 5.3×10⁸ kg/s at 15°), and stationary eddies and
time variations of the zonal mean are near zero. The APE paper's text gives
no eddy flux at 15° to compare with.

**Revisit:** M2 step 15's tuning pass, with the width of the Hadley cell and
the tropical eddy moisture flux as the first things to check.

### O-23: An ice-free aquaplanet runs warmer than Earth

**Seen (M2.9):** with Earth's orbit and CO₂ and Earth-like clouds and
humidity, the 3-D aquaplanet settles at ~297 K after 5 years (Earth: 288;
still warming ~0.1 K/year), with the tropics at ~301 K (Earth's oceans
~300) and the poles at ~279 K (Earth's ~250 under ice). Its planetary albedo
is Earth's (0.29), but its clear
sky is a dark ocean everywhere (albedo ~0.08, Earth's ~0.15 with ice, snow
and deserts), and without ice there is no ice–albedo feedback to cool the
high latitudes.

**Why it's acceptable for now:** it is what an ice-free Earth-like
aquaplanet should do, and the cloud parameters were calibrated for
realistic clouds, humidity and cloud radiative effects rather than to force
288 K. Forcing it now would leave the clouds wrong once land and ice
arrive.

Without ocean heat transport the slab also can't cool the tropics or warm
the poles.

**Plan:** M2 steps 10–12 add land, snow, sea ice and ocean heat transport
(which cools an aquaplanet's equator by ~8 K and moves its ice edge
10–15° poleward; Bhatnagar et al. 2026); M2 step 15's tuning pass then
checks the global mean against Earth's 288 K.

**Update (M2.10.1):** on four layers with recalibrated clouds the ice-free
slab aquaplanet settles at 290.7 K, about 3 K warmer than Earth (albedo 0.322,
OLR 230.7 W/m²). Ice and land will cool it further, so step 15 sets the final
temperature.

### O-24: The water cycle is too strong

**Seen (M2.9):** the aquaplanet evaporates and rains ~4.7 mm/day (latent
heat ~138 W/m²; Earth ~2.7–3 mm/day, ~85 W/m²) and its sensible heat flux
is small (~9 W/m²; Earth ~20). The surface absorbs Earth's ~160 W/m² of
sunlight, but loses only ~13 W/m² of net longwave (Earth: ~55; ~30 under
clear skies, so the clear sky is most of it), and the rest must leave by
evaporation. Partly that is the warm, moist aquaplanet
(O-23); partly the radiation: in the single column, downward longwave at
the surface rises ~8 W/m² per kelvin of warming at fixed relative humidity,
~30% more than Earth-like models (~6), mostly through the window
continuum, and evaporation rises ~5.5%/K (models: 2–3%/K). The single
column's clear-sky surface budget itself matches Earth's (downward 299 W/m²
and net loss 82 W/m² at 286 K; Earth ~314 and ~84 at 288 K).

**Why it's acceptable for now:** the budgets close, and the excess is
largest in the warm aquaplanet that land and ice will cool.

**Plan:** M2 step 15's tuning pass. Options: a boundary-layer sublayer in the
radiation (which would let the continuum act on the near-surface air
alone), or recalibrating the continuum against both the surface and
top-of-atmosphere budgets.

### O-25: Clouds are calibrated from layer-mean humidity

**Seen (M2.9):** the 3-D model's two cloud sheets (high at σ = 0.25, low at
0.8; now 0.3 and 0.9) take their cover from each layer's mean relative humidity (linear
between two bounds) plus anvils from convective rain. The bounds,
albedos, emissivities and the large-scale condensation thresholds are
calibrated so that layer humidities match the radiation's calibration
(0.31 upper, 0.63 lower) and cover (~0.64), planetary albedo (~0.30) and
cloud radiative effects (−70 shortwave, +15 longwave W/m²) are Earth-like
for an aquaplanet. Cover was very sensitive to the formula: SPEEDY's
quadratic, with the overcast bound at 0.8, gave a cold (277 K), 89%-cloudy
planet with almost no circulation, through a feedback loop (cloud cools the
surface, convection stops, the lower layer fills to its condensation
threshold, more cloud). The linear ramp is gentler near saturation.
The longwave effect is about half of Earth's (+28).

**Recalibrated (M2.10.1), on four layers:** the high sheet now reads the
moister σ 0.25–0.5 layer, so its humidity bounds rose to 0.45–0.85 and the low
sheet's lower bound to 0.5. At fixed APE CONTROL sea temperatures (240 days):
imbalance +14 W/m² (APE ~13), albedo 0.307, cloud effects −70 and +19 W/m²,
cover 0.62. Clear-sky absorbed minus outgoing radiation is +65 W/m², so the
net cloud effect has to be about −52 for that imbalance, which is what holds
the longwave effect below Earth's.

**Why it's acceptable:** a two-layer model can't resolve cloud layers; its
cloud scheme is necessarily a calibrated diagnosis, and its cloud feedback
under climate change follows the layers' humidity, so it is not a
prediction.

**Revisit if:** the climate's response to CO₂ or orbit matters
quantitatively (cloud feedback), or M2 step 15's benchmark shows the
cloud pattern off. Low stratocumulus under subsidence inversions (absent:
they need a boundary layer, O-22) would be the first addition.

### O-26: Strong easterlies aloft at the equator

**Seen (M2.9):** the aquaplanet's upper layer blows from the east at
~26 m/s over the equator (Earth: ~5). The jets (24 m/s at 45°) and the
lower layer's trades and westerlies are reasonable.

**Why it's acceptable for now:** it doesn't visibly affect the rain belt
or temperatures. A likely cause is that convection mixes heat and water
between the layers but not momentum (convective momentum transport, which
brakes upper-level winds in the tropics on Earth).

**Revisit if:** the equatorial winds start to matter (tropical waves,
monsoons, M2 step 15), by mixing momentum between the layers along with
the convective mass exchange.

### O-27: Shortwave water absorption depends on the layering

**Seen (M2.9.5):** the same column cut into four or six layers instead of
two (each at the two-layer profile's temperature, fixed relative
humidity) puts 7.7 and 12.3 W/m² less sunlight on the surface, though its
water agrees to 0.5% and its longwave to a few W/m². Water vapour's
near-infrared absorption saturates as `√w` layer by layer, so splitting a
layer raises it (`Σ√wₖ > √Σwₖ`); and each layer takes its share of all the
light rather than of what's left in the near-infrared, which over-counts
it further.

**Why it's acceptable for now:** two layers are calibrated with it in
place, and the climate stays on two layers until step 10.

**Plan:** M2 step 10, before recalibrating on three layers: carry the
near-infrared band separately through the adding method, with each slab's
optical depth from the cumulative water path (`τ ∝ √W_below − √W_above`),
so absorption doesn't depend on the layering. Check that more layers
converge.

**Fixed (M2.10.1):** sunlight now runs through the adding method in two
bands, the visible (ozone) and the near-infrared (water, from the
cumulative slant path). Four and six layers now put 0.03 and 0.01 W/m² less
sunlight on the surface than two. Recalibrated to the same column
absorption: ozone 0.03 → 0.046 (it now acts on the visible band only),
water 0.04 → 0.062; the reference column absorbs 49.3 W/m² against 49.4
before, and its response to halving the water is the same within 0.5 W/m².

### O-28: More than two layers make the vertical modes non-normal

**Seen (M2.9.5):** with two layers, the gravity-wave matrix `A` is
self-adjoint in the mass-weighted inner product (`δσₖ A` symmetric to
1e-11). With more, and stratification, it isn't (0.4% for three layers, 2%
for five), although the linearization matches the model's own tendencies
to 1e-12. The cause is the interface θ: Arakawa & Suarez (1983) conserve
energy with any interface value, but the arithmetic mean of the two layers
doesn't give a symmetric linear operator once there are several
interfaces. Using the symmetric part instead made the semi-implicit
solver inconsistent (spurious wind over the mountain test rose from 0.105
to 0.166 m/s), so the solver decomposes the general matrix.

**Why it's acceptable:** the eigenvalues stay real and distinct (or
exactly repeated, for an unstratified atmosphere), energy is still
conserved exactly by the spatial scheme, and all the layered tests pass.

**Revisit if:** the vertical modes misbehave with many layers (computational
modes, noise in θ between layers): Arakawa & Suarez's weighted interface θ
is the alternative.

### O-29: The layered core is ~12% slower than the two-layer one

**Seen (M2.9.5):** on two layers, with arithmetic identical to step 9's
(bit-for-bit results), the layered core takes 7.31 ms per moist step
against 6.60, and 3.09 per dry Held–Suarez step against 2.75, at 2.5k cells
([BENCHMARKS.md](BENCHMARKS.md)). Likely causes: small vectors allocated
per cell and step in the moist physics, radiation and interface humidity
ratios (`Vec` per column where step 9 had fixed arrays), per-layer vectors of
vectors, and copying the per-cell Exner and geopotential scratch into the
layers.

**Why it's acceptable for now:** it's a constant factor, well within the
performance budget's margin at this stage.

**Plan:** M2 step 16 (performance pass): profile first; then stack arrays
of `MAX_LAYERS` in the per-column physics, and reuse buffers.

### O-30: Terrain is smoothed, and spreads over coastal ocean

**Seen (M2.10.2.1):** a stratified atmosphere at rest over Earth's unsmoothed
terrain at 10k cells holds a persistent wind error of 5.7 m/s at the steepest
cells, measured by damping the winds with a one-day timescale for 15 days. It
is 15.6 m/s if released and left undamped. So the atmosphere sees the surface
height after ten conservative diffusion passes (`smooth_height`), which give:

| Passes | Highest cell | Persistent error |
|---|---|---|
| 0 | 5300 m | 5.70 m/s |
| 5 | 4982 m | 1.27 m/s |
| 10 | 4842 m | 0.87 m/s |
| 20 | 4468 m | 0.66 m/s |
| 40 | 3786 m | 0.58 m/s |

Two side effects:

- Terrain spreads onto all-ocean cells next to steep coasts: of 6802 at 10k
  cells, 236 end up above 100 m and 8 above 500 m (Greenland and Antarctic
  margins).
- The resting start isn't discretely balanced, because the layers sample a
  continuous profile. It radiates gravity waves of ~5 m/s in the first day,
  independent of how steep the terrain is.

**Why it's acceptable for now:** the remaining error is well below the winds
that matter, and it does not grow with resolution (0.83 m/s at 20k cells with
20 passes). The start's waves disperse and are damped in spin-up. Mountains
lose some height at the grid scale, as in any GCM; sub-grid effects are step
10.2.6's job.

**Checked with ocean tiles (M2.10.2.2):** holding all-ocean cells at sea
level during smoothing brings steep steps back at the coast, raising the
persistent wind error from 0.87 to 1.18 m/s at 10k cells, and loses the
filter's conservation. Leaving them lifted costs little. After 400 days on
Earth's terrain at 2.5k cells, the lifted cells' surface temperature against
sea-level ocean within 3° of latitude is:

| Lifted height | Cells | Mean height | Temperature difference |
|---|---|---|---|
| 50–100 m | 93 | 71 m | −0.56 K |
| 100–200 m | 53 | 134 m | −0.64 K |
| 200–500 m | 30 | 293 m | −1.06 K |
| above 500 m | 2 | 722 m | −2.08 K |

These cells all sit beside land, so continental air explains part of the
difference. We keep the lift.

**Revisit when:** the slab ocean (M2 step 11), which gives cells at height
sea ice and heat transport: check the same cells again.
