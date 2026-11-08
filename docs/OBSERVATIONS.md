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
| O-18 | [Two layers have no surface wind](#o-18-two-layers-have-no-surface-wind) | Dynamics | Scheduled: M2 step 10 |

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

