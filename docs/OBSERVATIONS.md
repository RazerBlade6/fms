# Observations

Known limitations, measured quirks and ideas we've chosen not to act on yet.
Scheduled work is in [ROADMAP.md](ROADMAP.md).

Each entry says what we saw, why it's acceptable and when to revisit, in a few
lines. Numbers never change, because code comments cite them. A resolved entry
shrinks to one line.

## Grid and numerics

### O-1: Operator errors peak at the 12 pentagons
*Scheduled: M2.14.* TRiSK reconstruction, Perot, vorticity, divergence, Laplacian
and kinetic energy converge in RMS, but their maxima (1–6%) stay fixed at or next
to the pentagons, as in Peixoto 2016. Williamson 2 shows ±1.5 m dipoles on a 3 km
depth after 5 days. The risk is grid imprinting in long runs. If visible:
Heikes–Randall optimization or Peixoto's corrected weights.

### O-2: Cell sizes vary by up to 1.76×
*Watching.* Largest over smallest is 1.46 at 20k cells and 1.76 at 1.6M, costing
~17% of the maximum `dt`. Revisit if `dt` is the bottleneck: density-weighted Lloyd
or spring dynamics.

### O-3: Grid frequencies need small prime factors
*Accepted.* Multilevel Lloyd refines by primes ≤ 7; others fall back to plain
Lloyd, O(f²) iterations. We choose the frequencies, and the UI should only offer
good ones.

### O-4: Terrain cells are assigned whole to climate cells
*Accepted.* A climate cell's child area is within ±25% of its own. Global sums
stay exact and smooth means are accurate. Exact overlap only matters if per-cell
conservation does.

### O-17: Zig-zag vertex modes limit explicit hyperdiffusion
*Accepted.* The vertex Laplacian reaches 3.7× the cell one, so those modes set the
stable strength. ∇⁸ would leave the 2–4 cell scales almost undamped, so we use ∇⁴.
The TRiSK vector Laplacian also doesn't converge pointwise (~12% RMS rotational).
Revisit if eddies look too damped: implicit diffusion, or filter the vertex modes.

### O-28: More than two layers make the vertical modes non-normal
*Accepted.* The arithmetic interface θ makes the gravity-wave matrix non-symmetric
(0.4% for three layers, 2% for five), so the solver decomposes the general matrix.
Eigenvalues stay real and energy is conserved. Revisit if modes misbehave with many
layers: Arakawa & Suarez's weighted interface θ.

### O-40: The pinned climate hash doesn't see convection's reference
*Accepted.* `a_world_without_land_is_the_aquaplanet_bit_for_bit` (5 days, 642 cells)
is unchanged by adding 1 mK to every `lift` result. It guards dynamics, radiation
and the in-layer profile. Changes to `lift` need their own check.

### O-41: `KineticHeating::rate` uses a parallel float sum
*Scheduled: M4.1.* Above 4,096 edges its result can depend on how rayon splits. It
only feeds the `dissipation` diagnostic.

## Planet and tooling

### O-5: Kepler's third law is not enforced
*Accepted.* Year length and semi-major axis are independent, for fantasy worlds. If
realism is wanted: an optional star mass.

### O-6: The equation of time is ignored
*Accepted.* Irrelevant for climate (about ±16 min on Earth).

### O-7: Daily-mean insolation uses one orbital state per day
*Accepted.* Negligible for any year of at least 20 days.

### O-9: Vector fields render as arrows only
*Watching.* Streamlines read better for large-scale circulation. Revisit if
general-circulation images are hard to read.

### O-10: Field names must be `'static` strings
*Scheduled: M4.2.* Switch `FieldSpec::name` to `Cow<'static, str>` so saves can keep
unknown fields.

## Performance

### O-8: Scheduler overhead is ~119 µs per step
*Scheduled: M2.16.* 20 trivial systems in 2 stages, ~4% at speed 2. Options: cache
stage plans, or run tiny stages serially.

### O-12: Shallow-water steps are overhead-bound
*Scheduled: M2.16.* About 70 µs per operator call, mostly thread scheduling. Options:
fuse operators into fewer passes, run stages over mesh partitions, precompute edge
factors. M2.11.1.5 didn't reach operator fusion.

### O-13: Semi-implicit steps are long but expensive
*Scheduled: M2.16.* 4.7–7.3× the explicit step, but only 1.35–1.9× faster in wall
time, because CG iterations grow with Δt (~36 at 7×). Options: multigrid or Chebyshev
preconditioning, looser tolerance (1e-6 halves iterations), fused iterations. Size the
step from a design maximum wind: flows accelerate, and Williamson 5 blew up otherwise.

### O-29: The layered core is ~12% slower than the two-layer one
*Scheduled: M2.16.* Likely per-cell small vectors and scratch copies. Plan: profile,
then stack arrays of `MAX_LAYERS` and reuse buffers.

## Dynamics

### O-14: The undamped core collapses to the grid scale
*Resolved at M2.6* by ∇⁴ hyperdiffusion and blow-up detection.

### O-15: Baroclinic growth is ~75% of local theory
*Accepted.* 0.78–0.83 /day against 1.10 from local QG theory, which is an upper bound
that ignores meridional structure (a factor of about 0.73). Tests check 0.6–1.0 of it.
A sharper check needs the spherical eigenvalue problem.

### O-16: No sponge layer
*Accepted.* There was no room for one on two layers, and the four-layer climate
hasn't been checked for a need. Revisit with more levels.

### O-18: Two layers have no surface wind
*Resolved at M2.10.1* by the thin bottom layer and a diagnosed 10 m wind.

### O-26: Strong easterlies aloft at the equator
*Watching.* About 26 m/s on the two-layer aquaplanet (Earth ~5), likely from no
convective momentum transport. Doesn't visibly affect rain or temperature. Revisit if
equatorial winds matter (M2.15).

## Radiation and energy

### O-19: Radiation is calibrated, not resolved
*Accepted.* Effective values set against Earth's global-mean column: CO₂ band growth 30
cm⁻¹ (spectroscopic ~10), in-layer lapse rate 6.5 K/km, the single column's 4.4 K/km,
a two-region k-distribution, a Rayleigh sheet that doesn't scale with column mass, and
vapour pressure capped at a quarter of air pressure. Gives 3.6 W/m² and 2.4 K per CO₂
doubling. Revisit if worlds far from Earth must be quantitatively right: correlated-k
from line-by-line spectra, and Rayleigh scaled with mass.

### O-20: No cloud or ice feedbacks yet
*Resolved:* clouds at M2.9, land snow at M2.10.2.5, sea-ice albedo at M2.11.5.

### O-21: Dissipated kinetic energy leaks out of the planet
*Resolved at M2.10.1:* it returns as heat, and the top of the atmosphere balances.

### O-27: Shortwave water absorption depends on the layering
*Resolved at M2.10.1* by the two-band adding method on the cumulative path.

## Moisture and clouds

### O-22: The subtropics rain out their own evaporation
*Improved at M2.10.1; rest at M2.15.* At fixed APE CONTROL SSTs, 25° rains 2.5–2.9
mm/day (Earth ~2). The circulation sets it, not the boundary layer: evaporation minus
rain stays 1.2–1.4 mm/day in every experiment, the Hadley cell ends near 25°, and eddies
carry water poleward across 15° at ¾ of the cell's return flow. Independent of
resolution. Start with Hadley width and the tropical eddy moisture flux.

### O-24: The water cycle is too strong
*Scheduled: M2.15.* The aquaplanet evaporates ~4.7 mm/day (Earth ~2.7–3). Downward
longwave at the surface rises ~8 W/m²/K (models ~6), mostly the window continuum.
Options: a radiation sublayer in the boundary layer, or recalibrating the continuum.

### O-25: Clouds are calibrated from layer-mean humidity
*Accepted.* Bounds, albedos and emissivities are tuned to Earth-like cover (~0.62) and
cloud effects (−70, +19 W/m²); the longwave effect is below Earth's +28. Cover is very
sensitive to the formula, so cloud feedback isn't a prediction. Revisit if the response
to CO₂ matters or M2.15 shows the pattern off; low stratocumulus would come first.

### O-35: The column cases find two gaps at the inversion
*Accepted; kept as ignored tests (`cargo gaps`).* BOMEX rains 11.8 mm/day (observed
none) because Betts–Miller has no inhibition; this may be part of O-22. DYCOMS entrains
41 mm/s against 3.8–5.9 from large-eddy runs, and has cover 0.31 instead of overcast,
because extrapolating the layer above leaves 0.36 K of the real 8.5 K jump. BOMEX
cover with deep convection off is 0.80 (large-eddy 0.13). Candidates at M2.15: an
inhibition trigger, and an inversion diagnosed from the column.

## Climate and terrain

### O-11: Land below sea level counts as ocean
*Watching.* Land and sea split at elevation 0; only the Caspian shows at climate
resolution. The two-layer ocean doesn't need a mask. M2.12's gyres do: take the
below-sea-level region connected to the open sea, leaving closed basins to hydrology.
Narrow straits make connectivity unreliable.

### O-23: An ice-free aquaplanet runs warmer than Earth
*Scheduled: M2.12, M2.15.* The ice-free four-layer slab settles at 290.7 K. Ice and
ocean transport cool it further; M2.15 sets the final temperature.

### O-30: Terrain is smoothed, and spreads over coastal ocean
*Watching; recheck at M2.11.8.* Ten passes cut the persistent σ-coordinate wind error on
Earth at 10k cells from 5.7 to 0.87 m/s, and it doesn't grow with resolution. Terrain
spreads onto all-ocean coastal cells (236 of 6802 above 100 m, 8 above 500 m), which run
0.6–2 K colder. Pinning them at sea level is worse (1.18 m/s, and loses conservation).

### O-31: Land has one soil layer
*Watching.* No deep-ground store carries summer heat into winter. Revisit at M2.15 if
continental interiors are too extreme or too early. The 10.2.8 close-out wrote annual
means only, so that comparison hasn't been made.

### O-32: Too little rain over the Amazon
*Scheduled: M2.15.* About 3.4 mm/day at 10k cells (Earth 6 or more), and Borneo 4.4 (8 or
more). Going from 2.5k to 10k cells doesn't fix it. Dry-soil surface air is also still
buoyant as if saturated, in `exchange`.

### O-33: The snow pack is simple
*Watching; M2.15, M4.* No heat capacity or insulation, frost lands as liquid in the
bucket, rain passes through, and nothing limits the pack on ice sheets (1627 kg/m² mean
after six years, 3013 at most). Antarctic land south of 75° is 236.9 K, still ~13 K above
the ~224 K noted from memory. M2.11.7 tested one cause at a time, against twin-run noise
of 0.01 K (annual), ~1 K (seasonal) and 1.6 K over Greenland:
- clear skies over ice sheets, a bound: −11.9 K;
- freeze-dry low cloud, adopted: −7.1 K from 244.0 K, and the 60–80°N winter −4.4 K;
- stability across the whole bottom layer: +0.7 K, the planet +2.4 K. The 1.2 km layer
  warms and radiates down; a shallow inversion needs more layers;
- ice-sheet skin conductance 7 → 1 W/m²/K (insulation): +0.5 K.

Untested: terrain ~800 m too low at 2.5k cells.

### O-34: More rain on the lee of a ridge than the windward side
*Scheduled: M2.15.* In idealized 2.5k-cell runs the lee gets more (3.5 against 2.5–2.9
mm/day): low-level westerlies are only 1–3 m/s, so rain follows storms and stationary
waves. The scheme responds correctly to the wind it's given. At 10k cells on Earth the
contrast has the right sign in the Pacific Northwest and Himalaya, none in Patagonia.

### O-36: Too much sea ice is expected until M2 step 12
*Scheduled: M2.12.* Only diffusion moves heat until then. Tests must not require Earth's
ice extent or SSTs, and excess ice isn't a reason to retune. On aquaplanets the
comparison with Bhatnagar et al.'s no-transport and diffusion-only cases holds in sign
and rough size, with no hemispheric lock-in. Revisit with full transport on.

### O-37: The ocean jump barely speeds up the ice
*Open: choose the jump setting before M2.11.8.* Against year 20 of a true-time tilted
aquaplanet at 2.5k cells, the default (three 3-year jumps) is −0.33 K and leaves ice at
6.0 and 7.7 million km² (target 4.2 and 4.0). Three 6-year jumps at the same cost: −0.05
K, 4.7 and 5.3. No schedule helps without tilt. The default was validated without ice,
so it needs an Earth check before changing.

### O-38: Perennial ice has no equilibrium in 20 years
*Watching; revisit at M2.12.* On a no-tilt aquaplanet the poles never melt: snow piles
up and floods into ice past the 10 m limit, as in the authors' code, and the top of the
atmosphere stays near −1 W/m². Energy still closes, and the ice edge is set by the much
thinner ice. A long-lived perennial world would need a snow limit or a thinner start.

### O-39: Ice sheets are seeded from the start temperature
*Watching; M2.11.7, M2.11.8, M2.15, M4.* A year of full snow cover alone made almost no
sheets, because Antarctica's pack melts out under seasonal albedo. So land starting below
267 K starts as a sheet. Earth at 2.5k cells: Antarctica −11.4 K, Greenland −7.3 K, the
planet −1.12 K, and snow-free land −0.64 K (noise ~0.1 K). The threshold was chosen on
Earth and has little margin: the coldest northern lowland starts at 268.1 K. Two Tibetan
cells are seeded too, and are melting away.
