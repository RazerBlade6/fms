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
*Scheduled: M2.15.* 20 trivial systems in 2 stages, ~4% at speed 2. Options: cache
stage plans, or run tiny stages serially.

### O-12: Shallow-water steps are overhead-bound
*Scheduled: M2.15.* About 70 µs per operator call, mostly thread scheduling. Options:
fuse operators into fewer passes, run stages over mesh partitions, precompute edge
factors. M2.11.1.5 didn't reach operator fusion.

### O-13: Semi-implicit steps are long but expensive
*Scheduled: M2.15.* 4.7–7.3× the explicit step, but only 1.35–1.9× faster in wall
time, because CG iterations grow with Δt (~36 at 7×). Options: multigrid or Chebyshev
preconditioning, looser tolerance (1e-6 halves iterations), fused iterations. Size the
step from a design maximum wind: flows accelerate, and Williamson 5 blew up otherwise.

### O-29: The layered core is ~12% slower than the two-layer one
*Scheduled: M2.15.* Likely per-cell small vectors and scratch copies. Plan: profile,
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
equatorial winds matter (M2.16).

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

### O-48: Budget-test energy residual moves with the radiation interval
*Watching.* The 45-day energy-closure tests on 642 cells leave -0.05 to -0.11 W/m²
depending on how often radiation runs (Ekman test, 21 steps a day: -0.089, -0.111,
-0.071, -0.050 for every 1, 2, 3, 6 steps), with no trend. Cause not found. A year of
Earth at 2,562 cells leaves +0.007. The test limit was raised from 0.1 to 0.15 at M2.13.2.

## Moisture and clouds

### O-22: The subtropics rain out their own evaporation
*Improved at M2.10.1; rest at M2.16.* At fixed APE CONTROL SSTs, 25° rains 2.5–2.9
mm/day (Earth ~2). The circulation sets it, not the boundary layer: evaporation minus
rain stays 1.2–1.4 mm/day in every experiment, the Hadley cell ends near 25°, and eddies
carry water poleward across 15° at ¾ of the cell's return flow. Independent of
resolution. Start with Hadley width and the tropical eddy moisture flux.

### O-24: The water cycle is too strong
*Scheduled: M2.16.* The aquaplanet evaporates ~4.7 mm/day (Earth ~2.7–3). Downward
longwave at the surface rises ~8 W/m²/K (models ~6), mostly the window continuum.
Options: a radiation sublayer in the boundary layer, or recalibrating the continuum.

### O-25: Clouds are calibrated from layer-mean humidity
*Accepted.* Bounds, albedos and emissivities are tuned to Earth-like cover (~0.62) and
cloud effects (−70, +19 W/m²); the longwave effect is below Earth's +28. Cover is very
sensitive to the formula, so cloud feedback isn't a prediction. Revisit if the response
to CO₂ matters or M2.16 shows the pattern off; low stratocumulus would come first.

### O-35: The column cases find two gaps at the inversion
*Accepted; kept as ignored tests (`cargo gaps`).* BOMEX rains 11.8 mm/day (observed
none) because Betts–Miller has no inhibition; this may be part of O-22. DYCOMS entrains
41 mm/s against 3.8–5.9 from large-eddy runs, and has cover 0.31 instead of overcast,
because extrapolating the layer above leaves 0.36 K of the real 8.5 K jump. BOMEX
cover with deep convection off is 0.80 (large-eddy 0.13). Candidates at M2.16: an
inhibition trigger, and an inversion diagnosed from the column.

## Climate and terrain

### O-11: Land below sea level counts as ocean
*Resolved at M2.12.5* by the sea mask: below-sea-level terrain is sea only in the
world ocean or a body of at least 0.005% of the sphere. Left open: inlets cut off
by straits narrower than a terrain cell become land if small (Sea of Marmara, Lake
Maracaibo, the Gulf of Ob at coarse terrain), which moves at most a quarter of one
climate cell; lake beds below sea level (Superior, Baikal) are land at height 0
until lakes arrive (M3.6); and Lake Superior's bed is within 5% of the threshold.

### O-23: An ice-free aquaplanet runs warmer than Earth
*Scheduled: M2.16.* The ice-free four-layer slab settles at 290.7 K. Ice and
ocean transport cool it further; M2.16 sets the final temperature. With currents
the tilted aquaplanet is heading for about 294 K.

### O-30: Terrain is smoothed, and spreads over coastal ocean
*Watching.* Ten passes cut the persistent σ-coordinate wind error on Earth at 10k cells
from 5.7 to 0.87 m/s, and it doesn't grow with resolution. Terrain spreads onto all-ocean
coastal cells (at 2.5k cells 92 above 100 m, at 10k 184), which run ~0.9 K colder. Pinning
them at sea level is worse (1.18 m/s, and loses conservation). Fewer passes (M2.11.8, at
2.5k): 3 give 1.39 m/s against 0.70, 1 gives 3.23 m/s, and none blows up over Tibet.

### O-31: Land has one soil layer
*Watching.* No deep-ground store carries summer heat into winter. Revisit at M2.16 if
continental interiors are too extreme or too early. The 10.2.8 close-out wrote annual
means only, so that comparison hasn't been made.

### O-32: Too little rain over the Amazon
*Scheduled: M2.16.* About 3.4 mm/day at 10k cells (Earth 6 or more), and Borneo 4.4 (8 or
more). Going from 2.5k to 10k cells doesn't fix it. Dry-soil surface air is also still
buoyant as if saturated, in `exchange`.

### O-33: The snow pack is simple
*Watching; M2.16, M4.* No heat capacity or insulation, frost lands as liquid in the
bucket, rain passes through, and nothing limits the pack on ice sheets (1627 kg/m² mean
after six years, 3013 at most). Antarctic land south of 75° is 236.9 K, still ~13 K above
the ~224 K noted from memory. M2.11.7 tested one cause at a time, against twin-run noise
of 0.01 K (annual), ~1 K (seasonal) and 1.6 K over Greenland:
- clear skies over ice sheets, a bound: −11.9 K;
- freeze-dry low cloud, adopted: −7.1 K from 244.0 K, and the 60–80°N winter −4.4 K;
- stability across the whole bottom layer: +0.7 K, the planet +2.4 K. The 1.2 km layer
  warms and radiates down; a shallow inversion needs more layers;
- ice-sheet skin conductance 7 → 1 W/m²/K (insulation): +0.5 K.

Terrain isn't the answer (M2.11.8): the region averages 1,860 m at 2.5k cells and 2,030 m
at 10k, both ~237 K, and 1–3 smoothing passes cool it only 1–1.5 K. The ~224 K is the
South Pole (2,835 m), while this mean includes the warmer ice shelves. Against ERA5
(O-42) the region is +1.9 K, and +0.6 K moved to ERA5's height (2,070 m), so most of
the "13 K" was the point comparison.

### O-34: More rain on the lee of a ridge than the windward side
*Scheduled: M2.16.* In idealized 2.5k-cell runs the lee gets more (3.5 against 2.5–2.9
mm/day): low-level westerlies are only 1–3 m/s, so rain follows storms and stationary
waves. The scheme responds correctly to the wind it's given. At 10k cells on Earth the
contrast has the right sign in the Pacific Northwest and Himalaya, none in Patagonia.

### O-36: Too much sea ice is expected until M2 step 12
*Resolved at M2.12.8.* Only diffusion moves heat until then. Tests must not require Earth's
ice extent or SSTs, and excess ice isn't a reason to retune. On aquaplanets the
comparison with Bhatnagar et al.'s no-transport and diffusion-only cases holds in sign
and rough size, with no hemispheric lock-in. On Earth (2.5k cells, year 20) the seasons
are in phase. Against NSIDC's area (O-42) the northern minimum is close (3.9 against 4.3
million km²), but the southern is 4.4 against 2.4, and winter peaks of ~20 million km² are
7 too large in the north and 6 in the south. Revisit with full transport on.
With currents (Earth, 2,562 cells, year 9 of the spin-up preset) northern ice is closer
to NSIDC: annual mean 10.7 against 9.4 M km² (13.4 before), winter peak 15.3 against
13.3 (19.3 before), summer minimum 2.4 against 4.3 (3.7 before). But southern ice has
almost gone; see O-47.

### O-37: The ocean jump barely speeds up the ice
*Resolved at M2.11.8:* default kept. Against year 20 of true time at 2.5k
cells, with freeze-dry: on the tilted aquaplanet the default (three 3-year jumps) is
−0.26 K with southern ice 6.8 million km² (target 4.6), and three 6-year jumps −0.03 K
and 6.2. On Earth the default is +0.44 K, 6-year jumps +0.66 K, and plain true time at
year 6 +0.24 K. Jumps always bring the deep ocean closer. No schedule helps without
tilt. Earth hasn't settled at year 20 either (−0.68 W/m² at the top).

### O-38: Perennial ice has no equilibrium in 20 years
*Watching; revisit at M2.12.* On a no-tilt aquaplanet the poles never melt: snow piles
up and floods into ice past the 10 m limit, as in the authors' code, and the top of the
atmosphere stays near −1 W/m². Energy still closes, and the ice edge is set by the much
thinner ice. A long-lived perennial world would need a snow limit or a thinner start. With currents (M2.12.4) the polar ice thins over six years where it
used to pile up, but the cap hasn't opened.

### O-39: Ice sheets are seeded from the start temperature
*Watching; M2.11.7, M2.11.8, M2.16, M4.* A year of full snow cover alone made almost no
sheets, because Antarctica's pack melts out under seasonal albedo. So land starting below
267 K starts as a sheet. Earth at 2.5k cells: Antarctica −11.4 K, Greenland −7.3 K, the
planet −1.12 K, and snow-free land −0.64 K (noise ~0.1 K). The threshold was chosen on
Earth and has little margin: the coldest northern lowland starts at 268.1 K. Two Tibetan
cells are seeded too, and are melting away. At 265 K a Greenland edge cell and two on
Ellesmere Island are lost; at 269 K the Antarctic Peninsula and Svalbard are gained, but
so is interior Alaska (M2.11.8).

### O-42: Earth against observations
*Watching; M2.12, M2.16.* First comparison with observations (`scripts/compare_earth.py`;
year 20 of true time at 2.5k cells; ERA5 temperature and wind, GPCP rain, NSIDC ice,
CERES radiation and cloud):
- Global −1.7 K and rain +25% (3.37 against 2.69 mm/day). RMS over latitude and month:
  T 3.5 K, rain 1.1 mm/day, absorbed sunlight 26 W/m², outgoing longwave 18 W/m².
- **Clouds are the largest error.** Albedo 0.338 against 0.291, so 225 W/m² absorbed
  against 241. Shortwave cloud effect −64.5 against −45.7 W/m², longwave +16.5 against
  +27.9. Over tropical oceans similar cover (0.72 against 0.66) reflects twice as much
  (−88 against −45): too bright and too low. Polar cloud is too thin (0.3–0.46 against
  0.5–0.77), and the Sahara is cloudy (0.47 against 0.23).
- North of 55°N −6 to −9 K, with seasons 8–20 K too wide; the Arctic Ocean −8 K (−5 K
  against NCEP). Greenland −5 K at 950 m (ERA5 1,800 m). 45–60°S +1 to +3 K.
- The subtropics are too wet: 3.1–3.7 mm/day at 20–31° against 1.8–2.3, the Sahara 3.9
  against 0.2 and Tibet 5.0 against 1.3. The Amazon is too dry (4.4 against 6.3).
- The circulation is weak: the bottom layer's trades are 1.7–2.7 m/s and the southern
  westerlies 3.5, against 3–4 and 6.8 at 10 m in ERA5, though they should be stronger.
- Antarctica: see O-33. Sea ice: see O-36.
- Across the M2.11.8 runs the T RMS spans 3.2–3.65 K (default 3.46). Less smoothing
  warms the Arctic Ocean ~2 K (O-30) but leaves Tibet 5–6 K colder at ERA5's height
  and doubles the terrain wind error, so the defaults stand.

Heat transport (M2.12) should warm the north and shrink the ice. Clouds, the subtropics'
rain and the winds belong to M2.16's tuning pass; the bright clouds are likely the main
reason the planet is cold.
- After M2.12 (currents on, year 9 of the preset, so not the same spin-up as the year-20
  figures above): global 287.89 K against 288.28 observed; RMS temperature 4.09 K (3.46
  for the same-length run without currents), absorbed sunlight 27.5 W/m², rain still
  +27%, albedo 0.330 against 0.291; 60-70°N -3.6 K, north of 70°N still -7 K. Details in
  O-47.

### O-43: The ocean carries about a quarter of the paper's heat
*Scheduled: M2.16.* On the tilted aquaplanet (642 cells, 20–40 years) the ocean's
transport peaks at 0.7–0.8 PW against Bhatnagar et al.'s 2.25–3.25. Eddies (0.5–0.65
PW against ~0.7) and diffusion (under 0.2) are close; Ekman transport is ~0.4 PW. The
wind stress is 4–5× weaker than Earth's (trades −0.011 to −0.014 N/m², westerlies
+0.015 to +0.026), and Ekman transport is linear in it. With no tilt the equator cools
3–4 K against no transport (paper ~8) and 55–65° warms 3.5 K (paper 5–6), with the ice
edge moving from 57° to 67°. Same at 2,562 cells. Not tuned here: fix the winds first
(O-42).

### O-44: With currents the ocean takes decades to settle
*Watching; M2.13.5, M2.16, M4.5.* The deep layer adjusts with an e-folding time of
about 15–20 years: on the tilted aquaplanet in true time the top of the atmosphere is
still +2 W/m² at year 20 and +0.8 at year 40. Jumping both layers overshoots the
surface, because it extrapolates the mixed layer's fast response. The default with
currents now jumps the deep layer alone (20 years at the ends of years 4 and 7): year 9
is then 293.68 K against 293.75 at true-time year 40, with every 10° band within 0.3 K.
Jumping from year 2, or again in the year after a jump, extrapolates a transient. Open:
tested on one planet at 642 cells, not on Earth or with land; the tropical deep layer
ends up to 1 K cold; and lasting sea ice isn't sped up, since water under ice is pinned
at freezing (no-tilt ice is 20 million km² at year 10 against 12 at true-time year 30).
On the still ocean a deep-only jump cools the deep layer, so that default is unchanged.
On Earth at 2,562 cells the preset's year 9 is still drifting 0.1-0.2 K a year in global
mean surface temperature (top of the atmosphere +0.08 W/m² with gyres, +0.31 without).
A true-time run of ~40 years to check the preset on Earth is deferred to M2.16, after
restart files (M2.13.5) and the performance pass (M2.15).

### O-45: Marginal seas are too small to be basins
*Watching; M2.12.6.* A basin needs 12 interior vertices. Earth at 2,562 cells is
one basin (the Mediterranean has 8 and is dropped); at 10,242 the Mediterranean
and Hudson Bay join it. The Black Sea, Baltic, Red Sea and Caspian never have
enough: they get Ekman, eddy and diffusive transport but no gyre. The world ocean
borders 12–29 landmasses, each a separate coast for the streamfunction. Revisit
when the gyre solver shows how small a gyre it resolves (O-46); an area rule wouldn't
move with resolution.

### O-46: Gyres: straits, the circumpolar current and weak winds
*Watching; winds and gyre depth are for the Earth tuning step (M2.16).* Resolved in
M2.12.7: the island rule replaced Ψ = 0 on every coast, removing the spurious closed
Southern Ocean gyre (27-31 Sv). With the analytic Earth-like wind, at 2,562 and 10,242
cells: Drake Passage 36.9 and 40.9 Sv eastward (observed 130-170; limited by friction and
by how narrow the passage is on the grid; in a test channel a barrier with a 10 degree
gap still passes 90-94% of the open flow), Indonesian throughflow 14.4 and 14.9 (observed
~15), Bering 2.4 and 2.5 (observed ~1), North Atlantic gyre 12.0 and 13.2, North Pacific
20.0 and 21.7. Strength rises 5-10% per doubling of resolution (upwind numerical
friction is about r at 2,562 cells).

Coarse-grid straits, from the half-ocean rule for wet cells: most of the Kuroshio passes
through the Sea of Japan (~18 Sv; real ~2.5), ~2 Sv leaks through Central America, Cook
Strait carries 13 Sv at 10,242 cells, and Madagascar is split at 2,562. In the climate
run the strongest gyre heating anywhere (+64 W/m2) is in the East China Sea from this.

In the Earth climate run (2,562 cells) the gyres are 2-5 times weaker than with the
analytic wind (North Atlantic 4.9 Sv, North Pacific 6.2, Drake ~13-25) because the
model's surface stress is weak (O-43). They carry at most +0.03 PW north and -0.09 PW
south; global mean surface temperature changes by -0.05 K and southern sea ice from 1.45
to 1.14 M km2. Run time +1%.

### O-47: Currents overwarm the Southern Ocean
*Scheduled: M2.16.* Earth at 2,562 cells, year 9, against the 6-year pre-transport run
(`output/e_default`, older code and not settled, so differences under ~1 K are noise).
With Ekman and eddy transport 50-80°S goes from +1-2 K to +6-9 K against ERA5, and
Antarctic sea ice from a 19.7 M km² winter peak to 2.5-3.0 (observed 14.5), summer
minimum 0.2 (observed 2.4). The tropics go from about right to -2 to -2.5 K near the
equator. The ocean carries 0.73 PW south at 39°S, of which eddies are 0.58: not
unrealistic in itself, so the pre-transport Southern Ocean was probably right for the
wrong reason (clouds, O-42). 50-70°N warms 2-3 K, as wanted. The gyres add at most 0.03
PW north and 0.09 south. Untested suspects: the eddy diffusivity (2000 m²/s, the
paper's aquaplanet value), the missing barrier of the circumpolar current to eddy
mixing, the deep-layer jump, and cloud brightness over the Southern Ocean. Tune with the
winds (O-43) at M2.16.
