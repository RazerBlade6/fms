# Design

How each built subsystem works today. Unbuilt subsystems are in
[PLAN.md](PLAN.md), and known limits (O-n) are in
[OBSERVATIONS.md](OBSERVATIONS.md). Numbers are defaults; the code is the
authority.

## Fields and scheduler (`fms-fields`, `fms-sched`)

- A **field** is a typed array with one value per cell, edge or vertex of one
  grid level. Names are `namespace.field` in snake_case, and the registry rejects
  anything else. Metadata: units, valid range, colour hint, buffering, whether
  it's saved, and the quantity (`Scalar`, or `NormalComponent` for an edge field
  holding a vector's component along each edge normal).
- Fields are registered up front. Each has a lock, but a conflicting access
  **panics**: it means a system under-declared its read/write set. Systems reach
  fields only through a context that checks every access.
- A **system** declares its reads, writes and a rate (a multiple of `dt`, or an
  event). The **scheduler** places each due system in the earliest stage after
  every earlier-added system it conflicts with, and runs each stage in parallel.
  The plan doesn't depend on thread count.
- Time is an integer step count. Day, month and year derive from it; fractional
  years alternate like leap years. Each system gets `interval_seconds`, the
  simulated time since it last ran.

Only insolation runs through the scheduler so far. The climate run steps
`LayeredSemiImplicit` directly.

## Grids (`fms-grid`)

- Both levels are subdivided icosahedra optimized into an SCVT by **multilevel
  Lloyd**: a coarse mesh (frequency ≤ 32) is fully optimized, then refined by
  small prime factors and re-optimized. Cells = 10·f² + 2.
- Cells are struct-of-arrays, with CSR neighbours, ordered along a Hilbert curve.
- **Operators** are TRiSK on a C-grid: normal velocities at edges, scalars at
  cells, vorticity at vertices. Vector fields are displayed through Perot's
  cell-centre reconstruction.
- **Fine → coarse**: each fine cell belongs wholly to the coarse cell containing
  its centre, found by a greedy Delaunay walk and precomputed. Intensive fields
  take area-weighted means; extensive ones sum.
- **Coarse → fine**: barycentric interpolation on the coarse Delaunay triangle.

## Planet (`fms-planet`)

- Orbit from Kepler's equation; insolation from latitude, declination and hour
  angle, instantaneous and daily mean. Stellar flux is `1361 W/m² · L / a²`.
- `t = 0` is the northern vernal equinox, at midnight at longitude 0. +z is the
  north pole, +x is (0°, 0°). Rotation is prograde:
  `1/T_sid = 1/T_sol + 1/T_orbit`.
- Seasons change at the true equinoxes and solstices, so they differ in length
  on an eccentric orbit.

## Terrain (`fms-worldgen`, `fms-sim`)

Heightmaps resample onto the fine grid as area-weighted means. Each climate cell
gets the mean and spread of surface height, the mean slope and the land
fraction. A terrain cell is sea if it is at or below zero and its connected body
of such cells is the world ocean (the largest) or covers at least 0.005% of the
sphere (about 25,500 km² on Earth); smaller bodies below sea level are land. The
Caspian and Black Sea stay water; the Qattara Depression and the Dead Sea are
land (O-11). Surface height is `max(elevation, 0)`. The atmosphere sees it after
ten conservative diffusion passes (O-30).

## Dynamical core (`fms-atmos`)

**Vertical.** σ layers (`σ = p/p_s`, top at p = 0) at any interfaces, up to 16.
The climate uses σ 0.25/0.5/0.85: four layers, the bottom ~1.4 km thick.
Differencing follows Arakawa & Suarez (1983), which conserves total energy
exactly in space and keeps an isentropic atmosphere at rest over mountains.

**Prognostics.** Surface pressure and mass-weighted θ and specific humidity at
cells; normal wind at edges, per layer.

**Time stepping.** Semi-implicit leapfrog with the RAW filter (ν 0.05, α 0.53).
Gravity-wave terms, linearized about an isothermal 300 K reference (or warmer if
the start is), are treated implicitly. They split by vertical mode into one
Helmholtz problem per layer, solved by Jacobi-preconditioned conjugate gradient
to 1e-8. The step comes from a design maximum wind (120 m/s in the climate run)
at Courant 0.6, rounded to divide a day. Forcing and diffusion are lagged a step.

**Dissipation.** ∇⁴ hyperdiffusion on winds and θ at the strongest strength the
step allows with a 2× margin (`τ` about 5 h at 20k cells). Kinetic energy that
diffusion, drag or mixing removes returns as heat.

**Humidity transport.** Flux-corrected transport (Zalesak 1979) per layer, using
the step's own mass fluxes. Water is conserved to round-off; the filter's trace
of negative humidity is refilled from the column.

**Safety.** `safety::check` stops a run at the first non-finite value, Courant
violation or out-of-bounds pressure or θ, and writes a diagnostic dump.

## Column physics (`column_physics`)

Surface fluxes and tiles update every step. Radiation runs every 3 h, with
sunlight averaged over the coming interval and the latest moist physics' clouds.
The atmosphere and tiles share the same stored fluxes and precipitation, so
energy and water budgets close by construction, up to the time-stepping residual.

`SingleColumn` is a fixed-humidity radiative–convective column for checking the
radiation against Earth's global mean. `column_cases` runs BOMEX and DYCOMS in
milliseconds.

## Radiation (`radiation`)

**Longwave.** Two-stream through the layers over a black surface, in grey
intervals:
- CO₂: an opaque band at 667 cm⁻¹ whose width grows with the log of the
  pressure-broadened amount (band-wing model, `l` = 30 cm⁻¹).
- Window, 800–1200 cm⁻¹: self-broadened water continuum.
- Water band cores and edges: a four-term k-distribution each, on the
  pressure-weighted water column.

Each layer's emission is linear in optical depth between its edges. An interface
takes the temperature of the layer below along the lapse rate; the top layer's
top is at its own mean temperature.

**Shortwave.** Two bands through the adding method, at `μ` = 0.5: the visible
(ozone at the top) and the near-infrared (35% of sunlight; water, with
transmissivity `exp(−a√u)` on the cumulative slant path, so absorption doesn't
depend on the layering). Rayleigh scattering is one non-absorbing 6% sheet at
σ = ½. Clouds and the surface reflect band by band.

**Tiles.** Radiation sees one surface per cell (mean albedo per band, mean σT⁴)
and splits absorbed sunlight by `1 − α` per band and emission by `T⁴`.

## Moisture and clouds (`moisture`)

Within a layer, relative humidity has the Manabe–Wetherald shape scaled to hold
the layer's mean humidity.

- **Convection** (Betts–Miller, 2 h): a surface-air parcel is lifted along the
  pseudo-adiabat to give a reference profile. If a layer above the lowest is
  colder than its reference, every layer from the highest such one down relaxes
  towards it. Deep convection rains and heats by exactly the latent heat. Shallow
  convection redistributes heat among the cloud layers and moves water only from
  layers above their reference humidity to those below.
- **Large-scale condensation** (4 h) above layer humidity thresholds of 0.45,
  0.7 and 0.9 at σ-centres ¼, ¾ and 0.925, interpolated between. Falling rain
  re-evaporates in proportion to each layer's dryness.
- **Rain or snow** by surface air temperature, a ramp over 271.15–275.15 K.
- **Clouds**: two sheets at σ 0.3 and 0.9, each from the relative humidity of the
  layer it lies in. Cover ramps linearly (0.45–0.85 high, 0.5–1.0 low); deep
  convection adds anvils to the high sheet, `0.15 √(mm/day)` up to 0.6. Low
  cloud thins in air below 3 g/kg, in proportion down to 0.15 of its cover
  ("freeze-dry", Vavrus & Waliser 2008). Albedo 0.4 and 0.5; longwave
  emissivity 0.8 and 1.

## Boundary layer (`surface_layer`)

The bottom layer is well mixed. The 10 m wind follows a log profile below it.
Stress, sensible heat and evaporation use bulk formulas at 10 m with neutral
coefficients from the roughness lengths, Charnock roughness over open water,
Louis stability in unstable air and `1/(1 + 10 Ri)` in stable air, and a gust
from `w*`.

The layer above is entrained at `w_e = A w*³/(h Δb)` (A = 0.2), exchanging dry
static energy, water and momentum. Where surface air saturates below the top and
stays buoyant, shallow cumulus carry water and momentum up with a mass flux of
`0.03 ρ w*`.

## Surface tiles and land (`column_physics`)

Each cell has ocean, sea-ice and land tiles, weighted by land fraction and ice
cover, each with its own temperature, albedo, roughness and fluxes.

- **Land energy.** A skin with no heat capacity over a 1 m soil store, conductance
  7 W/m²/K, albedo 0.2, roughness 0.1 m. Each step a bracketed secant solve
  balances the skin implicitly with the soil. The soil takes exactly the net flux
  the atmosphere was given.
- **Longwave lag.** Land emits at its current skin temperature; the difference
  from the last radiation call heats the bottom layer. On a radiation step the
  skin is solved just before the new sunlight arrives.
- **Soil water.** A 150 mm Manabe bucket, evaporating freely above 75% full.
  Overflow runs off and leaves the planet's water until rivers (M3).
- **Snow.** A pack per land tile, brightening land to 0.65 at 10 kg/m². A skin
  reaching freezing stays there and its surplus melts the pack into the bucket.
  Evaporation takes snow first. The pack counts as a heat deficit, `−L_f × snow`.
- **Ice sheets.** Land fully snow-covered for a whole year becomes a sheet (albedo
  0.82, roughness 1 mm), and reverts when cover is lost. Land starting below 267 K
  starts as a sheet under 1000 kg/m² (O-39).
- **Sub-grid hills.** Upslope wind rains out a share of its moisture (Smith 1979,
  `P = ε ρ q w`), with `ε` reaching 0.5 at 500 m of height spread and only above
  50% bottom-layer humidity. The spread also adds a drag (a simplified Lott &
  Miller) whose energy heats the boundary layer.

## Ocean (`fms-ocean`)

A 50 m mixed layer over a 150 m deep layer in each cell with ocean; it knows
nothing about the atmosphere. Heat moves by:
- Ekman transport: the wind stress on open
  water drives a mass flux in the mixed layer across each edge between ocean
  cells. With `w = exp(−f²/ε²)` and `ε` = 10⁻⁵ s⁻¹ it is the Ekman transport
  `(1 − w)/f · (−k×τ)`, a frictional flow `(w/ε) τ`, and near the equator a
  northward Sverdrup transport `w B (k·∇×τ)`, where `B = β/(β² + (ε/R)²)` is
  `1/β` on Earth and vanishes on a planet that doesn't rotate. The curl is
  taken at vertices whose three cells all have ocean. The deep layer carries
  the same mass back, and each cell's divergence wells up from below or sinks.
  Temperatures are carried upstream, in sub-steps that keep them within their
  range, and heat is conserved to round-off. `--no-ocean-sverdrup` gives
  `M = (ετ − f k×τ)/(ε² + f²)` instead;
- eddy transport, after Gent & McWilliams:
  a second mass flux in the mixed layer, `−ρ κ s` with `κ` = 2000 m²/s, where
  `s` is the slope of the isotherms: the gradient of the two layers' mean
  temperature over their difference across 500 m, floored at 0.25 K. The slope
  is limited to 0.002. It is returned and advected together with the Ekman
  flux, upstream of their sum;
- horizontal diffusion in both layers at 2000 m²/s, as fluxes across edges, each
  opened by the smaller of its cells' ocean fractions;
- gyres, below: the wind-driven flow round each basin, in both layers;
- convective adjustment when the mixed layer is colder than the deep layer.

The diffusivity and `κ` are Earth's. On another planet both are multiplied by
Earth's rotation rate over the planet's, at most 4.

The deep layer starts at the mixed layer's temperature.
`OceanParams::still()` is the ocean without currents (`--no-ocean-ekman
--no-ocean-eddies --no-ocean-gyres`), whose deep layer starts no warmer than 4 °C.
`OceanParams::slab(10.0)` is the old motionless 10 m slab, kept for pinned
tests.

The climate run writes the northward heat transport of each of these (Ekman,
eddies, gyres, diffusion) across each latitude, in PW, from every step of the final
year. It also maps the gyre streamfunction and the heating each gyre gives
(`ocean.gyre_streamfunction`, `ocean.gyre_heating`, `ocean.gyre_mixed_heating`,
`ocean.gyre_heat_gained`).

**Basins**, for the gyres: a cell is wet if at least half ocean, a vertex
is coast if any of its three cells isn't wet, and a basin is a connected set of
the other vertices, kept if it has at least 12. An aquaplanet is one basin with
no coast.

**Gyres** (`fms-ocean/gyres`). A steady Stommel balance for the mass-transport
streamfunction Ψ (kg/s) on the vertices, `r∇²Ψ + β ∂Ψ/∂x = curl τ`. Island rule:
each landmass's Ψ is one unknown constant, and its equation is the sum of its
vertices' equations (friction and β fluxes round it balance the wind's circulation
round it). The largest landmass is the reference, Ψ = 0; basins too small to keep
join the landmass they border; everything connected is one linear system, so the
Southern Ocean has a circumpolar current. An ocean with no land (aquaplanet) has no
gyre. The β term is in flux form (the flux across a dual edge is the difference of f
between its two cells) and taken upwind, which keeps the western boundary current
about one vertex spacing wide and adds a numerical friction of about βΔx/2. Friction
r = 2e-6 s⁻¹. β follows the signed rotation, so a retrograde planet's boundary
currents are on eastern coasts. The flow across an edge is the difference of Ψ
between its two vertices, so it is exactly non-divergent and zero at coasts. The
matrix and its ILU(0) factors are built once; each solve is BiCGSTAB, warm-started.
In the climate run the solver sees a 30-day running mean of the wind stress and
re-solves once a day. The flow is spread at one velocity over 700 m, so the mixed
and deep layers carry 50/700 and 150/700 of it in the same direction. It is
advected upstream together with Ekman and eddies, each layer upstream of its own net
flux. The gyre adds no vertical exchange: Ekman's flux returns in the deep layer and
sums to zero over depth, so the two don't count the same water twice.

## Sea ice (`fms-ocean/ice`)

Thermodynamic, after Codron (2012) and the authors' code. Each cell has one
fraction, mean thickness and surface temperature, with a linear profile down to
the freezing point (271.35 K) at the base.

- **Heat** changes every step; **area and thickness** change only at radiation
  calls, so stored sunlight stays balanced.
- **Ocean to ice**: `0.0057 ρ c u* (T − T_f)`, with `u*` from the wind stress
  over open water.
- **Basal growth** thickens ice to 1.5 m, then spreads it, more as it nears the
  10 m limit. **Melting** thins it to 1.5 m, then shares between thinning and
  area loss down to 0.15 m. Open water below freezing forms ice 0.15–1 m thick.
- **Snow** above 15 kg/m² is its own layer, solved with the interface temperature;
  below that the two temperatures are one. Wind blows a share (0.4 × open
  fraction) into the water. Surface melt takes snow, then ice, then warms the
  ocean. Snow below the waterline floods into ice.
- **Albedo** per band rises with thickness from the water's towards 0.57 visible
  and 0.31 near-infrared over a 0.3 m scale, blending to snow's 0.65 with cover.

## Climate run and spin-up (`fms-cli`)

`fms-climate climate` starts at rest from a radiative–convective equilibrium in
hydrostatic balance with the terrain. The ocean starts shaped by the year's
sunlight. By default the ocean is **jumped** ahead. With currents the deep layer alone
jumps, at the ends of years 4 and 7: each repeats that year's change for twenty
more, moving no temperature over 3 K, and two true-time years follow (O-44).
Without currents both layers jump at the ends of years 2–4, by three years each
(O-37). Pinned tests turn jumps off. `--ocean-jumps 0` is the true-time check.
