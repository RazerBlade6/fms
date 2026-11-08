# Design Plan

This is the agreed design for the Fantasy Map Simulator: a desktop app for generating,
editing and simulating fantasy planets. It covers terrain, water, climate and weather
over time. The app is built so it can grow into a civilisation simulator later.

## 1. Goals

- **Customizable planet**: radius, gravity, rotation speed, axial tilt, year length,
  orbital eccentricity, star luminosity, atmospheric pressure, CO₂ (greenhouse
  strength) and ocean fraction.
- **Terrain editor**: sculpt elevation, carve rivers, create lakes and seas.
- **Map modes**: climate (Köppen), biomes, heightmap, contours, wind, temperature,
  precipitation, ocean currents, sea ice, plates. Any field in the field registry
  can be shown as a map mode.
- **Auto-generated climates**: climate zones are not assigned by hand. They come
  out of the simulation and follow from the planet's physical setup.
- **Time simulation**: pressing Play advances time. Speed 1 = 1 sim-day per second.
  Speed 5 = the fastest the model can sustain (target is about 1 sim-year per second).
  A "day" is the planet's own solar day (one day/night cycle), not an Earth day.
- **Extensibility**: a future civilisation layer (soil fertility, aridity, agents)
  plugs in without changing the core.

## 2. Key decisions

| Area | Decision | Rationale |
|---|---|---|
| Language | Rust | rayon gives safe parallelism over cells, Cargo is simple, and ECS crates exist for the civ layer |
| Platform | Desktop, cross-platform (macOS first) | wgpu targets Metal, DX12 and Vulkan |
| UI and rendering | `eframe`/`egui` panels + custom `wgpu` renderer | Full control; the sim stays independent of any engine and runs headless |
| Geometry | Simulation runs on a sphere; flat maps are projections of it | Tilt, Coriolis and latitude need a real sphere; no polar distortion |
| Sphere grid | Icosahedral geodesic grid (Voronoi hexagons + 12 pentagons) | Cells are nearly uniform with no pole singularity; standard for finite-volume models (ICON, MPAS) |
| Climate fidelity | Simplified GCM ("GCM-lite"): actual fluid dynamics | Weather emerges from the physics instead of being scripted |
| Atmosphere | Layered primitive-equation model: any number of σ layers; two until M2 step 10, then three with a thin boundary layer | Two layers are the smallest model that produces Hadley and Ferrel cells, jets and storms on its own; a third resolves the surface wind and the boundary layer's water |
| Ocean | Two-layer dynamical slab (Ekman/Sverdrup, eddy and diffusive heat transport) + wind-driven gyres + thermodynamic sea ice | Ocean heat transport (western boundary currents) matters for coastal climates |
| Compute | CPU with rayon first, GPU (wgpu compute) later | Easier to debug the physics; data is laid out as struct-of-arrays so it can be ported |
| Time stepping | Semi-implicit GCM; base step derived from the CFL limit (~20–30 min at default resolution); slower subsystems run every N steps | Stable and affordable; resolves the day/night cycle |
| Day length | "1 day" = the planet's solar day (one day/night cycle) | Speed 1 always shows one day/night cycle per second |
| High speeds | Speed tiers: full GCM at speeds 1–2, cheaper climate mode at 3–5 | Speed 5 cannot afford full dynamics |
| Resolution | Two levels: fine terrain (~1–2.5M cells) and coarse climate grid (~10–40k cells) | Detailed terrain and rivers while keeping the GCM affordable |
| World gen | One-shot plate layout now; time-stepped plate sim later behind the same interface | Plausible mountain ranges in under 1 s; can be upgraded later |
| Geology during play | Terrain is static; erosion is an editor tool | Climate can settle, and a future civ layer gets a stable world |
| Hydrology | Mostly emergent, with manual overrides | Rain feeds runoff that is routed into rivers; depressions become lakes; painted rivers carve the terrain |
| Live edits | Edits apply while the sim runs, then the climate catches up via spin-up | Sandbox feel |
| Classification | Köppen–Geiger + Whittaker-style biomes | Familiar categories; both are derived fields |
| Planet range | Earth-like range + eccentricity; extremes clamped with warnings | Keeps the GCM numerically stable |
| Edit views | Brushes in both globe and flat views, acting on sphere cells | No projection distortion when painting |
| Determinism | Seeded RNG, fixed timestep, deterministic reductions; versioned binary saves | Reproducible worlds, exact save/load mid-run, reproducible bugs |
| Portable math | Transcendental functions (`sin`, `exp`, `atan2`, …) come from the pure-Rust `libm` port via `fms-math`, never std; enforced by `clippy.toml` | std calls the platform C library, which differs between macOS, Windows and Linux; without this, saves would diverge across machines |
| Undo | Command-based; each stroke stores a sparse diff of the cells it changed; terrain only | Cheap; the climate re-settles after undo instead of rewinding |
| Validation | Earth benchmark + idealized test planets, run as regression tests | The only honest way to tune a simplified GCM |
| First milestone | Headless GCM first | Prove the physics before building UI on top of it |

## 3. Architecture

### 3.1 Core concept: field registry + system scheduler

- A **field** is a named, typed array with one value per cell, edge or vertex of
  one grid level (`terrain` or `climate`). Edges carry normal velocities and
  fluxes; vertices carry vorticity. Examples: `terrain.elevation`,
  `atmos.surface_temp`, `atmos.wind_u_lower`, `climate.koppen`, and later
  `soil.fertility` and `civ.population`.
- **Naming convention:** `namespace.field`, both parts `snake_case`. The namespace
  is the owning system or crate (`terrain`, `atmos`, `ocean`, `hydro`, `climate`,
  `soil`, `civ`, …), so independently written crates can't collide. The registry
  rejects names that don't follow it.
  Types: f32, f64, u8, u16, u32, i32. Metadata: units, valid range (also used to
  detect blow-ups), colour-scale hint, single/double buffering, whether the
  field is saved, and the quantity. The quantity is `Scalar`, or for edge fields
  `NormalComponent`: a vector field stored as its component along each edge
  normal. Those are displayed by rebuilding cell-centre vectors with Perot's
  method.
- Fields are registered up front. After that, different fields can be read and
  written from many threads at once. Each field has a lock, but a conflicting
  access **panics** rather than blocks: it means a system under-declared its
  read/write set, which the scheduler is supposed to prevent.
- A **system** declares:
  - which fields it reads and which it writes
  - its tick rate as a multiple of the base step `dt`, or as an event (daily,
    monthly, on edit); see §4.6
  - which speed tiers it runs in
- The **scheduler** guarantees that every step gives the same result as running
  the due systems one after another **in the order they were added**. For
  parallelism, each due system is placed in the earliest stage after every
  earlier-added system it conflicts with (one writes a field the other reads or
  writes). Systems within a stage run in parallel. The plan doesn't depend on
  thread count, so results are bit-identical on 1 or 64 threads, and (with
  `fms-math`) on macOS, Windows and Linux. A golden-hash test in CI checks this.
- Systems access fields only through a context that checks every access
  against the declared read/write set. An undeclared access panics on the first
  try, not only when it happens to race with another system.
- Time is an integer step count, so it never drifts. Day, year, month and
  calendar quarter are derived from it. Fractional years alternate in length like
  leap years. Year fraction 0 is the northern vernal equinox by convention. Each
  system gets `interval_seconds`, the simulated time since it last ran, as its
  timestep.
- Seasons are astronomical and come from the orbit (`fms-planet`), not the
  calendar. They change at the equinoxes and solstices, so on an eccentric orbit
  they differ in length: Earth's northern spring, summer, autumn and winter are
  92.76, 93.65, 89.84 and 89.00 days.
- **Map modes** are renderers over fields: a colormap, contours, vectors or
  particles. A new field gets a basic map mode for free.
- **Civ layer (future)**: an ECS (`hecs` or `bevy_ecs`) for agents such as
  settlements and populations. Agents read and write fields through the same registry.
  Cells themselves are **not** entities; with over a million cells, an ECS would hurt cache performance.

### 3.2 Grids

- Both levels are subdivided icosahedra, optimized into a spherical centroidal
  Voronoi tessellation (SCVT) with **multilevel Lloyd**: a coarse mesh
  (frequency ≤ 32) is fully optimized, then refined by small prime factors and
  re-optimized at each level. Plain Lloyd needs roughly O(frequency²) iterations
  and is impractical above ~100k cells. Multilevel needs ~600 iterations per
  level: about 0.4 s for 20k cells and about 8 s for 1.6M cells (release build).
  Pick frequencies with small prime factors (e.g. 45, 400); cells = 10·f² + 2.
  Fine meshes will be cached to disk once persistence exists.
- Known limitations of the grids (cell-size spread, frequency constraints,
  pentagon errors) are tracked in [OBSERVATIONS.md](OBSERVATIONS.md).
- The two meshes are optimized independently, so they are not exactly nested.
  The mapping between them is found by point location (a greedy walk over the
  Delaunay graph) and precomputed once (~30 ms for 1.6M → 20k).
- **Fine → coarse**: each fine cell belongs wholly to the coarse cell containing
  its centre (~80 fine cells per coarse cell at 1.6M / 20k). Intensive fields use
  an area-weighted mean over the children; extensive fields (water, mass) use a sum,
  which conserves global totals exactly. Sub-grid terrain statistics (land
  fraction, mean and standard deviation of surface height, mean slope) are kept
  so the climate grid still feels mountains it can't resolve. Surface height is
  the land surface, or sea level over ocean.
- **Coarse → fine**: barycentric interpolation on the coarse Delaunay triangle,
  then corrected for altitude using the lapse rate for display and for hydrology.
- Cells are stored as struct-of-arrays, with neighbours in compressed sparse-row
  (CSR) adjacency lists. Cells are ordered along a space-filling curve so neighbours
  sit close together in memory.

### 3.3 Workspace layout

```
crates/
  math/       portable transcendental functions (libm), for cross-platform determinism
  geo/        icosahedral grids, two-level mapping, projections, picking
  fields/     field registry, SoA storage, serialization
  sched/      system scheduler (tick rates, R/W sets, determinism)
  planet/     planet parameters, orbit, insolation
  worldgen/   Generator trait, plate generator, noise, heightmap import
  atmos/      layered GCM (dynamics + physics parameterizations)
  ocean/      slab mixed layer, wind-driven currents, sea ice
  hydro/      runoff, flow routing, depression filling, lakes, basins
  climate/    rolling monthly stats, Köppen, biomes, high-speed climate mode
  sim/        orchestration, speed tiers, spin-up, save/load
  cli/        headless runner, PNG/field dumps, Earth benchmark
  app/        eframe/egui + wgpu renderer, brushes, map modes, inspector
```

Package names have an `fms-` prefix (`fms-geo`, `fms-sim`, …) to avoid clashing
with crates.io names. The binaries are `fms` (cli) and `fantasy-map` (app).

Dependencies point one way. The `app` crate depends on the sim crates; the sim
crates never depend on `app`. Everything except `app` runs headless.

## 4. Subsystems

### 4.1 Planet and insolation
- Orbit is computed from eccentricity and year length (Kepler's equation). Insolation
  at the top of the atmosphere depends on latitude, the sun's declination and the
  hour angle; the daily mean is used in the cheap climate mode.
- Rotation rate sets the Coriolis parameter, which in turn controls how many
  circulation cells form and how wide they are.
- Greenhouse strength is a parameterized function of CO₂ and pressure that sets
  longwave opacity.
- **Parameters and supported ranges.** Values outside these are clamped with a
  warning. NaN or infinite values are replaced by Earth's value.

  | Parameter | Range | Earth |
  |---|---|---|
  | Radius | 1,000–20,000 km | 6,371 km |
  | Surface gravity | 1–30 m/s² | 9.81 |
  | Solar day | 4–240 h | 24 h |
  | Year | 20–5,000 solar days | 365.2422 |
  | Obliquity | 0–60° | 23.44° |
  | Eccentricity | 0–0.3 | 0.0167 |
  | Solar longitude of perihelion | wrapped to 0–360° | 282.9° |
  | Star luminosity | 0.05–20 L☉ | 1 |
  | Semi-major axis | 0.05–50 AU | 1 |
  | Surface pressure | 0.1–10 bar | 1.013 |
  | CO₂ | 0–100,000 ppm | 420 |
  | Ocean fraction | 0–1 | 0.71 |

  Stellar flux `S₀ = 1361 W/m² · L / a²`. If it falls outside 500–3,000 W/m², the
  planet is allowed but flagged as likely frozen or runaway.
- **Conventions.** `t = 0` is the northern vernal equinox, at midnight at
  longitude 0. That matches the calendar's year fraction 0. The year length is
  the orbital period in solar days. Rotation is prograde, so the sidereal day
  comes from `1/T_sid = 1/T_sol + 1/T_orbit`. The Sun's position uses mean
  solar time; the equation of time is ignored.
- Known limitations of the planet model (Kepler's third law, equation of time)
  are tracked in [OBSERVATIONS.md](OBSERVATIONS.md).

### 4.2 Atmosphere (layered GCM)
- Primitive equations solved with a finite-volume method on the icosahedral grid.
  Each layer has wind, temperature and moisture (moisture mostly in the lowest
  layers).
- Physics parameterizations:
  - shortwave and longwave radiation (semi-grey; see below)
  - surface fluxes (bulk formulas)
  - large-scale condensation and simple convective adjustment
  - orographic lift using the sub-grid terrain statistics
  - boundary-layer friction
- **Vertical structure** (decided at M2.5, generalized at M2.9.5): sigma
  layers (`σ = p/p_s`, model top at p = 0) at any interfaces
  (`crates/atmos/src/levels.rs`). Two equal-mass layers (interface at σ = 0.5)
  are the classic Phillips two-level model, the climate's setting until M2
  step 10 moves it to three with a thin bottom layer. Prognostic: surface
  pressure, the wind in each layer, and potential temperature per layer
  (stored mass-weighted, `p_s θ`, so total mass and total θ are conserved to
  round-off). Vertical differencing follows Arakawa & Suarez (1983): layer
  Exner functions are mass-weighted layer means of `(p/p₀)^κ`, and the
  hydrostatic relation uses the same interface θ (the mean of the two
  layers') as vertical advection. Together these conserve total energy
  (kinetic + enthalpy + surface potential) exactly in the spatial
  discretization for any layers, and keep an isentropic atmosphere at rest
  over mountains exactly at rest. With more than one interface the
  arithmetic interface θ makes the vertical modes slightly non-normal
  (O-28); the semi-implicit solver decomposes them in general.
- **Radiation** (decided at M2.8): two-stream transfer through the layers
  and the surface, semi-grey in four kinds of spectral interval. CO₂ has an
  opaque band at 667 cm⁻¹ whose width grows with the log of the
  pressure-broadened CO₂ amount (the band-wing model), so its forcing is
  logarithmic by construction. The 800–1200 cm⁻¹ window has a self-broadened
  water continuum. Water's band cores and their weak edges each use a
  four-term k-distribution. Each layer's emission is linear in optical depth
  between its edges, and each interface takes the temperature of the layer
  below along the convective lapse rate.
  Water vapour comes from a fixed relative-humidity profile (Manabe &
  Wetherald 1967) until humidity is prognostic (M2 step 9). Shortwave, in
  two bands through the adding method: ozone in the visible, and
  near-infrared water-vapour absorption saturating with the cumulative
  slant path from the top (so it doesn't depend on the layering), a stand-in
  cloud deck at the layer interface, and the surface albedo. Calibrated in a
  single column against Earth's global-mean column; see `crates/atmos/src/radiation.rs`
  and [OBSERVATIONS.md](OBSERVATIONS.md) O-19, O-20.
- **Moisture** (decided at M2.9): each layer carries its mean specific
  humidity (mass-weighted, like θ), moved by flux-corrected transport
  (Zalesak 1979) around the leapfrog step, so it stays positive and water is
  conserved to round-off. Within a layer, relative humidity has the
  Manabe–Wetherald shape scaled to the layer's water, and humidity crosses
  the interface at the profiles' values there. Physics: Betts–Miller
  convection relaxing towards the moist adiabat of a surface-air parcel
  from the highest unstable layer down (deep, raining, or shallow, moving
  water only from layers above their reference to those below), large-scale condensation
  above calibrated layer humidities, re-evaporation of falling rain, rain or
  snow by surface air temperature, and evaporation from the slab by bulk
  formula. Two cloud sheets (high and low) from the layers' relative
  humidity plus convective anvils, acting on sunlight and longwave. See
  `crates/atmos/src/moisture.rs` and [OBSERVATIONS.md](OBSERVATIONS.md)
  O-22–O-26.
- **Semi-implicit time stepping** (as in SPEEDY): fast gravity waves are solved
  implicitly, so the timestep is limited by advection (winds) rather than wave speed.
  At the default resolution this gives a base step of ~20–30 min. See §4.6.
  The linearization is about an isothermal reference at 300 K (or warmer, if
  the start is): warmer and more stable than the atmosphere, as stability at
  long steps needs (Simmons, Hoskins & Burridge 1978).

### 4.3 Ocean
(Decided at M2.9.) A dynamical slab ocean after Codron (2012) and Bhatnagar
et al. (2026), plus a wind-driven gyre:
- **Two layers** per ocean cell: a mixed layer (50 m) exchanging heat with the
  atmosphere, over a deep layer (150 m), both prognostic temperatures. A
  fixed mixed-layer depth to start (regional depths are a later refinement).
- **Heat transport**, all as mass or heat fluxes across cell edges so the
  ocean's budget closes exactly: Ekman transport from the surface wind
  stress (frictional, blending into Sverdrup near the equator) with the deep
  layer's return flow and the vertical exchange it implies; Gent–McWilliams
  eddy transport along the isotherm slope; horizontal diffusion; convective
  adjustment. Coefficients follow the paper for Earth's rotation and scale
  with the deformation radius on other planets.
- **Gyres**: a barotropic Stommel–Munk streamfunction from the wind-stress
  curl over each basin, for western boundary currents (which shape coastal
  climates) and the ocean-current map. Beyond the paper, which has none.
- **Sea ice**: thermodynamic (no drift); a prognostic surface temperature with
  a linear profile, snow on ice, and a thickness-dependent albedo in visible
  and near-infrared light. Salinity only sets the freezing point; no
  density-driven overturning.

### 4.4 Hydrology (fine grid)
- Precipitation minus evaporation gives runoff, which is routed downhill along
  steepest-descent paths. Accumulated flow above a threshold becomes a river.
- Depressions are filled (priority-flood). A filled basin becomes a lake when
  inflow exceeds evaporation. Basins with no outlet produce salt lakes or pans.
- Painted rivers carve channels into the elevation field, so the result stays
  physically consistent.

### 4.5 Climate classification
- Rolling multi-year means of monthly temperature and precipitation per cell.
- Köppen–Geiger classification (Af, BWh, Cfb, …) and a Whittaker biome layer.
  Biomes drive surface albedo and roughness, which feed back into the GCM.

### 4.6 Time stepping and speed tiers

**Definition of a day.** "Day" always means the planet's own **solar day** (one
day/night cycle), not an Earth day. Speed 1 = one planetary day per second, at
any rotation rate.

**Base step (`dt`).** The scheduler runs one base step, and every system runs
every N base steps.
- `dt` is derived from the climate grid spacing, the planet radius and the maximum
  expected wind speed, using an advective CFL limit with a safety margin.
- It is then rounded down so that it divides the solar day evenly, and capped so
  that there are at least 24 steps per solar day (keeps the day/night cycle resolved).
- `dt` stays fixed for a given world. It only changes when the resolution or
  planet parameters change, which keeps saves and replays deterministic.
- Speed controls change how many base steps run per real second, never `dt`
  itself.

**System rates** (typical values at the default ~20k-cell climate grid on an
Earth-sized planet; each system declares its rate in its scheduler registration):

| System | Rate |
|---|---|
| Atmospheric dynamics (winds, pressure, heat/moisture transport) | every base step (~20–30 min, ~48–72 per day) |
| Surface fluxes, clouds and precipitation | every base step |
| Radiation | every ~3 h |
| Ocean mixed layer, currents, sea ice | every 1–6 h |
| Hydrology (runoff, flow routing, lakes) | daily |
| Climate statistics | accumulated daily; Köppen/biomes reclassified monthly |
| Future civ systems | daily / monthly / yearly |

Rough throughput: speed 1 at `dt` = 30 min is 48 steps/s, which is trivial. Speed 2
(~1 week/s) is about 340 steps/s, which a 2-layer model on 20k cells should sustain
on multiple CPU cores. Faster speeds use climate mode.

**Speed tiers.**

| Speed | Target | Mode |
|---|---|---|
| 1 | 1 day/s | Full GCM |
| 2 | ~1 week/s | Full GCM |
| 3 | ~1 month/s | Climate mode |
| 4 | ~3 months/s | Climate mode |
| 5 | ~1 year/s (best effort) | Climate mode |

**Climate mode**: during full-GCM runs and spin-up, the model learns a per-cell
seasonal climatology (means and variances). At high speeds the sim samples weather
from that climatology with random variation instead of solving the dynamics.
Switching back to speeds 1–2 restarts the full GCM from a state consistent with
the climatology.

### 4.7 Spin-up and live edits
- On world creation, load, or after edits, an accelerated spin-up runs a few
  sim-years in the background (oceans set the timescale). The UI shows progress.
- Edits apply right away. Fields that depend on terrain are marked dirty: the
  grid mapping, hydrology and sub-grid terrain statistics are recomputed for the
  affected area, then spin-up starts.

### 4.8 World generation
- `Generator` trait: given planet parameters and a seed, it produces fields
  (`elevation`, `plate_id`, `crust_type`, `crust_age`, …).
- **v1: one-shot plate layout**
  - Split the sphere into Voronoi plates. Each plate gets a rotation pole, an
    angular velocity, and a continental or oceanic type.
  - Classify each boundary from the plates' relative motion:
    - convergent: mountain ranges, subduction trenches, island arcs
    - divergent: rifts and mid-ocean ridges
    - transform
  - Elevation comes from distance to boundaries plus layered noise. Continental
    shelves are added, and sea level is set from the ocean fraction.
  - Runs in under 1 s.
- **Later: time-stepped plate sim** (in the style of Cortial et al. 2019,
  *Procedural Tectonic Planets*). It implements the same trait and outputs the
  same fields.
- Also: a blank-ocean start, and heightmap import (equirectangular PNG/GeoTIFF).
  The Earth benchmark uses heightmap import.

### 4.9 Editor
- Brushes in both the globe and flat views. Mouse positions are converted to cells
  on the sphere, so projections never distort a stroke.
- Tools: raise, lower, smooth, flatten, noise, erode (hydraulic and thermal),
  set sea level, paint river, fill lake.
- Undo/redo: each stroke is a command holding a sparse before/after diff of the
  cells it changed.

### 4.10 Visualization
- Globe (orbit camera) and flat projections (equirectangular first, more later).
- Map modes: heightmap, contours, temperature, precipitation, wind (animated
  particles), ocean currents, sea ice and snow, clouds and precipitation overlay,
  Köppen, biomes, plates, plus any registry field.
- HUD showing date, season and speed; Play/Pause; speeds 1–5.
- Cell inspector: click a cell to see its current values and annual climate graphs
  (temperature and precipitation by month).

### 4.11 Persistence
- Versioned binary save format holding planet parameters, all persistent fields,
  sim clock, RNG state and the full undo stack. Fields are tagged by name and
  version so older saves can be migrated.
- Determinism contract: same save + same build = bit-identical sim, on any
  thread count and on any of macOS, Windows and Linux. Parallel reductions use
  a fixed order, and transcendental math goes through `fms-math` (`libm`).

## 5. Validation

- **Analytic tests**: insolation by latitude and season; energy balance on a
  planet with no atmosphere.
- **Idealized planets**:
  - aquaplanet: a symmetric Hadley circulation forms
  - no-tilt planet: no seasons
  - fast rotation: more, narrower circulation cells
  - slow rotation: the Hadley cell widens
- **Earth benchmark**: real topography with Earth parameters. Checks that
  - trade winds blow east to west
  - mid-latitude jets form
  - the tropical rain belt migrates with the seasons
  - the Sahara and Atacama are deserts and the Amazon and Congo are rainforests
  - monsoons appear
  - Köppen agreement with observed maps exceeds a target threshold
- **Determinism tests**: running twice from the same save must give
  bit-identical results.
- These tests run before each push as regression checks.

## 6. Equations and numerical methods

Names only, as a reading list. Grouped by subsystem.

### 6.1 Orbit, rotation and insolation
- Kepler's equation (mean anomaly → eccentric anomaly), solved by Newton iteration
- Polar equation of an ellipse (orbital distance from true anomaly)
- Kepler's third law (orbital period, if deriving year length from star mass)
- Synodic vs. sidereal day relation (rotation rate Ω from solar day and year length)
- Solar declination formula (from obliquity and solar longitude)
- Hour angle and solar zenith angle formula
- Inverse-square law for stellar flux
- Daily-mean insolation formula with sunset hour angle (Berger 1978)

### 6.2 Atmospheric dynamics
- Primitive equations on the sphere:
  - horizontal momentum equations, in vector-invariant form
  - hydrostatic equation
  - continuity (mass conservation) equation
  - thermodynamic energy equation (first law of thermodynamics)
  - moisture conservation equation
- Ideal gas law
- Poisson's equation for potential temperature
- Coriolis parameter (f = 2Ω sin φ) and the β-plane approximation, for analysis
- Sigma (terrain-following) vertical coordinate
- Phillips two-level model (the classical 2-layer formulation), and its
  generalization to any number of sigma layers (Arakawa & Suarez 1983)
- Geostrophic balance, thermal wind equation, Rossby number, Rossby deformation
  radius (used as diagnostics and checks)
- Dry and moist adiabatic lapse rates

### 6.3 Numerics on the sphere
- Spherical centroidal Voronoi tessellation (SCVT) via Lloyd's algorithm (grid quality)
- TRiSK scheme (Thuburn–Ringler–Skamarock–Klemp), finite-volume C-grid operators on Voronoi meshes, as in MPAS
- Discrete gradient, divergence and curl operators; Gauss divergence theorem;
  Stokes' theorem
- Leapfrog time stepping with the Robert–Asselin–Williams (RAW) filter
- Semi-implicit gravity-wave treatment, which leads to a Helmholtz equation
- Preconditioned conjugate gradient method (Helmholtz solver)
- Courant–Friedrichs–Lewy (CFL) condition
- Hyperdiffusion (∇⁴ / del-4) and Rayleigh friction/damping (sponge layer)
- Positive-definite / flux-limited tracer advection (flux-corrected transport, or
  monotonic upwind schemes)
- Williamson et al. (1992) shallow-water test cases, and the Held–Suarez (1994)
  dynamical-core benchmark

### 6.4 Radiation
- Stefan–Boltzmann law
- Beer–Lambert law
- Grey-atmosphere two-stream approximation (Schwarzschild equations)
- Logarithmic CO₂ radiative forcing (Myhre et al. 1998), with pressure broadening
  for optical depth
- Planetary albedo and the energy balance model (equilibrium temperature check)

### 6.5 Moisture, clouds and precipitation
- Clausius–Clapeyron equation (Magnus–Tetens approximation for saturation vapour pressure)
- Relative and specific humidity relations
- Large-scale (stratiform) condensation
- Moist convective adjustment (Manabe 1965), or Betts–Miller relaxation
- Latent heat release
- Orographic lift (vertical velocity from wind blowing across terrain slope)

### 6.6 Surface and boundary layer
- Bulk aerodynamic formulas (momentum drag, sensible heat, latent heat/evaporation)
- Surface energy balance equation
- Manabe bucket model (soil moisture)
- Snow-albedo and ice-albedo feedbacks
- Sub-grid orographic drag (simplified Lott–Miller)

### 6.7 Ocean and sea ice
- Slab mixed-layer heat budget equation
- Wind-stress bulk formula
- Ekman transport
- Sverdrup balance
- Stommel / Munk western boundary current models
- Advection–diffusion equation (ocean heat transport)
- Freezing point of seawater (as a function of salinity)
- Semtner zero-layer sea-ice model, and Stefan's law of ice growth

### 6.8 Hydrology and erosion
- Water balance (P − E), runoff
- Steepest-descent flow direction (D8 analogue on a hex mesh) and flow accumulation
- Priority-flood depression filling (Barnes et al. 2014)
- Lake water balance (inflow vs. evaporation), which decides endorheic or salt lakes
- Stream power law (river incision, for the erosion tool)
- Particle/droplet hydraulic erosion; thermal erosion using the angle of repose (talus)

### 6.9 Climate statistics and classification
- Welford's online algorithm (running mean and variance)
- Köppen–Geiger criteria (Peel et al. 2007 version)
- Whittaker biome diagram
- Potential evapotranspiration (Thornthwaite or Priestley–Taylor), for aridity
- Climate-mode weather generator: AR(1) autoregressive process plus a Markov
  chain for wet/dry days (Richardson-type stochastic weather generator)

### 6.10 World generation and geometry
- Spherical Voronoi diagram (plates)
- Rigid rotation of plates about an Euler pole (v = ω × r)
- Airy isostasy (continental vs. oceanic elevation)
- Fractal noise: simplex noise, fractional Brownian motion, domain warping
- Rodrigues' rotation formula; spherical linear interpolation (slerp);
  barycentric interpolation
- Map projections: equirectangular, later Mollweide and orthographic
- Marching triangles (contour extraction on the mesh)
- Pairwise or Kahan summation (deterministic, accurate reductions)

## 7. Milestones

Milestones and their steps are tracked in [ROADMAP.md](ROADMAP.md).

## 8. Branching and release

- `dev` is the working branch, where all development happens. `main` is the release branch.
- CI (fmt, clippy, docs, and golden-hash tests on macOS, Windows and Linux) runs
  on every push to `dev`; the full test suite runs locally before a push.
- `main` is only ever fast-forwarded from `dev` when a release is shipped.
- CD runs on pushes to `main` and builds a GitHub Release, from M5 onwards
  ([roadmap](ROADMAP.md), M5 step 10).

## 9. Risks

| Risk | Mitigation |
|---|---|
| GCM numerical instability | Conservative CFL margin, diffusion/filtering, clamped planet parameters, idealized tests |
| High-speed tier looks unconvincing | Learn the climatology from enough GCM years; validate the climate mode's statistics against full-GCM runs |
| Performance at fine resolution | Space-filling-curve cell ordering, SoA, rayon, profiling; GPU port planned |
| Determinism vs. parallelism | Fixed-order reductions; avoid float atomics; test in CI |
| No visible app until milestone 5 | The CLI's PNG dumps and animated GIFs give visual feedback from milestone 1 |
