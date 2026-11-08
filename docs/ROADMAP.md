# Roadmap

Milestones and steps for the Fantasy Map Simulator. Design details live in
[PLAN.md](PLAN.md); step references like "M5 step 10" point here.

Every step ends with an automated test or a visual check (a PNG/GIF from the CLI).

## M1 — Foundations
- [x] 1. **Workspace and CI**: Cargo workspace with all crate skeletons; rustfmt and
   clippy; GitHub Actions running tests.
- [x] 2. **Icosahedral mesh**: subdivide the icosahedron, build the Voronoi dual (cell
   centres, vertices, edges), and optimize it into an SCVT using Lloyd's algorithm.
   Tests:
   - exactly 12 pentagons
   - Euler characteristic = 2
   - cell areas sum to 4πr²
- [x] 3. **Mesh geometry and layout**: cell areas, edge lengths, dual lengths, edge
   normals and tangents, and the TRiSK weights. Neighbours stored as CSR, cells
   reordered along a Hilbert curve.
- [x] 4. **Two-level grid**: map each fine cell to its coarse parent with area weights.
   Point → cell lookup (for picking and rasterization).
- [x] 5. **Field registry**: typed struct-of-arrays fields tagged by grid level;
   double buffering; field metadata (units, valid range, default colormap).
- [x] 6. **Scheduler**:
   - `System` trait with read/write sets and a rate as a multiple of `dt`
   - deterministic ordering, with independent systems run in parallel
   - sim clock and calendar (solar day, year, seasons)
- [x] 7. **Planet**:
   - parameter struct with validation and clamping
   - derived quantities: Ω, f, sidereal day
   - orbit (Kepler) and insolation (instantaneous and daily mean)
- [x] 8. **CLI**: run a world headless; dump any field to an equirectangular PNG with
   a colormap; render a year of insolation as a GIF.
- [x] 9. **Tests**:
   - insolation matches analytic values (equinox and solstice, polar night)
   - annual-mean symmetry at zero eccentricity

## M1.5 — Tooling and visibility
Added after the M1 review: tooling the GCM work depends on.
- [x] 1. **Portable math**: `fms-math` wraps `libm` for all transcendental
   functions; `clippy.toml` forbids the std versions; golden-hash tests prove
   results are bit-identical on macOS, Windows and Linux.
- [x] 2. **Wind visualization**: reconstruct a cell-centred vector field from edge
   normal velocities (Perot-style reconstruction); CLI renders vector fields as
   speed plus arrows or streamlines, and edge/vertex fields via cell averages.
   Needed to debug every dynamics step in M2.
- [x] 3. **Per-system timing**: wall-clock time for each system in the
   `StepReport`; `fms` prints a timing breakdown per run.
- [x] 4. **Benchmarks**: `criterion` suite for mesh build, geometry, two-level
   transfers, point location and scheduler overhead. Run locally (CI timing is
   too noisy) to catch performance regressions.
- [x] 5. **API docs in CI**: `cargo doc` with `-D warnings` (broken intra-doc
   links fail the build).
- [x] 6. **Astronomical seasons**: season and equinox/solstice dates from the
   orbit (solar longitude) rather than calendar quarters, so eccentric orbits get
   seasons of unequal length. The calendar quarter stays available as a
   convention.

## M2 — Headless GCM
- [x] 0. **Point probe**: `fms probe --lat --lon -f FIELD [-f …] --days D` writes a
   time series at one location to `output/` as CSV plus a simple line chart. Needed
   to check diurnal and seasonal cycles and energy drift in the steps below.
- [x] 1. **Earth heightmap import** (ETOPO-style equirectangular): resample to the fine
   grid, then aggregate to the coarse grid, keeping sub-grid terrain statistics.
- [x] 2. **Discrete operators**: TRiSK gradient, divergence, curl, kinetic energy and
   tangential-velocity reconstruction. Test: convergence against analytic
   spherical-harmonic fields; mass and energy conservation.
- [x] 3. **Shallow-water model (explicit)**: single layer. Pass Williamson test cases
   2 (steady geostrophic flow), 5 (flow over a mountain) and 6 (Rossby–Haurwitz wave).
- [x] 4. **Semi-implicit time stepping**: leapfrog with the RAW filter; implicit
   gravity-wave terms; Helmholtz solver (preconditioned conjugate gradient). Rerun
   the Williamson tests at the larger `dt`.
- [x] 5. **2-layer dry dynamical core**: sigma coordinates, hydrostatic balance,
   potential temperature per layer, a surface-pressure equation. Equal-mass
   layers (σ = 0.5); explicit RK4 first, then semi-implicit with a vertical-mode
   split. Tests: an atmosphere at rest over mountains stays at rest; mass and
   energy conserved; baroclinic instability grows at the rate linear theory
   predicts.
- [x] 6. **Dissipation and safety**: hyperdiffusion, global energy and mass
   diagnostics, automatic blow-up detection with a diagnostic dump. Moved ahead
   of Held–Suarez at the end of M2.5: without dissipation the core collapses to
   the grid scale ~8 days into a storm life cycle
   ([OBSERVATIONS.md](OBSERVATIONS.md) O-14). Done as ∇⁴ (order configurable),
   at the strongest strength the step allows (O-17). The sponge layer is
   deferred until there are more levels (O-16). Tests: a storm life cycle
   survives 30 days with mass and θ exact and the energy budget closed; rest
   over a mountain is untouched; blow-ups are caught before they become NaNs.
- [x] 7. **Held–Suarez benchmark**: Newtonian temperature relaxation plus Rayleigh
   friction. Mid-latitude jets and baroclinic eddies should appear. Use ∇⁴, not
   the paper's ∇⁸ ([OBSERVATIONS.md](OBSERVATIONS.md) O-17). Done with the
   forcing averaged over each layer's mass (`fms held-suarez`). 1200 days at
   10k cells: upper jets of 26 m/s at 49°, north and south within 0.5 m/s, and an
   energy budget that closes (Newtonian heating 0.38 W/m² balanced by friction
   and diffusion). Tests (k = 16, 450 days): jets, hemispheric agreement,
   steady and closed energy budget. The lower layer is not the surface (O-18).
- [x] 8. **Radiation**: grey two-stream scheme, shortwave absorption, surface albedo,
   CO₂/pressure-dependent optical depth. Replaces the Held–Suarez forcing.
   Test: global energy balance closes. Done as semi-grey two-stream radiation
   (CO₂ band-wing model, window, k-distribution water), with water vapour at
   fixed relative humidity and a temporary 10 m slab aquaplanet with bulk
   sensible and latent heat fluxes and convective adjustment (`fms climate`,
   PLAN §4.2). Single column: Earth-like climate (287 K, OLR 247 W/m²),
   3.6 W/m² and 2.4 K per CO₂ doubling, sane across 0–100,000 ppm, 0.1–10 bar,
   sunlight and gravity. 3-D (5 years, 2.5k cells): global mean ~287 K,
   planetary albedo ~0.29, equator 300 K and poles ~245 K, and the planet's
   energy budget closes to ~0.01 W/m².
   Calibrated rather than resolved: O-19; no cloud or ice feedbacks yet: O-20.
- [x] 9. **Moisture**: positive-definite humidity advection, Clausius–Clapeyron
   saturation, large-scale condensation, moist convective adjustment, latent
   heating, precipitation output. Replaces step 8's fixed-relative-humidity
   water vapour and immediate condensation; diagnose clouds to replace the
   fixed cloud deck ([OBSERVATIONS.md](OBSERVATIONS.md) O-20). Done as
   humidity per layer carried by flux-corrected transport (water exact to
   round-off), Betts–Miller convection towards the moist adiabat,
   large-scale condensation with re-evaporation, rain or snow, and two
   diagnosed cloud sheets acting on sunlight and longwave, plus Rayleigh
   scattering (PLAN §4.2). Aquaplanet (5 years, 2.5k cells): ~297 K
   (ice-free, O-23), albedo ~0.29, cloud ~0.63, layer humidities as
   calibrated, rain belt at the equator; water budget exact, energy budget
   to ~0.1 W/m² (time stepping).
   Known gaps: the subtropics rain out their own evaporation (O-22), the
   water cycle is too strong (O-24), clouds are calibrated (O-25).
- [x] 9.5. **N-layer vertical structure**: generalize the two-layer model to
    any number of σ layers with arbitrary interfaces: the Arakawa–Suarez
    (1983) energy-conserving vertical scheme, one vertical mode per layer in
    the semi-implicit solver, flux-corrected humidity transport across every
    interface, and radiation and moist physics per layer. The climate stays
    on two layers; step 10 moves it to three, with a thin bottom layer for
    the boundary layer, and recalibrates once. Test: two layers reproduce
    step 9 to round-off (Held–Suarez, baroclinic wave, energy and water
    budgets); with three to five layers, mass, Θ, water and energy are
    conserved and Held–Suarez and the moist aquaplanet run stably.
    Done (`levels.rs`, `LayeredModel`, `fms … --levels`): with step 9's
    semi-implicit reference, two layers reproduce step 9 bit for bit; with
    three and five, energy is conserved exactly in space, mass, Θ and water
    to round-off, and Held–Suarez (200 days) and the moist aquaplanet (60
    days) run stably at 7.5× the explicit step. Two changes beyond the
    refactor, both needed for more layers: the semi-implicit reference is
    now an isothermal 300 K atmosphere (Simmons et al. 1978; the start's own
    profile blew up three-layer Held–Suarez at day 28), and shallow
    convection only moves water from layers above their reference humidity
    to those below (it drew dry layers negative). The two-layer climate is
    unchanged (5 years: 297.0 K, albedo 0.286, cloud 0.63, 4.83 mm/day;
    energy closes to 0.03 W/m²). Found: shortwave water absorption depends
    on the layering (O-27); with more than one interface the vertical modes
    are non-normal (O-28).
- [ ] 10. **Boundary layer and land surface**, in two reviewed parts:
    - [ ] 10.1. **Boundary layer** (on the aquaplanet):
      - make the shortwave independent of the layering (O-27): the
        near-infrared band carried separately, each layer's absorption from
        the cumulative water path;
      - move the climate to three layers with a thin (~1.4 km) bottom one;
      - diagnose a 10 m surface wind from the bottom layer (log profile,
        roughness, stability), replacing the lower-layer wind
        ([OBSERVATIONS.md](OBSERVATIONS.md) O-18);
      - bulk surface drag and stability-dependent heat and moisture
        exchange in place of the Rayleigh friction;
      - return friction's and diffusion's dissipated energy as heat (O-21);
      - stability-dependent mixing out of the bottom layer that keeps
        evaporated water near the surface (O-22);
      - recalibrate once on three layers.

      Tests: surface sunlight converges with more layers; the subtropics
      rain ≲2 mm/day; Held–Suarez-like surface westerlies under the jets;
      the top of the atmosphere balances in equilibrium; energy and water
      budgets close.
    - [ ] 10.2. **Land surface**:
      - terrain into the climate run: surface geopotential, land fraction
        and sub-grid height spread per cell (check spurious σ-coordinate
        winds over steep terrain);
      - land and ocean tiles per cell, each with its own surface
        temperature and fluxes, weighted by land fraction (the ocean tile
        stays the 10 m slab until step 11);
      - land energy balance: a skin temperature over a thin soil heat
        store, solved implicitly;
      - Manabe bucket soil moisture with runoff;
      - snow mass, melt and snow albedo;
      - orographic rain and sub-grid orographic drag from the height
        spread.

      Developed on an idealized continent, then run on the Earth
      heightmap. Tests: energy and water budgets close (soil water and
      snow included); land swings more than ocean daily and seasonally.
- [ ] 11. **Two-layer slab ocean and sea ice** (the dynamical slab ocean of
    Codron 2012 and Bhatnagar et al. 2026, PLAN §4.3): a 50 m mixed layer
    over a 150 m deep layer per ocean cell, horizontal diffusion in both,
    and convective adjustment when the surface layer is colder than the deep
    one. Sea ice holds the mixed layer at freezing and grows or melts with
    the excess energy; a prognostic ice surface temperature with a linear
    profile, snow on ice with its own heat capacity, and an albedo that rises
    with thickness from the ocean's, separately in visible and near-infrared
    light. Ice-albedo feedback (O-20). Tests: energy conserved through
    freezing and melting; the ice edge and seasonal cycle are plausible.
- [ ] 12. **Ocean heat transport by currents**:
    - Ekman transport from the surface wind stress, `M = (ετ + f k×τ)/(ε² + f²)`
      at the edges, with the deep layer's return flow and the upwelling or
      downwelling their divergence implies; blended near the equator into
      the Sverdrup transport of the wind-stress curl. Needs step 10's
      surface wind (O-18).
    - Gent–McWilliams eddy transport, proportional to the isotherm slope
      between the layers.
    - A wind-driven barotropic gyre: a Stommel–Munk streamfunction from the
      wind-stress curl, β and friction, solved over each basin with the
      conjugate-gradient solver, whose currents carry mixed-layer heat:
      western boundary currents for coastal climates, and the current map.
    - Rotation-aware coefficients: cap the Sverdrup term for slow rotators
      and scale the eddy coefficients with the deformation radius (the
      paper's are calibrated for Earth).
    - Tests: the aquaplanet benchmark of Bhatnagar et al. 2026 (heat
      transport ~2–3 PW by component, an equatorial cold tongue, the ice
      edge, against a run with transport off); the ocean's heat budget
      closes; western intensification in a closed basin.
- [ ] 13. **Multi-rate scheduling**:
    - derive `dt` from the CFL limit, rounded so it divides the solar day. Size
      it from a design maximum wind, not the current one: flows accelerate
      ([OBSERVATIONS.md](OBSERVATIONS.md) O-13)
    - radiation every ~3 h, ocean every 1–6 h (§4.6)
- [ ] 14. **Idealized planets**: aquaplanet (symmetric Hadley cell), no-tilt (no
    seasons), fast and slow rotation (number and width of circulation cells).
    Automated assertions on zonal-mean wind and temperature. Also a
    **grid-imprinting check**: measure how strongly the pentagons' pattern shows up
    in the long-run Held–Suarez climate ([OBSERVATIONS.md](OBSERVATIONS.md) O-1), and
    apply a fix if it is visible. Already visible in Williamson case 2 as error
    dipoles at the pentagons (±1.5 m of 3 km after 5 days).
- [ ] 15. **Earth benchmark, part 1**: compare against a reference climatology:
    - trade winds and jets
    - the tropical rain belt and its seasonal migration
    - subtropical highs
    - monsoons
    - deserts in the right places

    Includes a parameter-tuning pass.
- [ ] 16. **Performance pass**: profile, parallelize the hot kernels with rayon, check
    the memory layout. Target: speed 2 (~340 steps/s) at the default resolution.
    Known costs going in ([OBSERVATIONS.md](OBSERVATIONS.md)):
    - per-call parallel overhead dominates the dynamics (O-12): fuse operators
    - the Helmholtz solve takes ~36 CG iterations per step (O-13): better
      preconditioner (multigrid or Chebyshev), looser tolerance
    - scheduler overhead (O-8)

## M3 — Hydrology and climate
- [ ] 1. **Downscaling**: interpolate climate-grid fields onto the fine grid, with a
   lapse-rate correction for elevation.
- [ ] 2. **Rolling statistics**: per-cell monthly temperature and precipitation means
   and variances (Welford), averaged over multiple years.
- [ ] 3. **Köppen–Geiger classification**, plus an agreement metric against the
   Peel et al. 2007 Earth map.
- [ ] 4. **Biomes**: Whittaker classification. Biome sets surface albedo, roughness
   and soil water capacity, which feed back into the GCM.
- [ ] 5. **Flow routing**: steepest-descent directions, priority-flood depression
   filling, flow accumulation, river extraction above a threshold.
- [ ] 6. **Lakes and basins**: lake water balance; endorheic basins and salt lakes;
   lakes become water surfaces in the GCM.
- [ ] 7. **Earth benchmark, part 2**: hit the Köppen agreement threshold; major rivers
   and basins appear in plausible places.

## M4 — Speed tiers, determinism and persistence
- [ ] 1. **Determinism audit**:
   - fixed-order reductions and pairwise summation
   - a separate seeded RNG stream for each system
   - test: bit-identical replay across different thread counts
- [ ] 2. **Save format**: versioned binary format with field tags and migration hooks.
   Round-trip tests; mid-run save/load continues bit-identically.
- [ ] 3. **Climatology learning**: per-cell, per-month mean, variance and
   autocorrelation, collected during full-GCM runs.
- [ ] 4. **Climate mode**: stochastic weather generator (AR(1) process plus wet/dry
   Markov chain) driven by the learned climatology.
- [ ] 5. **Tier switching**: hand off cleanly between the GCM and climate mode in both
   directions, starting the GCM from a state consistent with the climatology.
- [ ] 6. **Spin-up**: accelerated spin-up (asynchronous ocean coupling) on a background
   thread, with progress reporting and cancellation.
- [ ] 7. **Validation**: climate-mode statistics match full-GCM statistics within
   tolerance.

## M5 — GUI
- [ ] 1. **App shell**: eframe + wgpu. The sim runs on its own thread and sends
   double-buffered field snapshots to the UI over a channel.
- [ ] 2. **Renderers**: globe (hex mesh, orbit camera) and flat equirectangular view.
- [ ] 3. **Map mode framework**: colormaps and legends; every field automatically gets
   a map mode; climate-grid fields are shown downscaled.
- [ ] 4. **Contours**: marching triangles on the mesh, with configurable intervals.
- [ ] 5. **Time controls**: play/pause, speeds 1–5, date and season HUD, sim-rate
   readout.
- [ ] 6. **Inspector**: click to pick a cell; show current values and annual
   temperature and precipitation graphs.
- [ ] 7. **Live overlays**: wind particles, ocean current vectors, clouds and
   precipitation, snow and sea-ice extent.
- [ ] 8. **Planet panel**: edit planet parameters, which re-derives `dt` and triggers
   spin-up; shows warnings when values are clamped.
- [ ] 9. **Save/load UI.**
- [ ] 10. **CD pipeline**: a GitHub Actions workflow triggered by pushes to `main`.
    It builds release binaries for macOS, Windows and Linux and publishes them as
    a GitHub Release. Releases are cut on every push to `main` from M5 onwards.

## M6 — Editor
- [ ] 1. **Generator trait**: blank-ocean world; heightmap import UI.
- [ ] 2. **Plate generator**:
   - Voronoi plates, each with an Euler pole and a crust type
   - boundary classification from relative motion
   - elevation synthesis from boundaries plus fractal noise
   - sea level set from the ocean fraction
- [ ] 3. **Brush framework**:
   - picking in both views
   - falloff curves
   - stroke → command holding a sparse diff
   - undo/redo stack
- [ ] 4. **Sculpt brushes**: raise, lower, smooth, flatten, noise.
- [ ] 5. **Erosion tools**: hydraulic (droplet) and thermal (talus).
- [ ] 6. **Water tools**: set sea level, paint river (carves a channel), fill lake.
- [ ] 7. **Live-edit pipeline**: track dirty regions, recompute only the affected parts
   (grid aggregation, sub-grid stats, hydrology), then start spin-up.

## M7 — Later
- [ ] Port the hot loops to wgpu compute shaders
- [ ] Time-stepped plate tectonics (Cortial et al. 2019 style), behind the `Generator` trait
- [ ] Civ layer: ECS agents, soil fertility, aridity indices, resources from plate
  and crust data
- [ ] Tidal locking, moons, extreme obliquity
