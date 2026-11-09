# Roadmap

Milestones and steps for `fms-climate`. Design details live in
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
   `StepReport`; `fms-climate` prints a timing breakdown per run.
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
- [x] 0. **Point probe**: `fms-climate probe --lat --lon -f FIELD [-f …] --days D` writes a
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
   forcing averaged over each layer's mass (`fms-climate held-suarez`). 1200 days at
   10k cells: upper jets of 26 m/s at 49°, north and south within 0.5 m/s, and an
   energy budget that closes (Newtonian heating 0.38 W/m² balanced by friction
   and diffusion). Tests (k = 16, 450 days): jets, hemispheric agreement,
   steady and closed energy budget. The lower layer is not the surface (O-18).
- [x] 8. **Radiation**: grey two-stream scheme, shortwave absorption, surface albedo,
   CO₂/pressure-dependent optical depth. Replaces the Held–Suarez forcing.
   Test: global energy balance closes. Done as semi-grey two-stream radiation
   (CO₂ band-wing model, window, k-distribution water), with water vapour at
   fixed relative humidity and a temporary 10 m slab aquaplanet with bulk
   sensible and latent heat fluxes and convective adjustment (`fms-climate climate`,
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
    Done (`levels.rs`, `LayeredModel`, `fms-climate … --levels`): with step 9's
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
    - [x] 10.1. **Boundary layer** (on the aquaplanet):
      - make the shortwave independent of the layering (O-27): the
        near-infrared band carried separately, each layer's absorption from
        the cumulative water path;
      - move the climate to four layers (σ 0.25/0.5/0.85): a thin (~1.4 km)
        bottom one, and a top split at 250 hPa so that the Hadley outflow
        isn't averaged with the stratosphere;
      - diagnose a 10 m surface wind from the bottom layer (log profile,
        roughness, stability), replacing the lower-layer wind
        ([OBSERVATIONS.md](OBSERVATIONS.md) O-18);
      - bulk surface drag and stability-dependent heat and moisture
        exchange in place of the Rayleigh friction;
      - return friction's and diffusion's dissipated energy as heat (O-21);
      - stability-dependent mixing out of the bottom layer that keeps
        evaporated water near the surface (O-22);
      - recalibrate once on four layers.

      Tests: surface sunlight converges with more layers; the subtropics
      rain ≲2 mm/day; Held–Suarez-like surface westerlies under the jets;
      the top of the atmosphere balances in equilibrium; energy and water
      budgets close.

      Done, with one target unmet. Four layers (σ 0.25/0.5/0.85); surface
      sunlight agrees to 0.03 W/m² across layerings; Held–Suarez has +3 to
      +4 m/s under its jets; five slab years settle at 290.7 K with the top
      of the atmosphere balanced (O-23). At fixed APE CONTROL sea
      temperatures: Hadley cell 80×10⁹ kg/s at 15° (45 on three layers),
      imbalance +14 W/m² (APE ~13), albedo 0.31, rain 3.4 mm/day. Unmet:
      the subtropics rain 2.5–2.9 mm/day, set by the circulation rather
      than the boundary layer (O-22, for step 15).
    - [ ] 10.2. **Land surface**, developed on an idealized continent, then
      run on the Earth heightmap:
      - [x] 10.2.1. **Terrain in the climate run**: surface geopotential
        from the world's terrain; a resting start in hydrostatic balance
        with it (the planet's surface pressure is the sea-level value); an
        idealized continent (`--continent`). Test: a stratified atmosphere
        at rest over the continent stays nearly at rest, measuring the
        spurious σ-coordinate winds; smooth the terrain if they demand it.
        Done, with ten conservative diffusion passes over the terrain the
        atmosphere sees: the persistent wind error on Earth at 10k cells
        falls from 5.7 to 0.87 m/s, 0.48 over the continent, 0.01 flat
        (O-30). 1.2 years of full physics run stably on the continent and
        on Earth's terrain, energy closing to under 0.01 W/m². The surface
        is still the slab until the next step.
      - [x] 10.2.2. **Land and ocean tiles**: land fraction into the
        physics; each cell has a land and an ocean tile with its own surface
        temperature, albedo, roughness, radiation and fluxes, weighted by
        land fraction (land starts as a slab; the ocean tile stays the 10 m
        slab until step 11). Tests: no land reproduces the aquaplanet bit
        for bit; energy and water budgets close with tiles.
        Done: land is a 1 m wet slab, albedo 0.2, roughness 0.1 m. No land
        matches the pre-tile aquaplanet's state hash exactly; with a
        continent, energy closes to −0.02 W/m² over 90 days and water to
        6×10⁻¹⁴. Lifted coastal ocean stays lifted (O-30).
      - [ ] 10.2.3. **Land energy balance**: a skin temperature over a thin
        soil heat store, solved implicitly. Tests: the implicit solve
        against an analytic case; land swings more than ocean daily and
        seasonally.
      - [ ] 10.2.4. **Soil water**: Manabe bucket limiting evaporation, with
        runoff. Test: the water budget closes with soil water and runoff.
      - [ ] 10.2.5. **Snow**: snow mass, melt and snow albedo. Test: energy
        and water budgets close through melting; snow comes and goes with
        the seasons.
      - [ ] 10.2.6. **Sub-grid orography**: orographic rain and drag from
        the height spread. Test: more rain on the windward side of the
        ridge than the lee; the drag slows the flow over rough terrain and
        its energy returns as heat.
      - [ ] 10.2.7. **Single-column cases** for the boundary layer and
        convection (BOMEX trade cumulus, DYCOMS stratocumulus): the column
        physics as a per-column function forced by the cases' prescribed
        large-scale tendencies, checked against their observations and
        large-eddy simulations in seconds per run.
      - [ ] 10.2.8. **Earth run and close-out**: the full model on the
        Earth heightmap, with terrain land, ocean, soil water and snow
        together; energy and water budgets close over a year.
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

## M3 — Hydrology, climate and vegetation
- [ ] 0. **`.fmsmap` v0**: the container, planet settings, elevation, land and
   optional vegetation ([FMSMAP.md](FMSMAP.md)); a heightmap-to-`.fmsmap` importer;
   `fms-climate` runs from a file.
- [ ] 1. **Downscaling**: interpolate climate-grid fields onto the fine grid, with a
   lapse-rate correction for elevation.
- [ ] 2. **Rolling statistics**: per-cell monthly temperature and precipitation means
   and variances (Welford), averaged over multiple years.
- [ ] 3. **Köppen–Geiger classification**, plus an agreement metric against the
   Peel et al. 2007 Earth map.
- [ ] 4. **Evolving vegetation**: a vegetation amount per land tile, starting from
   the imported layer and following its layer state, moving towards what the
   climate supports. It feeds back on albedo, roughness and soil water capacity.
- [ ] 5. **Flow routing**: steepest-descent directions, priority-flood depression
   filling, flow accumulation, river extraction above a threshold.
- [ ] 6. **Lakes and basins**: lake water balance; endorheic basins and salt lakes;
   lakes become water surfaces in the GCM.
- [ ] 7. **Earth benchmark, part 2**: hit the Köppen agreement threshold; major rivers
   and basins appear in plausible places.
- [ ] 8. **Export**: climate results into `.fmsmap` (monthly climate, vegetation, ice,
   rivers), passing through layers the climate doesn't understand.

## M4 — Determinism, persistence and long runs
- [ ] 1. **Determinism audit**:
   - fixed-order reductions and pairwise summation
   - a separate seeded RNG stream for each system
   - test: bit-identical replay across different thread counts
- [ ] 2. **Resume saves**: a climate-only versioned binary format with field tags and
   migration hooks. Round-trip tests; mid-run save/load continues bit-identically.
- [ ] 3. **Climatology learning**: per-cell, per-month mean, variance and
   autocorrelation, collected during full-GCM runs.
- [ ] 4. **Long runs**: alternate a few GCM years with a century or more of
   vegetation, for runs of about 10,000 years (PLAN §4.7).
- [ ] 5. **Spin-up**: accelerated spin-up (asynchronous ocean coupling) on a background
   thread, with progress reporting and cancellation.
- [ ] 6. **Validation**: an unstable drawn desert fades and a stable one persists;
   long-run results match runs with the GCM stepped throughout, over shorter spans.

## M5 — Viewer
View only; the viewer code is shared with `fms-geo`.
- [ ] 1. **App shell**: eframe + wgpu. The sim runs on its own thread and sends
   double-buffered field snapshots to the UI over a channel.
- [ ] 2. **Renderers**: globe (hex mesh, orbit camera) and flat equirectangular view.
- [ ] 3. **Map mode framework**: colormaps and legends; every field automatically gets
   a map mode; climate-grid fields are shown downscaled.
- [ ] 4. **Contours**: marching triangles on the mesh, with configurable intervals.
- [ ] 5. **Time controls**: play/pause, live speeds, starting a long run, date and
   season HUD, sim-rate readout.
- [ ] 6. **Inspector**: click to pick a cell; show current values and annual
   temperature and precipitation graphs.
- [ ] 7. **Live overlays**: wind particles, ocean current vectors, clouds and
   precipitation, snow and sea-ice extent.
- [ ] 8. **Planet panel**: shows the planet's settings, read-only.
- [ ] 9. **Opening worlds**: open `.fmsmap` files and heightmaps, and resume runs.
- [ ] 10. **CD pipeline**: a GitHub Actions workflow triggered by pushes to `main`.
    It builds release binaries for macOS, Windows and Linux and publishes them as
    a GitHub Release. Releases are cut on every push to `main` from M5 onwards.
- [ ] 11. **Presets**: about a dozen shipped worlds, as `.fmsmap` files.

## Later
- [ ] Port the hot loops to wgpu compute shaders
- [ ] Tidal locking, moons, extreme obliquity
- [ ] Orbital drift and ice sheets, for runs towards 100,000 years

## Moved to other modules
The earlier M6 and part of M7, kept here until those projects start
(PLAN §1.1):
- **`fms-editor`**: brush framework (picking in both views, falloff curves, undo
  and redo as sparse diffs); sculpt brushes (raise, lower, smooth, flatten, noise);
  water tools (set sea level, paint river, fill lake); drawing vegetation; a
  blank-ocean generator and heightmap import.
- **`fms-geo`**: a plate generator (Voronoi plates with Euler poles and crust types,
  boundary classification, elevation synthesis, sea level from the ocean fraction);
  time-stepped plate tectonics (Cortial et al. 2019 style); hydraulic and thermal
  erosion.
- **`fms-civ`**: ECS agents, soil fertility, aridity indices, resources from plate
  and crust data, and the biome → typical climate table for maps with no climate
  run.
