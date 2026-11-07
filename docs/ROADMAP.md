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
- [ ] 6. **Astronomical seasons**: season and equinox/solstice dates from the
   orbit (solar longitude) rather than calendar quarters, so eccentric orbits get
   seasons of unequal length. The calendar quarter stays available as a
   convention.

## M2 — Headless GCM
- [ ] 1. **Earth heightmap import** (ETOPO-style equirectangular): resample to the fine
   grid, then aggregate to the coarse grid, keeping sub-grid terrain statistics.
- [ ] 2. **Discrete operators**: TRiSK gradient, divergence, curl, kinetic energy and
   tangential-velocity reconstruction. Test: convergence against analytic
   spherical-harmonic fields; mass and energy conservation.
- [ ] 3. **Shallow-water model (explicit)**: single layer. Pass Williamson test cases
   2 (steady geostrophic flow), 5 (flow over a mountain) and 6 (Rossby–Haurwitz wave).
- [ ] 4. **Semi-implicit time stepping**: leapfrog with the RAW filter; implicit
   gravity-wave terms; Helmholtz solver (preconditioned conjugate gradient). Rerun
   the Williamson tests at the larger `dt`.
- [ ] 5. **2-layer dry dynamical core**: sigma coordinates, hydrostatic balance,
   potential temperature per layer, a surface-pressure equation.
- [ ] 6. **Held–Suarez benchmark**: Newtonian temperature relaxation plus Rayleigh
   friction. Mid-latitude jets and baroclinic eddies should appear.
- [ ] 7. **Dissipation and safety**: hyperdiffusion, sponge layer, global energy and
   mass diagnostics, automatic blow-up detection with a diagnostic dump.
- [ ] 8. **Radiation**: grey two-stream scheme, shortwave absorption, surface albedo,
   CO₂/pressure-dependent optical depth. Replaces the Held–Suarez forcing.
   Test: global energy balance closes.
- [ ] 9. **Moisture**: positive-definite humidity advection, Clausius–Clapeyron
   saturation, large-scale condensation, moist convective adjustment, latent
   heating, precipitation output.
- [ ] 10. **Land surface**: surface energy balance, bulk fluxes, boundary-layer friction,
    bucket soil moisture, snow cover and albedo, orographic precipitation using the
    sub-grid terrain statistics.
- [ ] 11. **Slab ocean and sea ice**: mixed-layer heat budget coupled to the
    atmosphere; Semtner zero-layer sea ice with ice-albedo feedback.
- [ ] 12. **Ocean currents**: currents driven by wind stress through Ekman/Sverdrup
    balance, with western boundary intensification and heat advection.
- [ ] 13. **Multi-rate scheduling**:
    - derive `dt` from the CFL limit, rounded so it divides the solar day
    - radiation every ~3 h, ocean every 1–6 h (§4.6)
- [ ] 14. **Idealized planets**: aquaplanet (symmetric Hadley cell), no-tilt (no
    seasons), fast and slow rotation (number and width of circulation cells).
    Automated assertions on zonal-mean wind and temperature.
- [ ] 15. **Earth benchmark, part 1**: compare against a reference climatology:
    - trade winds and jets
    - the tropical rain belt and its seasonal migration
    - subtropical highs
    - monsoons
    - deserts in the right places

    Includes a parameter-tuning pass.
- [ ] 16. **Performance pass**: profile, parallelize the hot kernels with rayon, check
    the memory layout. Target: speed 2 (~340 steps/s) at the default resolution.

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
