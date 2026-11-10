# Roadmap

What we've done and what we're doing next. How things work is in
[DESIGN.md](DESIGN.md), and open issues (O-n) are in
[OBSERVATIONS.md](OBSERVATIONS.md). Every step ends with an automated test or
a visual check.

**Next: M2.11.8.**

## M1 — Foundations
- [x] 1. Cargo workspace, rustfmt, clippy and CI
- [x] 2. Icosahedral SCVT mesh
- [x] 3. Mesh geometry, TRiSK weights, CSR neighbours, Hilbert-curve ordering
- [x] 4. Two-level grid and point location
- [x] 5. Field registry
- [x] 6. Scheduler, sim clock and calendar
- [x] 7. Planet parameters, orbit and insolation
- [x] 8. CLI: headless runs, field PNGs, insolation GIF
- [x] 9. Insolation tests against analytic values

## M1.5 — Tooling and visibility
- [x] 1. Portable math (`fms-math`) and golden-hash tests
- [x] 2. Wind visualization (Perot reconstruction, arrows)
- [x] 3. Per-system timing
- [x] 4. Criterion benchmarks
- [x] 5. API docs in CI
- [x] 6. Astronomical seasons

## M2 — Headless GCM
- [x] 0. Point probe (CSV and line chart at one location)
- [x] 1. Earth heightmap import with sub-grid statistics
- [x] 2. TRiSK discrete operators
- [x] 3. Explicit shallow-water model (Williamson cases 2, 5, 6)
- [x] 4. Semi-implicit leapfrog, RAW filter, Helmholtz solver
- [x] 5. Two-layer dry dynamical core
- [x] 6. Hyperdiffusion, energy and mass diagnostics, blow-up detection
- [x] 7. Held–Suarez benchmark
- [x] 8. Semi-grey radiation on a temporary slab aquaplanet
- [x] 9. Moisture: humidity transport, convection, condensation, clouds
- [x] 9.5. N-layer vertical structure
- [x] 10. Boundary layer and land surface
  - [x] 10.1. Boundary layer: four layers, 10 m wind, surface drag,
    entrainment, dissipative heating, layering-independent shortwave
  - [x] 10.2. Land surface
    - [x] 10.2.1. Terrain in the climate run; idealized continent
    - [x] 10.2.2. Land and ocean tiles
    - [x] 10.2.3. Land energy balance: skin over a soil store
    - [x] 10.2.4. Soil water bucket and runoff
    - [x] 10.2.5. Snow on land
    - [x] 10.2.6. Sub-grid orographic rain and drag
    - [x] 10.2.7. Single-column BOMEX and DYCOMS cases
    - [x] 10.2.8. Earth run and close-out
- [ ] 11. Two-layer slab ocean and sea ice
  - [x] 11.1. Two-layer ocean in `fms-ocean`, no ice
  - [x] 11.1.5. Early performance pass (2.6× at 2.5k cells)
  - [x] 11.1.6. Test tiers and a shorter ocean spin-up
  - [x] 11.2. Two-band surface albedo
  - [x] 11.3. Sea ice as a third tile
  - [x] 11.4. Snow on sea ice
  - [x] 11.5. Aquaplanet ice-albedo runs
  - [x] 11.6. Ice-sheet surface on land
  - [x] 11.7. Antarctic experiments; freeze-dry low cloud adopted (O-33)
  - [ ] 11.8. **Earth run and close-out**: the accelerated run at 2.5k cells,
    monthly ice extent per hemisphere, budgets, O-30's coastal cells again, and
    the ocean-jump setting (O-37).
- [ ] 12. **Ocean heat transport by currents** (PLAN §4.3):
  - Ekman transport blending into Sverdrup near the equator, with the deep
    layer's return flow and the upwelling it implies
  - Gent–McWilliams eddy transport
  - Wind-driven barotropic gyres carrying mixed-layer heat, solved per basin
    (needs the ocean mask of O-11)
  - Rotation-aware coefficients for other planets
  - Tests: the Bhatnagar et al. 2026 aquaplanet (heat transport ~2–3 PW by
    component, equatorial cold tongue, ice edge, against transport off); the
    ocean's heat budget closes; western intensification in a closed basin
- [ ] 13. **Multi-rate scheduling**: `dt` from the CFL limit and a design
  maximum wind (O-13), rounded to divide the solar day; radiation every ~3 h and
  ocean every 1–6 h (PLAN §4.6).
- [ ] 14. **Idealized planets**: aquaplanet, no tilt, fast and slow rotation,
  with assertions on zonal-mean wind and temperature; grid-imprinting check in
  long-run Held–Suarez, with a fix if visible (O-1).
- [ ] 15. **Earth benchmark, part 1**, against a reference climatology, with a
  tuning pass: trade winds and jets, the tropical rain belt and its migration,
  subtropical highs, monsoons, deserts.
- [ ] 16. **Performance pass**: target speed 2 (~340 steps/s) at the default
  resolution. Solver, operator fusion, scheduler overhead (O-8, O-12, O-13,
  O-29).

## M3 — Hydrology, climate and vegetation
- [ ] 0. `.fmsmap` v0 ([FMSMAP.md](FMSMAP.md)), a heightmap importer, and
  `fms-climate` runs from a file
- [ ] 1. Downscaling onto the fine grid with a lapse-rate correction
- [ ] 2. Rolling monthly temperature and precipitation statistics
- [ ] 3. Köppen–Geiger classification, with agreement against Peel et al. 2007
- [ ] 4. Evolving vegetation, feeding back on albedo, roughness and soil water
- [ ] 5. Flow routing, depression filling, river extraction
- [ ] 6. Lakes and basins; lakes become water surfaces in the GCM
- [ ] 7. Earth benchmark, part 2: Köppen threshold; rivers and basins plausible
- [ ] 8. Export climate results into `.fmsmap`

## M4 — Determinism, persistence and long runs
- [ ] 1. Determinism audit: fixed-order reductions, per-system RNG streams,
  bit-identical replay across thread counts (O-41)
- [ ] 2. Versioned resume saves; mid-run save/load continues bit-identically
- [ ] 3. Climatology learning during full-GCM runs
- [ ] 4. Long runs: a few GCM years alternating with a century or more of
  vegetation, to ~10,000 years (PLAN §4.7)
- [ ] 5. Accelerated spin-up on a background thread, with progress and cancel
- [ ] 6. Validation: unstable drawn deserts fade, stable ones persist; long runs
  match GCM-stepped runs over shorter spans

## M5 — Viewer
View only; the viewer code is shared with `fms-geo`.
- [ ] 1. App shell (eframe + wgpu), sim on its own thread
- [ ] 2. Globe and equirectangular renderers
- [ ] 3. Map modes, colormaps and legends for every field
- [ ] 4. Contours
- [ ] 5. Time controls and HUD
- [ ] 6. Cell inspector with annual graphs
- [ ] 7. Live overlays: wind, currents, clouds, precipitation, snow, sea ice
- [ ] 8. Read-only planet panel
- [ ] 9. Open `.fmsmap` files and heightmaps, resume runs
- [ ] 10. CD: pushes to `main` build a GitHub Release for all three platforms
- [ ] 11. About a dozen preset worlds

## Later
- [ ] Port the hot loops to wgpu compute shaders
- [ ] Tidal locking, moons, extreme obliquity
- [ ] Orbital drift and ice sheets, for runs towards 100,000 years
