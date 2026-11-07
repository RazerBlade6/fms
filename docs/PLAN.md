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
| Atmosphere | 2-layer primitive-equation model | Smallest model that produces Hadley and Ferrel cells, jets and storms on its own |
| Ocean | Slab mixed layer + wind-driven surface currents + thermodynamic sea ice | Ocean heat transport (western boundary currents) matters for coastal climates |
| Compute | CPU with rayon first, GPU (wgpu compute) later | Easier to debug the physics; data is laid out as struct-of-arrays so it can be ported |
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
| Undo | Command-based; each stroke stores a sparse diff of the cells it changed; terrain only | Cheap; the climate re-settles after undo instead of rewinding |
| Validation | Earth benchmark + idealized test planets, run as regression tests | The only honest way to tune a simplified GCM |
| First milestone | Headless GCM first | Prove the physics before building UI on top of it |

## 3. Architecture

### 3.1 Core concept: field registry + system scheduler

- A **field** is a named, typed array with one value per cell, stored on a specific
  grid level (`terrain` or `climate`). Examples: `elevation`, `surface_temp`,
  `precip`, `wind_u_lower`, `koppen`, and later `soil_fertility` and `population`.
- A **system** declares:
  - which fields it reads and which it writes
  - its tick rate (every step, daily, monthly, on edit)
  - which speed tiers it runs in
- The **scheduler** orders systems by their read/write dependencies and runs
  independent systems in parallel. It also enforces determinism: fixed order and
  ordered reductions.
- **Map modes** are renderers over fields: a colormap, contours, vectors or
  particles. A new field gets a basic map mode for free.
- **Civ layer (future)**: an ECS (`hecs` or `bevy_ecs`) for agents such as
  settlements and populations. Agents read and write fields through the same registry.
  Cells themselves are **not** entities; with over a million cells, an ECS would hurt cache performance.

### 3.2 Grids

- Both levels are subdivided icosahedra. The climate grid is a coarser subdivision
  of the same icosahedron, so the mapping between the two levels is simple
  and precomputed.
- **Fine → coarse**: averaged by area. Sub-grid terrain statistics (mean,
  variance, slope) are kept so the climate grid still feels mountains it can't resolve.
- **Coarse → fine**: interpolated, then corrected for altitude using the lapse
  rate for display and for hydrology.
- Cells are stored as struct-of-arrays, with neighbours in compressed sparse-row
  (CSR) adjacency lists. Cells are ordered along a space-filling curve so neighbours
  sit close together in memory.

### 3.3 Workspace layout

```
crates/
  geo/        icosahedral grids, two-level mapping, projections, picking
  fields/     field registry, SoA storage, serialization
  sched/      system scheduler (tick rates, R/W sets, determinism)
  planet/     planet parameters, orbit, insolation
  worldgen/   Generator trait, plate generator, noise, heightmap import
  atmos/      2-layer GCM (dynamics + physics parameterizations)
  ocean/      slab mixed layer, wind-driven currents, sea ice
  hydro/      runoff, flow routing, depression filling, lakes, basins
  climate/    rolling monthly stats, Köppen, biomes, high-speed climate mode
  sim/        orchestration, speed tiers, spin-up, save/load
  cli/        headless runner, PNG/field dumps, Earth benchmark
  app/        eframe/egui + wgpu renderer, brushes, map modes, inspector
```

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

### 4.2 Atmosphere (2-layer GCM)
- Primitive equations solved with a finite-volume method on the icosahedral grid.
  Upper and lower layers each have wind, temperature and moisture (moisture mostly
  in the lower layer).
- Physics parameterizations:
  - shortwave and longwave radiation (grey, two-band)
  - surface fluxes (bulk formulas)
  - large-scale condensation and simple convective adjustment
  - orographic lift using the sub-grid terrain statistics
  - boundary-layer friction
- Timestep is limited by the CFL condition; roughly 10–50 steps per sim-day,
  depending on resolution and planet radius.

### 4.3 Ocean
- Mixed-layer heat budget, with a mixed-layer depth that can vary by region.
- Wind-driven surface currents from wind stress and Coriolis (Ekman/Sverdrup-style
  balance), with heat carried along by the currents.
- Thermodynamic sea ice: grows and melts, and its albedo feeds back into climate.

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

### 4.6 Speed tiers
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
- Determinism contract: same save + same build = bit-identical sim. Parallel
  reductions use a fixed order.

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
- These tests run in CI as regression checks.

## 6. Milestones

1. **Foundations**: Cargo workspace; `geo` (icosahedral grids, two-level mapping,
   projections), `fields`, `sched`, `planet`; a CLI that renders any field to PNG.
   Test: insolation matches the analytic results.
2. **Headless GCM**: `atmos` + `ocean` running on an Earth heightmap. Tests:
   the idealized planets and the circulation checks of the Earth benchmark.
3. **Hydrology and climate**: `hydro`, Köppen, biomes, and biome feedback into the GCM.
   Test: the Earth Köppen comparison.
4. **Speed tiers and determinism**: climate mode, spin-up, save/load, replay tests.
5. **GUI**: globe and flat renderer, map modes, play controls with speeds 1–5,
   inspector, wind particles.
6. **Editor**: world generation (plate layout, blank start, import), brushes in
   both views, undo/redo, live-edit spin-up.
7. **Later**:
   - port the hot loops to wgpu compute
   - time-stepped plate sim
   - civ layer (ECS agents, soil fertility, aridity, resources from plate and
     crust data)
   - tidal locking and moons

## 7. Risks

| Risk | Mitigation |
|---|---|
| GCM numerical instability | Conservative CFL margin, diffusion/filtering, clamped planet parameters, idealized tests |
| High-speed tier looks unconvincing | Learn the climatology from enough GCM years; validate the climate mode's statistics against full-GCM runs |
| Performance at fine resolution | Space-filling-curve cell ordering, SoA, rayon, profiling; GPU port planned |
| Determinism vs. parallelism | Fixed-order reductions; avoid float atomics; test in CI |
| No visible app until milestone 5 | The CLI's PNG dumps and animated GIFs give visual feedback from milestone 1 |
