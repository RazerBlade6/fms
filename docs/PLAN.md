# Plan

The canonical design of `fms-climate`: what we set out to build and the
decisions that shape it. Update it whenever a decision overrules it. Progress
is in [ROADMAP.md](ROADMAP.md), and how built subsystems work is in
[DESIGN.md](DESIGN.md).

`fms-climate` is a climate-centred simulator of theoretical planets. Given a
planet and its terrain, it simulates the atmosphere, ocean, ice, land surface,
rivers and lakes, and vegetation over time.

## 1. Goals

- **Customizable planet**: radius, gravity, rotation, tilt, year length,
  eccentricity, luminosity, surface pressure, CO₂ and ocean fraction (§4.1).
- **Imported worlds**: `.fmsmap` files ([FMSMAP.md](FMSMAP.md)), global
  heightmaps, and about a dozen presets.
- **Emergent climates**: climate zones come out of the simulation. Frozen or
  uninhabitable worlds are valid results.
- **Evolving vegetation**: a drawn desert the climate doesn't support fades.
- **A viewer, not an editor**: map modes for climate (Köppen), vegetation,
  terrain, wind, temperature, precipitation, currents, sea ice and any
  registered field. Play shows a run live.
- **Time**: speed 1 = one planetary solar day per second. Long runs reach
  about 10,000 years, and at most about 100,000 (§4.7).
- **Exports**: climate results go back into `.fmsmap`.

### 1.1 The `fms` family

Separate modules exchanging `.fmsmap` files. Each uses only the layers it
simulates and passes the rest through. Coupling goes through files
(geo → climate → geo, or civ → climate).

| Module | Role | Interface |
|---|---|---|
| `fms-editor` | Make custom maps at any level of detail | GUI, edits |
| `fms-climate` | This project: climate and hydrosphere | Viewer only |
| `fms-geo` | Plate tectonics and erosion (design later) | Viewer only |
| `fms-civ` | Civilisation simulation and editing (design later) | GUI, edits |

Typical paths: editor → geo → climate → civ, or editor (with biomes) → civ,
which then uses a biome → typical-climate table and warns about it.

All modules stay in this Cargo workspace until the format is stable.
`fms-editor`'s design is in [editor/PLAN.md](editor/PLAN.md).

## 2. Key decisions

| Area | Decision |
|---|---|
| Language | Rust |
| Platform | Desktop, cross-platform (macOS first) |
| UI and rendering | `eframe`/`egui` panels + custom `wgpu` renderer, shared by view-only modules |
| Geometry | Simulate on a sphere; flat maps are projections |
| Sphere grid | Icosahedral geodesic grid (Voronoi hexagons + 12 pentagons) |
| Climate fidelity | Simplified GCM with real fluid dynamics; weather emerges |
| Atmosphere | Primitive equations on any number of σ layers; four by default, the bottom one a thin boundary layer |
| Ocean | Two-layer dynamical slab with Ekman/Sverdrup, eddy and diffusive heat transport, wind-driven gyres, and thermodynamic sea ice |
| Compute | CPU with rayon first, GPU later; struct-of-arrays data |
| Time stepping | Semi-implicit, base step from the CFL limit (~20–30 min at default resolution); slower systems run every N steps |
| Day length | "1 day" = the planet's solar day |
| Long runs | Full GCM live; long runs alternate a few GCM years with a century or more of vegetation |
| Resolution | Fine terrain (~1–2.5M cells) and coarse climate grid (~10–40k cells) |
| Terrain | Fixed during a run; tectonics and erosion belong to `fms-geo` |
| Hydrology | Emergent: runoff routed into rivers, depressions become lakes |
| Classification | Köppen–Geiger as a derived field; vegetation is simulated |
| Planet range | Earth-like range + eccentricity; extremes clamped with warnings |
| Determinism | Seeded RNG, fixed timestep, deterministic reductions, versioned saves |
| Portable math | Transcendentals come from `libm` via `fms-math`, enforced by `clippy.toml` |
| Validation | Earth benchmark + idealized planets, as regression tests |
| Interchange | `.fmsmap` files with planet settings and named layers |
| Layer states | Locked, initial (default) or absent, per layer |
| First milestone | Headless GCM before any UI |

## 3. Architecture

### 3.1 Fields and systems

Simulation state is a registry of named, typed **fields** per grid level.
**Systems** declare what they read and write and how often they run. A
**scheduler** runs them in parallel with results identical to running them in
order. **Map modes** are renderers over fields, so a new field gets one for
free.

### 3.2 Grids

Two subdivided-icosahedron grids: a fine terrain grid and a coarse climate grid,
linked by a precomputed fine → coarse mapping. Sub-grid terrain statistics let
the climate grid feel mountains it can't resolve.

### 3.3 Workspace layout

```
crates/
  math/         portable transcendental functions (libm)
  grid/         icosahedral grids, two-level mapping, projections, picking
  fields/       field registry, SoA storage, serialization
  sched/        system scheduler
  planet/       planet parameters, orbit, insolation
  worldgen/     heightmap import, idealized terrain
  atmos/        layered GCM (dynamics + physics)
  ocean/        slab ocean, wind-driven currents, sea ice
  hydro/        runoff, flow routing, depression filling, lakes, basins
  climatology/  rolling monthly stats, Köppen, vegetation
  sim/          orchestration, spin-up, long runs, save/load
  cli/          the `fms-climate` headless runner
  viewer/       eframe/egui + wgpu renderer, map modes, inspector
```

A crate for `.fmsmap` joins these in M3 step 0. Packages have an `fms-` prefix.
The binaries are `fms-climate` and `fms-viewer`. The sim crates never depend on
the viewer.

## 4. Subsystems

### 4.1 Planet and insolation

Supported ranges; values outside are clamped with a warning, and NaN or infinite
values become Earth's.

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

`S₀ = 1361 W/m² · L / a²`. Outside 500–3,000 W/m² the planet is allowed but
flagged as likely frozen or runaway.

### 4.2 Atmosphere

A layered primitive-equation GCM with semi-grey radiation, prognostic moisture,
diagnosed clouds, a thin boundary layer, and a surface of ocean, sea-ice and
land tiles. Built; see [DESIGN.md](DESIGN.md).

### 4.3 Ocean

A dynamical slab after Codron (2012) and Bhatnagar et al. (2026), plus
wind-driven gyres. Where the paper and its code
(`docs/sources/generic-pcm-r3423/`) differ, we follow the code. Built; see
[DESIGN.md](DESIGN.md). Heat moves as fluxes across cell edges, so the budget
closes:

- **Ekman transport** from the surface wind stress, blending into the Sverdrup
  transport of the wind-stress curl near the equator, with the deep layer's
  return flow and the vertical exchange it implies.
- **Gent–McWilliams** eddy transport along the isotherm slope between layers.
- **Gyres**: a Stommel streamfunction from the wind-stress curl, β and linear
  friction, with each landmass's value found by the island rule, carrying heat
  in both layers. Stommel rather than Stommel–Munk: the grid can't resolve
  either boundary layer, and linear friction needs no extra boundary condition.
- **Rotation**: the Sverdrup term is capped for slow rotators; the eddy
  coefficients scale with the deformation radius.

Regional mixed-layer depths are a later refinement. Salinity only sets the
freezing point: no density-driven overturning.

### 4.4 Hydrology (fine grid)

Precipitation minus evaporation gives runoff, routed along steepest descent.
Accumulated flow above a threshold is a river. Depressions fill (priority-flood)
and become lakes when inflow exceeds evaporation; basins with no outlet make
salt lakes or pans. Imported rivers and lakes follow their layer states (§4.9).

### 4.5 Climate classification and vegetation

Rolling multi-year monthly means give a Köppen–Geiger classification. Vegetation
is simulated per land tile: it starts from the imported layer (or is generated),
moves towards what the climate supports over decades to centuries, and feeds
back on albedo, roughness and soil water.

### 4.6 Time stepping and speeds

Every system runs every N base steps. `dt` comes from the grid spacing, radius
and maximum expected wind (advective CFL with a margin), is rounded to divide
the solar day, gives at least 24 steps per solar day, and stays fixed for a
given world. Speed controls change steps per real second, never `dt`.

| System | Rate |
|---|---|
| Dynamics, surface fluxes, clouds, precipitation | every base step |
| Radiation | every ~3 h |
| Ocean surface heat, sea ice | every base step |
| Ocean transport | with radiation, ~3 h |
| Hydrology | daily |
| Climate statistics | accumulated daily; Köppen monthly; vegetation on long-run steps |

Speed 1 at `dt` = 30 min is 48 steps/s; speed 2 (~1 week/s) is about 340 steps/s.

### 4.7 Spin-up and long runs

An accelerated spin-up runs a few sim-years in the background. Long runs
alternate a few GCM years with a century or more of vegetation: about 30 cycles
of 3 years covers 10,000 years. With a fixed orbit and no ice-sheet dynamics a
planet settles in a few thousand years, so about 10,000 is the design target.

### 4.8 Importing worlds

`.fmsmap` files, global heightmaps (equirectangular PNG/GeoTIFF), and presets
(shipped `.fmsmap` files). Also a blank-ocean start and an idealized continent.
Terrain is not edited here.

### 4.9 Layer states

Each imported layer the climate could also simulate is locked, initial (default:
started from, then free to evolve) or absent (generated). See
[FMSMAP.md](FMSMAP.md).

### 4.10 Viewer

View only, shared with `fms-geo`. Globe and flat projections; map modes for
terrain, contours, temperature, precipitation, wind particles, currents, sea ice
and snow, clouds, Köppen, vegetation and any field; date, season and speed HUD;
Play/Pause; click-to-inspect with annual graphs.

### 4.11 Persistence

Results for other modules go into `.fmsmap`. Resume saves are a climate-only
versioned binary format (planet, persistent fields, clock, RNG), with fields
tagged by name and version for migration. Same save + same build = bit-identical
sim on any thread count and any of macOS, Windows and Linux.

## 5. Validation

- **Analytic**: insolation by latitude and season; energy balance with no
  atmosphere.
- **Idealized planets**: aquaplanet Hadley circulation; no tilt, no seasons; fast
  rotation gives more, narrower cells; slow rotation widens the Hadley cell.
- **Earth benchmark**: trades, mid-latitude jets, rain-belt migration, monsoons,
  Sahara and Atacama deserts, Amazon and Congo rainforests, and Köppen agreement
  above a target.
- **Vegetation**: a drawn desert fades if unsupported and persists if supported.
- **Determinism**: two runs from one save are bit-identical.

## 6. Risks

| Risk | Mitigation |
|---|---|
| GCM numerical instability | CFL margin, diffusion, clamped parameters, idealized tests |
| Long runs drift from the live GCM | Validate against GCM-stepped runs over shorter spans |
| `.fmsmap` changes often at first | Versioned files and layers; one workspace until stable |
| Performance at fine resolution | Cell ordering, SoA, rayon, profiling; GPU port planned |
| Determinism vs. parallelism | Fixed-order reductions; tests in CI |
| No visible app until M5 | CLI PNG dumps and GIFs from M1 |
