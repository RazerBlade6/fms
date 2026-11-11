# `fms` plan

The plan for the `fms` family as a whole. Each module has its own plan in its
folder. The file the modules exchange is described in [FMSMAP.md](FMSMAP.md).

## 1. The `fms` family

Separate modules exchanging `.fmsmap` files. Each uses only the layers it
simulates and passes the rest through. Coupling goes through files
(geo → climate → geo). `fms-civ` reads the world and never changes it.

| Module | Role | Interface |
|---|---|---|
| `fms-editor` | Make custom maps at any level of detail | GUI, edits |
| `fms-climate` | Climate and hydrosphere | Viewer only |
| `fms-geo` | Plate tectonics and erosion | Viewer only |
| `fms-civ` | Civilisation simulation and editing | GUI, edits |

Typical paths: editor → geo → climate → civ, or editor (with biomes) → civ,
which then uses a biome → typical-climate table and warns about it.

All modules stay in this Cargo workspace until the format is stable.

| Module | Plan |
|---|---|
| `fms-climate` | [climate/PLAN.md](climate/PLAN.md) |
| `fms-geo` | [geo/PLAN.md](geo/PLAN.md) |
| `fms-civ` | [civ/PLAN.md](civ/PLAN.md) |
| `fms-editor` | [editor/PLAN.md](editor/PLAN.md) |

`fms-climate` is the module being built now. Its progress is in
[climate/ROADMAP.md](climate/ROADMAP.md).

## 2. Family-wide decisions

| Area | Decision |
|---|---|
| Language | Rust |
| Platform | Desktop, cross-platform (macOS first) |
| UI and rendering | `eframe`/`egui` panels + custom `wgpu` renderer, shared by view-only modules |
| Geometry | Simulate on a sphere; flat maps are projections |
| Sphere grid | Icosahedral geodesic grid (Voronoi hexagons + 12 pentagons) |
| Determinism | Seeded RNG, fixed timestep, deterministic reductions, versioned saves |
| Portable math | Transcendentals come from `libm` via `fms-math`, enforced by `clippy.toml` |
| Interchange | `.fmsmap` files with planet settings and named layers |
| Layer states | Locked, initial (default) or absent, per layer |
