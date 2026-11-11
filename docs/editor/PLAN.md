# `fms-editor` plan

The canonical design of `fms-editor`: what it is for and the decisions that
shape it. Settled on 2026-11-11; nothing is built yet. The family it belongs to
is in [../PLAN.md](../PLAN.md) §1.1, and the file format it writes is in
[../FMSMAP.md](../FMSMAP.md).

`fms-editor` is the GUI for making and changing maps. Its first job is to
produce good `.fmsmap` inputs for the simulators. Its second, later job is to be
a map-making tool in its own right.

## 1. Goals

- **Simulator inputs first**: terrain, biomes, rivers and lakes that
  `fms-geo`, `fms-climate` and `fms-civ` can run from.
- **Cartography later**: labels, borders, styling and image export. Only
  reserved for now (§6).
- **A complete planet**: the editor sets every planet setting, so an exported
  map is runnable as it is.
- **Round trips**: results from a simulator come back into the same project for
  more editing.

Painted climate layers (temperature, rainfall) are a stretch goal.

## 2. Key decisions

| Area | Decision |
|---|---|
| Layers authored | `terrain.elevation`, biomes, rivers, lakes |
| Planet settings | All set in the editor; Earth defaults, same clamped ranges as `fms-climate` |
| Data model | Per-cell rasters for elevation and biomes; vector objects for rivers and lakes |
| Edit surface | Globe and flat projections, both views of the same sphere data |
| Raster edits | Destructive, with a stroke-history undo |
| Resolution | Grid frequency chosen per map; can be resampled up later |
| Sea | Elevation at or below zero, in the world ocean or a body of at least 0.005% of the sphere; smaller hollows are land, as in `fms-climate`. A tool shifts global sea level |
| Rivers and lakes | Never alter stored elevation; carved at export (§3.2) |
| Biomes | One fixed list shared by the family, defined in `FMSMAP.md` |
| Starting a map | Blank ocean, heightmap, `.fmsmap`, or (after v1) a generated world |
| Generator | Spherical noise shaped by a quick plate sketch; real tectonics is `fms-geo`'s |
| Saving | Native project file; Export writes a baked `.fmsmap` |
| Layer states | A lock toggle per layer, default initial, whole layers only; also a lock toggle for CO₂, which `fms-geo` can otherwise evolve ([../geo/PLAN.md](../geo/PLAN.md) §5.3) |
| Rendering | Shared render crate, used by `fms-viewer` and `fms-editor` (§4) |

## 3. Data

### 3.1 Rasters

Elevation and biomes are arrays on one icosahedral grid, at the frequency the
user chose for the map. Brushes change them directly. Brush sizes are in km.

- The sea floor is sculpted with the same brushes as land.
- Resampling up interpolates and invents no detail. Roughness comes from the
  noise brushes (§5).
- Heightmaps (PNG, GeoTIFF) come in through the `worldgen` importer.

### 3.2 Rivers and lakes

Rivers are lines and lakes are outlines with a surface level, both with points
on the sphere. They stay editable as objects, and the stored elevation is never
changed by them.

On export the editor cuts a downhill channel along each river and a basin under
each lake into the exported elevation, and also writes the river and lake
layers. The simulated hydrology then follows the drawn courses if the climate
gives them water. The editor shows the cut as a live preview.

### 3.3 Project file and export

The project file holds everything the editor needs: rasters, vector objects,
planet settings, lock toggles, undo history and (later) generator settings.
Other modules never read it.

Export bakes a `.fmsmap`: rasters as layers, rivers and lakes carved and
rasterised, each layer's state from its lock toggle.

### 3.4 Importing simulator results

The user picks which layers of a `.fmsmap` to pull into the open project, for
example `fms-geo`'s eroded elevation while keeping their own rivers.

- An imported raster replaces the project's own, as one undoable step.
- Layers the editor doesn't understand are kept and passed through on export.
- Each imported layer remembers the terrain it was computed from. After a
  terrain edit it is marked stale, and Export asks whether to drop or keep
  stale layers. Nothing is deleted silently.

## 4. Application

A new `fms-editor` binary in the family workspace, built on `eframe`/`egui` and
`wgpu`. The renderer, projections, camera and picking live in a shared library
crate used by both `fms-viewer` and `fms-editor`. The editor adds tools,
panels, and partial GPU updates for brushed cells.

The viewer is not built yet, so its renderer should be written as that library
from the start.

## 5. Stages

Building starts no earlier than `fms-climate` M3 step 0, which creates the
`.fmsmap` crate.

| Stage | Contents |
|---|---|
| First milestone | Open a blank or imported map, sculpt with undo, save a project, export a `.fmsmap` that `fms-climate` runs |
| v1 | Sculpt brushes (raise, lower, smooth, flatten, set height), biome paint and fill, layer import, planet panel |
| Generator | Noise and plate-sketch worlds from a seed |
| Rivers and lakes | Vector objects, carve at export, live preview |
| Stamps | Mountain-range, crater and island stamps; noise brushes |
| Cartography | Its own design session first (§6) |
| Painted climate | Stretch goal |

## 6. Reserved for later

- **Cartography**: the vector model and project file must be able to hold
  labels, borders and styling. Political borders come from `fms-civ`'s layers
  and name tables in `.fmsmap`; label placement and styling stay in the editor's
  project file ([../civ/PLAN.md](../civ/PLAN.md) §9.1).
- **The biome list**: drafted when `fms-climate` reaches vegetation (M3 step
  4), so both modules agree on it.
- **Plate drawing**: a tool that exports plate and crust layers; `fms-geo`
  already reads them when present.
