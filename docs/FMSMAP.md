# The `.fmsmap` file

Every module in the `fms` family reads and writes `.fmsmap` files; they are the
only way modules talk to each other (see [PLAN.md](PLAN.md) §1.1). A file holds
what other modules can use: a planet and its layers.

**Status:** design. Version 0, with only the layers `fms-climate` needs, is built
in M3 step 0 ([ROADMAP.md](ROADMAP.md)). The full layer
list is decided when `fms-editor` starts. The binary layout is decided in M3
step 0.

## Contents

- **Header:** a magic number and the format version.
- **Planet settings:** radius, gravity, rotation, orbit, CO₂ and surface
  pressure. A module can't simulate a planet from layers alone.
- **Grids:** each grid a layer uses, given by its frequency and mesh version. The
  meshes are rebuilt deterministically from those, and a changed mesh is detected
  rather than silently scrambling old maps.
- **Layers:** named arrays, one value per cell of a grid.

Resume state (winds, ocean temperatures, and so on) is **not** in the file. A
module that can resume a run keeps its own save file for that.

## Layers

Each layer records:

| Field | Meaning |
|---|---|
| Name | `namespace.name` in snake_case, as in the field registry (`terrain.elevation`) |
| Type and units | The element type and physical units |
| Grid | Which grid it lives on (terrain or climate level) |
| Version | The layer's own format version, for migration |
| Written by | The module that wrote it, and its version |
| State | Locked, initial or absent (below) |
| Computed from | Optional: the layers a derived layer was computed from, each with a content hash |
| Terrain | Optional: a content hash of the `terrain.elevation` the layer was made for |

*Computed from* lets a module tell that a derived layer is stale. `fms-editor`
uses it to flag climate layers after the terrain under them is edited
([editor/PLAN.md](editor/PLAN.md) §3.4).

*Terrain* lets a module tell that a layer no longer sits on the land under it. A
module that changes terrain in place (erosion) updates it on the layers that
still fit. A module that moves terrain (`fms-geo` with tectonics) leaves it,
first stamping any layer that lacks one with the terrain it arrived on. Readers
warn about or ignore a layer whose *Terrain* does not match
([geo/PLAN.md](geo/PLAN.md) §6.3).

## Layer states

Each layer a module could also simulate says how that module should treat it:

| State | Meaning |
|---|---|
| **Initial** (default) | The module starts from it and lets it evolve |
| **Locked** | The module keeps it exactly as it is |
| **Absent** | The file has no such layer, so the module generates it |

A drawn desert is an initial vegetation layer, so it fades in `fms-climate` if the
climate doesn't support it. Someone who wants a world exactly as drawn takes it
from `fms-editor` straight to `fms-civ`, without running a simulation; locked
exists for cases where only some layers should change. States apply to whole
layers.

The CO₂ planet setting also carries a state, locked or initial (default),
because `fms-geo` can evolve it ([geo/PLAN.md](geo/PLAN.md) §5.3). `fms-geo`
also rewrites the ocean-fraction setting to match the terrain it exports.

## Rules for modules

- Declare which layers are required and which optional. Fail with a clear message
  naming anything required that's missing.
- Pass through every layer you don't understand, unchanged. A chain of modules
  must never lose another module's data.
- Update *written by* on each layer you change.

## Version 0 layers

| Layer | Required by `fms-climate` | Notes |
|---|---|---|
| Planet settings | Yes | Header section |
| `terrain.elevation` | Yes | m above sea level, on the terrain grid |
| `land.vegetation` | No | Initial vegetation; generated if absent |

## Layers decided for `fms-editor`

Settled in [editor/PLAN.md](editor/PLAN.md); names and encodings are fixed when
each is built.

| Layer | Notes |
|---|---|
| Biomes | One fixed list of about 12–16 biomes shared by the whole family. `fms-climate` maps it to and from its vegetation, and `fms-civ`'s biome → climate table keys on it. The list is drafted at `fms-climate` M3 step 4 |
| Rivers | Rasterised from the editor's drawn rivers. The exported elevation already has their channels cut |
| Lakes | Rasterised from the editor's drawn lakes, with surface level. The exported elevation already has their basins cut |

Vector objects, undo history and other editing data stay in the editor's
project file, not here.

## Layers decided for `fms-geo`

Settled in [geo/PLAN.md](geo/PLAN.md); names and encodings are fixed when each
is built.

| Layer | Notes |
|---|---|
| Plates | Plate id and plate velocity |
| Crust | Type, age and thickness; boundary type |
| Rock types | One fixed list of about 8–10 classes shared by the whole family, plus sediment thickness |
| Volcanism and hotspots | Per-cell rasters |
| Seismic hazard | Per-cell raster |

Mineral resources are decided with `fms-civ`.

## Climate without a climate run

`fms-civ` needs temperature and rainfall. When a map comes straight from
`fms-editor` with painted biomes and no climate layers, it falls back on a
biome → typical climate lookup table, and warns that results may be exaggerated
or wrong without real climate data. Running `fms-climate` briefly replaces the
table's values with simulated ones.
