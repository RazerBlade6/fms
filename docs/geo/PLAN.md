# `fms-geo` plan

The canonical design of `fms-geo`: what it is for and the decisions that shape
it. Settled on 2026-11-12; nothing is built yet. The family it belongs to is in
[../PLAN.md](../PLAN.md) §1.1, and the file format it reads and writes is in
[../FMSMAP.md](../FMSMAP.md).

`fms-geo` simulates a planet's geology: plates move, crust is made and
destroyed, mountains rise and wear down, and the result is terrain
`fms-climate` can run on.

## 1. Goals

- **Whole worlds from a seed**: about 1 Gyr of history, starting from scattered
  proto-continents.
- **Drawn maps too**: what geo does with existing terrain follows the elevation
  layer's state (§6.1).
- **Emergent tectonics**: rifting, plate reorganisation and supercontinent
  cycles come out of the forces on the plates, not from scripts.
- **Watchable drift**: a run plays live in the shared viewer and can be exported
  at any moment. No rewind; the same seed and settings reproduce any epoch.
- **Terrain for climate**: valleys at terrain resolution, with sea level, ocean
  fraction and CO₂ consistent with the exported world.
- **Geology for other modules**: plates, crust, rock types, volcanism and
  seismic hazard as layers.
- **View only**: seed and settings fully determine a world; there is no
  steering.

## 2. Key decisions

| Area | Decision |
|---|---|
| Tectonic model | Rigid plates rotating on the sphere, driven by slab pull, ridge push, mantle drag and collision resistance. No mantle simulation |
| Regime | Earth-like moving plates only |
| Crust storage | Each plate carries its own crust samples, which rotate rigidly with it; they are resampled onto the fixed icosahedral grid for boundaries, erosion, viewing and export |
| Elevation | Not stored; it follows from crust thickness, density and sea-floor age by isostasy |
| Start state | Scattered proto-continents covering a chosen fraction, on seeded plates |
| History | About 1 Gyr by default, in steps of about 1 Myr; user-set |
| Run time | About an hour for a full history at default resolution |
| Grids | A mid grid (~160k cells, ~55 km) for plates and the coarse erosion budget; the terrain grid (1–2.5M cells) for the fine pass |
| Surface processes | Rivers and hillslopes, sediment deposition, glacial carving |
| Sea level | A fixed volume of water fills the basins |
| Climate estimate | A small energy-balance model for temperature plus rule-based winds and rain shadows. Real `fms-climate` layers are used only while the terrain still matches them |
| CO₂ | Set by the balance of volcanic outgassing and rock weathering, when its state allows (§5.3) |
| Planet settings | Radius and gravity change the geology; orbital and stellar settings feed only the climate estimate |
| Drawn maps | Plates and crust are inferred from the elevation unless the file supplies them |
| Outputs | Elevation, plates and crust, rock types, volcanism, seismic hazard. No river or lake layers |
| Other modules' layers | Passed through unchanged, and marked stale when tectonics has moved the terrain |
| Interface | View only; settings fixed at start |
| Persistence | Results in `.fmsmap`; its own exact resume save |
| Determinism | The same guarantee as `fms-climate` |
| Routing code | One crate shared with `fms-climate`'s hydrology |
| Validation | Earth statistics and idealized cases |
| First milestone | Headless, plates first |

## 3. Tectonics

### 3.1 Plates and crust

Plates are rigid; each has a rotation about a pole. Each crust sample carries
type (continental or oceanic), thickness, age, rock type and sediment
thickness. Each step: plates rotate, crust is resampled onto the mid grid, and
overlaps and gaps between plates locate the boundaries. Because samples ride on
their plate, coasts and ranges drift for a billion years without blurring.

### 3.2 Forces

Each step every plate's rotation comes from a balance of torques: slab pull from
its subducting edges, ridge push, drag from the mantle beneath it, and
resistance where continents collide. A *vigour* setting scales the driving
forces. A plate under enough tension rifts in two; plates welded by a collision
merge. Exact force laws are fixed when built.

### 3.3 Boundaries

| Boundary | What happens |
|---|---|
| Spreading ridge | The gap fills with new oceanic crust of age zero; the sea floor sinks as it ages |
| Subduction | The denser plate sinks, making a trench and a volcanic arc. Arcs build new continental crust; some sediment is carried down |
| Continental collision | Continental crust cannot sink, so it thickens into ranges and plateaus, and the plates suture |
| Rift | A continent under tension thins and splits, opening a new ocean |
| Transform | Plates slide past each other; earthquakes, little relief |
| Arc or fragment arrival | Island arcs and continental fragments weld onto the continent they reach |

### 3.4 Hotspots

Seeded mantle plumes stay fixed while plates move over them, leaving island
chains and flood basalts.

### 3.5 Isostasy and sea level

Elevation is how high each column of crust floats. Collisions thicken crust and
raise it; erosion thins it and it rebounds; ocean floor deepens with age.
Gravity limits how tall mountains stand.

The planet has a fixed volume of water (a setting), and sea level is wherever
that fills the basins: young shallow oceans flood the continents and old deep
ones drain them. Exported elevation is relative to that level, the
ocean-fraction planet setting is rewritten to match, and what counts as sea
follows the family rule ([../editor/PLAN.md](../editor/PLAN.md) §2).

## 4. Surface processes

### 4.1 Two levels

Every step, on the mid grid: a coarse erosion and sediment budget that wears
ranges down and fills basins and shelves. Before each export, on the terrain
grid: a fine pass that carves the last few Myr of valleys at terrain
resolution, driven by the uplift from the mid-grid state. Landscapes forget
their valleys within about 10 Myr, so older ones need not be carried. The live
view shows the mid-grid state; an export waits for the fine pass.

### 4.2 Processes

- **Rivers**: incision by stream power, from drainage area, rainfall, slope,
  rock erodibility and gravity.
- **Hillslopes**: diffusion.
- **Sediment**: eroded material is carried downstream and laid down as
  floodplains, deltas, shelves and basins; it becomes the sediment-thickness
  layer and sedimentary rock.
- **Glacial carving** (a late stage): fjords, U-shaped valleys and scoured
  shields where the climate estimate's long-term average ice cover sits. Ice
  only carves: it does not lower sea level or weigh down the crust, since steps
  of ~1 Myr cannot resolve single ice ages.

### 4.3 Routing

Steepest descent, depression filling and flow accumulation on the icosahedral
grid live in one small crate used by both geo's erosion and `fms-climate`'s
hydrology; whichever module reaches it first builds it. geo writes no river or
lake layers: the valleys and basins are in the elevation, and `fms-climate`
finds the rivers from real runoff.

## 5. Climate estimate and CO₂

### 5.1 The estimate

Erosion needs rainfall and glaciers need temperature, so geo estimates both.
Temperature comes from a small energy-balance model: sunlight from the planet
settings, greenhouse warming from CO₂, ice-albedo feedback, cooling with
altitude. Wind belts come from rotation and tilt. Rain comes from moisture
carried inland from the sea and dropped on windward slopes. The estimate is
refreshed every few steps as the terrain moves, is tuned once during
development against `fms-climate` runs at several CO₂ levels, and is never
written out as climate layers.

### 5.2 Real climate layers

When the file carries `fms-climate` results computed for its terrain, they
replace the estimate in runs that do not move the terrain: erosion-only runs,
including redoing the fine pass with real rainfall, which is the main use of
geo → climate → geo. Once tectonics moves the terrain, geo uses its own
estimate. There is no per-file calibration.

### 5.3 CO₂

Volcanoes at ridges, arcs and hotspots add CO₂; weathering of rock removes it,
faster where it is warm, wet and steep. Carbon settles faster than a step, so
each step CO₂ is the value at which the two balance.

The CO₂ planet setting carries a state like a layer
([../FMSMAP.md](../FMSMAP.md)). *Initial* (default) means a run with tectonics
writes the evolved value on export; *locked* means geo never changes it. Runs
without tectonics never change it. Outgassing and weathering rates are viewable
either way. CO₂ stays within the family's range; a snowball or hothouse world is
a valid result.

## 6. Files

### 6.1 Run modes and layer states

| Elevation | What geo does |
|---|---|
| Absent | Generates a world from the seed |
| Initial | Starts from the drawn terrain, with either a full run (tectonics and erosion) or erosion only, chosen per run |
| Locked | Derive only: the terrain is untouched and geology layers are inferred for it |

To start from drawn terrain geo inverts isostasy to get crust thickness, takes
crust type from elevation, and partitions plates from the seed with boundaries
biased toward trenches, ridges and ranges. If the file already has plate and
crust layers (from an earlier geo run), they are used instead. Inferred geology
is a plausible guess and is flagged as one. The grid and the water volume come
from the file. geo's own layers are always treated as initial.

### 6.2 Layers written

Names and encodings are fixed when each is built.

| Layer | Notes |
|---|---|
| Elevation | `terrain.elevation`, relative to sea level, with the fine pass applied |
| Plates | Plate id and plate velocity; the velocity lets a later run warm-start |
| Crust | Type, age and thickness; boundary type |
| Rock type | One fixed family-wide list of about 8–10 classes, e.g. basalt, granite, sandstone and shale, limestone, metamorphic, volcanic, loose sediment. Each class has an erodibility. Shared like the biome list |
| Sediment thickness | Laid down by deposition (§4.2) |
| Volcanism and hotspots | Per-cell rasters, since the format has no point features |
| Seismic hazard | Per-cell raster |

geo also updates the planet settings: ocean fraction always, CO₂ when its state
is initial.

### 6.3 Other modules' layers

Every layer geo does not simulate passes through unchanged. After a run with
tectonics those layers (painted biomes, drawn rivers and lakes, climate
results) no longer sit on the land they were made for, so they are marked stale
through the terrain marker in [../FMSMAP.md](../FMSMAP.md), and other modules
warn or ignore them. Erosion-only runs leave them valid.

### 6.4 Saves

An exported `.fmsmap` carries enough geo layers for a later run to pick up
approximately. For exact continuation geo keeps its own resume save (plates,
crust samples, clock, RNG), written at pause or at the end. The same save and
the same build continue bit-for-bit.

## 7. Application

A headless `fms-geo` binary and map modes in the shared `fms-viewer`, with the
crates in the family workspace.

Settings, fixed at the start of a run: seed, history length, starting plate
count, continental fraction, water volume, tectonic vigour, erosion strength,
terrain grid frequency.

The viewer plays the drift live with a clock in Myr, map modes for each geo
layer, pause, and export.

Determinism is as in `fms-climate`: seeded RNG, fixed step, deterministic
reductions, `fms-math` transcendentals; bit-identical on any thread count and on
macOS, Windows and Linux.

## 8. Validation

Regression tests on what a generated world should share with Earth: a
two-peaked elevation histogram, the sea-floor age distribution, the plate-size
distribution, plate speeds, the length of the supercontinent cycle, and
river-network scaling laws. Idealized cases with known answers: one subduction
zone, one collision, isostatic rebound, a single eroding block.

## 9. Stages

Building starts no earlier than `fms-climate` M3 step 0, which creates the
`.fmsmap` crate. It is headless, with PNG dumps, until the viewer stage.

| Stage | Contents |
|---|---|
| Plates | Rigid plates carrying crust samples, resampling, forces |
| Boundaries | Spreading, subduction, collision, rifting, arcs, hotspots |
| Isostasy and sea level | Elevation from crust, fixed water volume |
| Coarse erosion | Mid-grid erosion and sediment budget, rainfall rules |
| Fine pass and export | The routing crate, valleys on the terrain grid, `.fmsmap` out |
| Drawn maps | Import, plate inference, erosion-only and derive-only runs, real climate layers |
| Viewer and saves | Map modes, live drift, resume save |
| Glacial carving and CO₂ | Energy-balance temperature, ice cover, the CO₂ balance |

## 10. Reserved for later

- **Mineral resources**: decided in the `fms-civ` design session.
- **Other regimes**: stagnant-lid worlds, and the impact craters they would
  need.
- **Coasts and reefs.**
- **A recorded timeline**: scrubbing back and exporting past epochs.
- **Steering a run**, and interior planet settings such as heat and age.
- **Plate drawing** in `fms-editor`; geo already reads plate and crust layers
  when present.
- **Ice volume**: sea-level change and crust loading by ice sheets.
- **Layered stratigraphy.**
- **Running present-day Earth forward** as a demonstration.
- **Calibrating the climate estimate** per file.
