# Headless CLI (`fms`)

`fms` builds a world, runs it with no GUI, and renders fields as
equirectangular PNGs and GIFs. Output goes to `output/` by default (gitignored).
Use a release build: grid construction is ~10× slower in debug.

```sh
cargo build --release -p fms-cli
./target/release/fms --help
```

## Commands

| Command | What it does |
|---|---|
| `fms fields` | List every field: grid, element, units, size, description |
| `fms dump -f NAME [-f NAME …]` | Run, then write `output/NAME.png` for each field |
| `fms dump --all` | Same, for every field |
| `fms animate [-f NAME]` | Run and write `output/NAME.gif`; defaults to a year of `planet.daily_insolation` |
| `fms patterns` | Write analytic test fields to `output/patterns/`: solid-body and Rossby–Haurwitz winds, vorticity, stream function and an edge scalar. Checks that rendering of edge, vertex and vector fields is right |
| `fms probe --lat LAT --lon LON -f NAME [-f NAME …]` | Record fields at one location over time; write a CSV and a line chart to `output/probe_<LOCATION>.{csv,png}` |
| `fms williamson --case N` | Run Williamson shallow-water test case 2, 5 or 6 on the climate grid; write a GIF of the surface height and final-state PNGs to `output/williamsonN/` |
| `fms baroclinic` | Run the two-layer dry core on a sheared jet seeded with small waves; write a GIF of surface pressure and final-state PNGs to `output/baroclinic/` |

### Cell, edge and vertex fields

Fields are drawn per cell, so other elements are resampled first. The legend
title says which:

| Field | Drawn as |
|---|---|
| Cell | as is |
| Vertex | kite-area-weighted average onto cells (`… AT CELLS`) |
| Edge scalar | length-weighted average of each cell's edges (`… AT CELLS`) |
| Edge normal component (winds, currents) | vector rebuilt at cell centres (Perot's method); colour shows speed (`… SPEED`) and arrows show direction |

### Time

The world starts at the northern vernal equinox, at midnight at longitude 0.

- `dump --days D --hours H` runs `D` days plus `H` hours before rendering.
- `animate --days D --frames N` spreads `N` frames evenly over `D` days. The
  default span is one planet year.

The legend's date stamp is the time the fields were computed, which is the last
completed step. It ends with the northern hemisphere's astronomical season
(`N SPRING` …), which changes at the equinoxes and solstices.

After simulating, `dump` and `animate` print a timing report to stderr:
throughput (steps/s and simulated days/s), then each system's run count, total
time, mean time per run and share of all system time.

### Images

| Option | Default | Meaning |
|---|---|---|
| `--width W` | 1024 | Width in pixels; the map is `W × W/2`, plus a legend strip |
| `--range MIN,MAX` | data range | Fixed colour range (GIFs use one range for all frames) |
| `--smooth` | off | Blend between cell centres instead of drawing each cell flat. Prettier, but hides the true resolution. Ignored for categorical fields |
| `--out-dir DIR` | `output` | Where files go |
| `--fps N` | 12 | GIF frame rate (`animate` only) |
| `--vectors FIELD` | none | Overlay arrows from a vector field on any map, e.g. wind over temperature |
| `--no-arrows` | off | Draw no arrows, even for vector fields |
| `--arrow-spacing PX` | 28 | Approximate distance between arrows. Arrows thin out towards the poles so they stay evenly spaced on the sphere, and their length scales with speed (the legend gives the full-length speed) |

Colours come from each field's colour-scale hint. Elevation uses separate ocean
and land ramps, with sea level as a hard boundary. Non-finite values are magenta.

### World size

| Option | Default | Cells |
|---|---|---|
| `--terrain-frequency` | 160 | 256k (400 gives 1.6M) |
| `--climate-frequency` | 45 | 20k |
| `--steps-per-day` | 48 | 30-minute base step on Earth |

### Terrain

Without a heightmap the world is a blank ocean, 4000 m deep everywhere.
`--heightmap FILE` (accepted by every command) makes a global equirectangular
heightmap the terrain:

| Option | Meaning |
|---|---|
| `--heightmap FILE` | A GeoTIFF (`.tif`, `.tiff`) or greyscale PNG (`.png`) covering the whole globe, north up. Values are metres above sea level, negative below |
| `--heightmap-range MIN,MAX` | PNG only: the elevations (m) of black and white. 8- and 16-bit PNGs both work |

For Earth, `scripts/fetch-etopo.sh` downloads ETOPO 2022 (60 arc-second, ice
surface, about 470 MB) to `data/`, which is gitignored:

```sh
./scripts/fetch-etopo.sh
fms dump -f terrain.elevation --heightmap data/ETOPO_2022_v1_60s_N90W180_surface.tif
```

GeoTIFFs must be in latitude/longitude. Both pixel-is-area and pixel-is-point
(grid-registered) files work, and a repeated column at +180° is dropped. Plain
TIFFs with no georeferencing are read like a PNG. Files much finer than the grid
are box-averaged while reading, down to about 4 pixels per terrain-cell spacing,
which keeps memory low. Each terrain cell then takes the area-weighted mean of
the pixels inside it, or a bilinear sample at its centre if the file is so
coarse that no pixel centre falls inside.

The importer sets `terrain.elevation`. These fields are derived from it:

| Field | Grid | Meaning |
|---|---|---|
| `terrain.slope` | terrain | Surface slope (m per m; the sea surface is flat) |
| `terrain.land_fraction` | climate | Fraction of the cell above sea level |
| `terrain.surface_height` | climate | Mean surface height: land, or sea level over ocean |
| `terrain.surface_height_stddev` | climate | Spread of surface height inside the cell (sub-grid mountains) |
| `terrain.mean_slope` | climate | Mean slope of the terrain cells inside |

### Shallow-water test cases

`fms williamson` runs the single-layer shallow-water model (the core the GCM is
built from) on the climate grid, using the standard test cases of Williamson
et al. (1992). It ignores the planet: the cases fix Earth's radius, rotation
and gravity.

| Case | Flow | Standard length | What to look for |
|---|---|---|---|
| 2 | Steady zonal flow in geostrophic balance | 5 days | Nothing should change; `sw.height_error` shows the error against the exact solution |
| 5 | Zonal flow over a 2 km conical mountain | 15 days | A Rossby wave train downstream of the mountain (30° N, 90° W) |
| 6 | Rossby–Haurwitz wave, wavenumber 4 | 14 days | The pattern moves east at ~11.6°/day without changing shape |

| Option | Default | Meaning |
|---|---|---|
| `--case N` | required | 2, 5 or 6 |
| `--alpha DEG` | 0 | Case 2: tilt of the flow relative to the poles (45 sends it over the pentagons) |
| `--days D` | the case's standard length | How long to run |
| `--frames N`, `--fps N` | 61, 12 | GIF frames (including the start) and speed |
| `--semi-implicit` | off | Step with semi-implicit leapfrog (implicit gravity waves, RAW filter) instead of RK4. The step is then limited only by the wind and the Coriolis term: 5–7× longer |
| `--courant C` | 0.5 (RK4), 0.6 (semi-implicit) | Time step as a fraction of the stability limit. For RK4 the limit is gravity waves plus wind; semi-implicit uses the case's peak wind over the run |

Outputs: `sw.surface_height.gif` (with wind arrows), and PNGs of the final
`sw.surface_height`, `sw.vorticity`, `sw.wind` and, for case 2,
`sw.height_error`. The run prints the time step, speed, and the relative change
in total mass (round-off) and energy (a tiny change from the time stepping).
Semi-implicit runs also print the mean conjugate-gradient iterations per step.
Use `--climate-frequency` for resolution, and a small `--terrain-frequency`
(e.g. 90) to skip building the unused fine grid.

### Baroclinic storms

`fms baroclinic` runs the two-layer dry dynamical core (semi-implicit, with
∇⁴ hyperdiffusion) on the climate grid. It starts from a zonal jet in
gradient and thermal-wind balance (faster aloft, colder poles), seeds zonal
waves 4–8 at 45° N and S, and lets baroclinic instability grow them into
mid-latitude storms. They become visible after about 5 days, peak around
day 17, and then slowly decay. There is no surface friction yet (that comes
with Held–Suarez), so the eddies drive lower-layer westerlies of 40+ m/s,
and polar surface pressure falls to ~840 hPa to balance them.

| Option | Default | Meaning |
|---|---|---|
| `--shear U` | 30 | Upper minus lower layer wind at the equator (m/s). Unstable above ~11 m/s |
| `--stability K` | 40 | Upper minus lower layer potential temperature (K) |
| `--seed K` | 0.5 | Amplitude of each seeded wave (K) |
| `--days D` | 30 | Simulated days |
| `--frames N`, `--fps N` | 121, 12 | GIF frames and speed |
| `--max-wind U` | 120 | Design maximum wind for the time step (m/s) |
| `--diffusion-order N` | 4 | Hyperdiffusion order: 4 (∇⁴), 6 or 8 |
| `--e-folding H` | auto | Damping time of the shortest resolved wave (hours). Default: the strongest the step allows with a 2× margin (about 5 h at 20k cells). See [OBSERVATIONS.md](OBSERVATIONS.md) O-17 |
| `--inviscid` | off | No hyperdiffusion: fronts and vortices collapse to the grid scale |

Outputs in `output/baroclinic/`:
- `dyn.surface_pressure.gif`, with lower-layer wind arrows;
- PNGs of the final `dyn.surface_pressure`, `dyn.theta_lower`,
  `dyn.vorticity_lower`, `dyn.wind_upper` and `dyn.wind_lower`;
- `diagnostics.csv`, four samples per day. Columns: total mass, θ, kinetic
  and total energy, the energy removed by diffusion so far, the relative
  drifts, the energy-budget residual (energy change not explained by
  diffusion, i.e. time-stepping error), max wind, and the surface-pressure
  and θ ranges.

The state is checked after every step. The run stops if any value isn't
finite, a wind passes the step's Courant limit (`--max-wind` / 0.6, about
200 m/s by default), or surface pressure or θ leaves generous physical
bounds. The error says what failed, where and when. A dump goes to
`output/baroclinic/blowup/`:
- `report.txt`: the problem, its location, the run settings, the last
  healthy state's extremes, and the last 12 diagnostics samples;
- `diagnostics.csv`;
- PNGs of the last healthy state.

The GIF up to that point is written as usual.

### Point probes

`fms probe` reads each field at the cell nearest the requested point, on that
field's own grid. It reads the cell's value without interpolating, so you see
exactly what the model holds there. Non-cell fields are resampled to that cell
the same way maps are, so a probe and a map always agree. Vector fields give
two columns, `NAME.east` and `NAME.north`.

| Option | Default | Meaning |
|---|---|---|
| `--lat`, `--lon` | required | Location in degrees; south and west are negative |
| `--days D`, `--hours H` | 1 day | How long to record |
| `--every N` | 1 | Record every N base steps, starting with the first |
| `--width W` | 1024 | Chart width in pixels |
| `--out-dir DIR` | `output` | Where files go |

Files are named after the location, e.g. `probe_33.9S_18.4E.csv`. The CSV has
one row per sample: `step`, `time_s`, `day` (solar days since `t = 0`), `year`,
`day_of_year` (from 1), `hour`, then one column per field, headed
`NAME [units]`. Values are written in full precision (Rust's shortest
round-trip form). A sample taken after step `s` holds the fields computed for
time `s · dt`. The chart has one panel per column on a shared time axis
(hours for runs of up to 3 days, otherwise days). Each panel title gives the
grid cell that was read and where its centre is. `fms probe` also prints each
column's min, max and mean.

Sampling every N steps always lands on the same time of day when N divides
the day (e.g. `--every 48` at 48 steps per day). That's useful for seasonal
curves, but a diurnal signal then shows up as a fixed offset. Use `--every 1`
or a daily-mean field to see the whole cycle.

## Examples

```sh
fms dump --all --hours 12                                   # every field at noon on day 1
fms dump -f planet.daily_insolation --days 92 --smooth      # June solstice
fms animate                                                 # a year of daily-mean sunlight
fms animate -f planet.toa_insolation --days 1 --frames 48   # one day of the terminator moving
fms patterns --smooth                                       # rendering self-test
fms dump -f terrain.elevation -f terrain.surface_height_stddev \
    --heightmap data/ETOPO_2022_v1_60s_N90W180_surface.tif   # Earth's terrain on both grids
fms probe --lat 0 --lon 0 -f planet.toa_insolation --days 2 # two diurnal cycles at the equator
fms probe --lat 60 --lon -150 -f planet.daily_insolation --days 365.2422 --every 48  # a year, daily
fms williamson --case 5 --terrain-frequency 90             # flow over a mountain, 15 days
fms baroclinic --terrain-frequency 90                      # a storm life cycle on a jet, 30 days
```
