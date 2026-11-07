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
fms probe --lat 0 --lon 0 -f planet.toa_insolation --days 2 # two diurnal cycles at the equator
fms probe --lat 60 --lon -150 -f planet.daily_insolation --days 365.2422 --every 48  # a year, daily
```
