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
| `fms dump --all` | Same, for every cell field |
| `fms animate [-f NAME]` | Run and write `output/NAME.gif`; defaults to a year of `planet.daily_insolation` |

### Time

The world starts at the northern vernal equinox, at midnight at longitude 0.

- `dump --days D --hours H` runs `D` days plus `H` hours before rendering.
- `animate --days D --frames N` spreads `N` frames evenly over `D` days. The
  default span is one planet year.

The legend's date stamp is the time the fields were computed, which is the last
completed step.

### Images

| Option | Default | Meaning |
|---|---|---|
| `--width W` | 1024 | Width in pixels; the map is `W × W/2`, plus a legend strip |
| `--range MIN,MAX` | data range | Fixed colour range (GIFs use one range for all frames) |
| `--smooth` | off | Blend between cell centres instead of drawing each cell flat. Prettier, but hides the true resolution. Ignored for categorical fields |
| `--out-dir DIR` | `output` | Where files go |
| `--fps N` | 12 | GIF frame rate (`animate` only) |

Colours come from each field's colour-scale hint. Elevation uses separate ocean
and land ramps, with sea level as a hard boundary. Non-finite values are magenta.

### World size

| Option | Default | Cells |
|---|---|---|
| `--terrain-frequency` | 160 | 256k (400 gives 1.6M) |
| `--climate-frequency` | 45 | 20k |
| `--steps-per-day` | 48 | 30-minute base step on Earth |

## Examples

```sh
fms dump --all --hours 12                                   # every field at noon on day 1
fms dump -f planet.daily_insolation --days 92 --smooth      # June solstice
fms animate                                                 # a year of daily-mean sunlight
fms animate -f planet.toa_insolation --days 1 --frames 48   # one day of the terminator moving
```
