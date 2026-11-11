# Earth reference values

Observed Earth climate for comparing a run against:

    python3 scripts/compare_earth.py output/<run>

Full citations are in [REFERENCES.md](../../REFERENCES.md).

| File | Holds | Source |
|---|---|---|
| `zonal_era5.csv` | Monthly zonal means in 1° bands: skin and 2 m temperature (K), sea-level pressure (hPa), 10 m eastward wind (m/s) | ERA5, 1991–2020 (Hersbach et al. 2020) |
| `zonal_gpcp.csv` | Monthly zonal-mean precipitation (mm/day), 2.5° rows | GPCP v2.3, 1991–2020 (Adler et al. 2018) |
| `zonal_ceres.csv` | Monthly zonal means, 1° rows: incoming sunlight, absorbed sunlight, outgoing longwave, shortwave and longwave cloud effects (W/m²), cloud fraction | CERES EBAF TOA Ed4.2.1, 2001–2020 (Loeb et al. 2018) |
| `regions.csv` | Annual means over the boxes in `regions.json`: all of the above, ERA5's terrain height, NCEP's skin temperature, the albedo and the area share | All of these, with ERA5's land mask |
| `sea_ice.csv` | Monthly sea-ice area and extent per hemisphere (million km²) | NSIDC Sea Ice Index v4, G02135, 1991–2020 (Fetterer et al.) |
| `regions.json` | The region boxes, shared by the build and compare scripts | — |

## Rebuilding

Rarely needed. `scripts/fetch-reference.sh` downloads NCEP, GPCP and NSIDC into
`data/reference/raw` (gitignored). Two files need a free account and are
downloaded by hand into the same folder:

- `era5.grib`: from the Copernicus Climate Data Store, "ERA5 monthly averaged
  data on single levels from 1940 to present", product type "monthly averaged
  reanalysis", 1991–2020, all months, GRIB. Variables: skin temperature, 2 m
  temperature, mean sea-level pressure, 10 m u and v wind, total
  precipitation, sea-ice area fraction, land–sea mask, geopotential (7.3 GB at
  0.25°; only the first four and the last two are used).
- `CERES_EBAF-4.nc`: from NASA's CERES ordering page (Earthdata login), EBAF
  TOA Edition 4.2.1, monthly 1°, netCDF, with incoming solar, all-sky and
  clear-sky shortwave and longwave, and total cloud area.

Then run `.venv/bin/python scripts/build_reference.py`. The project virtualenv
needs numpy, netCDF4 and eccodes.

## Caveats

- **Reanalyses are models.** Over ice sheets and sea ice ERA5's and NCEP's
  temperatures differ by a few kelvin (the Arctic Ocean, 262 against 259 K;
  the Sahara, 300 against 297 K). `regions.csv` keeps NCEP's skin temperature
  as a measure of that spread.
- **Land and ocean are ERA5's mask** (nearest point for GPCP and CERES). The
  model uses land fraction > 0.5, so coastal boxes differ a little.
- **CERES clear sky** is measured over the cloud-free parts of each region, not
  computed with the clouds removed as the model does. Its cloud effects can
  differ by a few W/m² for that reason alone (Loeb et al. 2018). Its cloud
  fraction over snow and ice and in polar night is less certain. The global
  means match the paper's Table 5 (July 2005–June 2015) to within 0.5 W/m².
- **Sea-ice area, not extent.** The model's ice cover is an area. NSIDC leaves
  the Arctic pole hole out of area, so `build_reference.py` adds it back as ice
  (0.31 million km² to 2007 and 0.029 from 2008; user guide, Table 8).
- **Calendar months.** The model's months start at the March equinox, so
  `compare_earth.py` interpolates the reference to the model's mid-months.
- **Height adjustment.** The "ΔT at ERA5 z" column moves the model's
  temperature to ERA5's terrain height at 6.5 K/km. That is a rough guide over
  ice sheets, where surface inversions break the lapse rate.
- **The model's surface air** is its bottom layer brought down dry-adiabatically.
  Over ice sheets that misses the inversion, so compare skin temperatures there.
- **The model's wind** in the comparison is its bottom layer's, about 1 km
  deep, which should be stronger than ERA5's 10 m wind.
- **Not yet included:** surface energy fluxes, ERA5 cloud cover, winds and
  temperatures aloft (ERA5 pressure levels, M2.15), and sea-level pressure in
  the model's output.
- **CO₂.** The reference period's mean is about 380 ppm; the model's default is
  420 ppm, which is worth a few tenths of a kelvin.
