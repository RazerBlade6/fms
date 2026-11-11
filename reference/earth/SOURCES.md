# Earth reference values

Observed Earth climate, 1991–2020, for comparing a run against. Compare with:

    python3 scripts/compare_earth.py output/<run>

Rebuild the tables (rarely needed) with `scripts/fetch-reference.sh` and then
`.venv/bin/python scripts/build_reference.py`. The project virtualenv needs
numpy and netCDF4. Full citations are in [REFERENCES.md](../../REFERENCES.md).

| File | Holds | Source |
|---|---|---|
| `zonal_ncep.csv` | Monthly zonal means of skin and 2 m air temperature (K), one row per grid latitude (1.9°) | NCEP/NCAR Reanalysis 1, monthly long-term means (Kalnay et al. 1996) |
| `zonal_gpcp.csv` | Monthly zonal-mean precipitation (mm/day), 2.5° rows | GPCP v2.3 monthly long-term means (Adler et al. 2018) |
| `regions.csv` | Annual means over the boxes in `regions.json`: skin and 2 m temperature, precipitation, NCEP's terrain height, and area share | Both of the above, with NCEP's land mask |
| `sea_ice.csv` | Monthly sea-ice area and extent per hemisphere (million km²) | NSIDC Sea Ice Index v4, G02135 (Fetterer et al.) |
| `regions.json` | The region boxes, shared by the build and compare scripts | — |

## Caveats

- **NCEP Reanalysis 1 is old and coarse.** Its terrain is a smooth T62
  spectral field, so regional heights differ from ETOPO's (Antarctica south of
  75°S: 2,343 m in NCEP). Its surface temperatures over the ice sheets and sea
  ice come from its own model and are uncertain by a few kelvin. ERA5 would be
  better, but it needs an account.
- **Land and ocean are NCEP's mask.** For GPCP, its nearest point is used. The
  model uses land fraction > 0.5, so coastal boxes differ a little.
- **Sea-ice area, not extent.** The model's ice cover is an area. NSIDC leaves
  the Arctic pole hole out of area, so `build_reference.py` adds it back as ice
  (0.31 million km² to 2007 and 0.029 from 2008; user guide, Table 8).
- **Calendar months.** The model's months start at the March equinox, so
  `compare_earth.py` interpolates the reference to the model's mid-months.
- **Height adjustment.** The "ΔT at ref z" column moves the model's temperature
  to the reference's terrain height at 6.5 K/km. That is a rough guide over ice
  sheets, where surface inversions break the lapse rate.
- **The model's surface air** is its bottom layer brought down dry-adiabatically.
  Over ice sheets that misses the inversion, so compare skin temperatures there.
- **Not yet included:** top-of-atmosphere radiation (CERES EBAF needs a NASA
  Earthdata login, and the Loeb et al. 2018 paper could not be fetched), cloud
  cover, winds and sea-level pressure.
- **CO₂.** The reference period's mean is about 380 ppm; the model's default is
  420 ppm, which is worth a few tenths of a kelvin.
