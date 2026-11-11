#!/bin/bash
# Download the raw Earth reference data into data/reference/raw (gitignored).
# From the repository root:
#   scripts/fetch-reference.sh
# Then build the committed tables with scripts/build_reference.py.
# Sources and what each file holds: reference/earth/SOURCES.md.
set -euo pipefail
cd "$(dirname "$0")/.."
RAW=data/reference/raw
mkdir -p "$RAW"

get() { # url file
  if [ -s "$RAW/$2" ]; then echo "have $2"; return; fi
  echo "get  $2"
  curl -sSfL --retry 3 -o "$RAW/$2.part" "$1"
  mv "$RAW/$2.part" "$RAW/$2"
}

# NCEP/NCAR Reanalysis 1, monthly long-term means 1991–2020, T62 Gaussian grid.
NCEP=https://downloads.psl.noaa.gov/Datasets/ncep.reanalysis.derived/surface_gauss
get $NCEP/skt.sfc.mon.ltm.1991-2020.nc ncep_skt.nc
get $NCEP/air.2m.mon.ltm.1991-2020.nc ncep_air2m.nc
get https://downloads.psl.noaa.gov/Datasets/ncep.reanalysis/surface_gauss/land.sfc.gauss.nc ncep_land.nc
get https://downloads.psl.noaa.gov/Datasets/ncep.reanalysis/surface_gauss/hgt.sfc.gauss.nc ncep_hgt.nc

# GPCP v2.3 precipitation, monthly long-term means 1991–2020, 2.5° grid.
get https://downloads.psl.noaa.gov/Datasets/gpcp/precip.mon.ltm.1991-2020.nc gpcp_precip.nc

# NSIDC Sea Ice Index v4 (G02135), monthly extent and area per hemisphere.
for h in north south; do
  H=$(echo ${h:0:1} | tr a-z A-Z)
  for m in 01 02 03 04 05 06 07 08 09 10 11 12; do
    get https://noaadata.apps.nsidc.org/NOAA/G02135/$h/monthly/data/${H}_${m}_extent_v4.0.csv \
      nsidc_${H}_${m}.csv
  done
done
echo done
