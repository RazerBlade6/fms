#!/bin/sh
# Download the ETOPO 2022 global relief model (60 arc-second, ice surface) into
# data/, for `fms --heightmap` and the Earth benchmark. About 470 MB.
#
# Source: NOAA NCEI, https://www.ncei.noaa.gov/products/etopo-global-relief-model
# (public domain). Re-running is safe: an existing complete file is kept.
set -eu

URL="https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO2022/data/60s/60s_surface_elev_gtif/ETOPO_2022_v1_60s_N90W180_surface.tif"
SIZE=465969062

cd "$(dirname "$0")/.."
mkdir -p data
OUT="data/$(basename "$URL")"

if [ -f "$OUT" ] && [ "$(wc -c < "$OUT" | tr -d ' ')" = "$SIZE" ]; then
    echo "$OUT is already complete"
    exit 0
fi

echo "downloading $URL"
curl --fail --location --continue-at - --output "$OUT" "$URL"

if [ "$(wc -c < "$OUT" | tr -d ' ')" != "$SIZE" ]; then
    echo "error: $OUT is not the expected $SIZE bytes" >&2
    exit 1
fi
echo "saved $OUT"
