#!/usr/bin/env bash
# Run calibration variants of the climate example in parallel, a few threads
# each (one run at the default resolution can't use all cores; several at
# once get more done). Each line of the variants file is a name followed by
# environment settings, e.g.
#
#   base
#   rhc60   FMS_RHC_LOW=0.6
#   bl95    FMS_CONDENSE_BL=0.95 FMS_RHC_LOW=0.6
#
# Usage: scripts/tune.sh VARIANTS [YEARS] [K]
#   THREADS (default 4) threads per run, JOBS (default cores / THREADS) runs
#   at once, and any FMS_* settings in the environment apply to every run
#   (e.g. FMS_SST=control FMS_AVERAGE=60 for fixed-SST runs of 0.3 years).
# Logs go to output/tuning/<name>.log; a summary line per run is printed.
set -euo pipefail

variants=${1:?usage: scripts/tune.sh VARIANTS [YEARS] [K]}
years=${2:-0.3}
k=${3:-16}
threads=${THREADS:-4}
cores=$(sysctl -n hw.ncpu 2>/dev/null || nproc)
jobs=${JOBS:-$(( cores / threads > 0 ? cores / threads : 1 ))}
out=output/tuning
mkdir -p "$out"

cargo build -q --release -p fms-atmos --example climate
bin=target/release/examples/climate

run() {
    local name=$1; shift
    env RAYON_NUM_THREADS="$threads" "$@" "$bin" "$k" "$years" > "$out/$name.log" 2>&1 \
        || echo "$name: failed (see $out/$name.log)"
}

start=$(date +%s)
running=0
while read -r name settings; do
    [[ -z "$name" || "$name" == \#* ]] && continue
    # shellcheck disable=SC2086
    run "$name" $settings &
    running=$((running + 1))
    if (( running >= jobs )); then
        wait -n
        running=$((running - 1))
    fi
done < "$variants"
wait
echo "all runs in $(( $(date +%s) - start )) s ($jobs at a time, $threads threads each)"

while read -r name _; do
    [[ -z "$name" || "$name" == \#* ]] && continue
    log="$out/$name.log"
    summary=$(grep -m1 "^final" "$log" || grep -m1 "BLOW-UP" "$log" || echo "no summary")
    psi=$(grep -m1 "largest |psi|" "$log" | sed 's/.*: //' || true)
    echo "$name: $summary; psi $psi"
done < "$variants"
