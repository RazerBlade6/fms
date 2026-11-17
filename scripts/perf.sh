#!/usr/bin/env bash
# Timing and profiling runs of `fms-climate run` for the performance pass
# (M2.15). Everything goes under output/m2.15/ (or $OUT): one folder per run,
# every command in RUNS.txt, every timing as a row of results.tsv.
#
#   scripts/perf.sh conditions                   # the machine's state, for BENCHMARKS.md
#   scripts/perf.sh spinup  <cells>              # a full run from rest; saves the restart file
#   scripts/perf.sh bench   <cells> [label]      # a timed run from the restart file
#   scripts/perf.sh sweep   <cells> [threads…]   # bench at each thread count (default 1 4 8 15)
#   scripts/perf.sh profile <cells> [label]      # the bench run under samply
#   scripts/perf.sh view    <label>              # open a saved profile in the browser
#   scripts/perf.sh results                      # the table of timings so far
#
# <cells> is 642, 2562, 10242 or 20252. Options to `fms-climate run` follow `--`.
#
# Environment:
#   WORLD=earth|aqua   Earth's heightmap (default) or a blank all-ocean world
#   THREADS=N          rayon threads (default: all cores)
#   YEARS=N            length of a bench or profile run (default 1; 2 with COLD=1)
#   COLD=1             bench or profile from rest, without a restart file or ocean jumps
#   THROTTLE=1         background priority and 4 threads, for spin-ups while the laptop is in use
#   RATE=N             samply's sampling rate in Hz (default 1000)
#   OUT=DIR            where everything goes (default output/m2.15)
set -euo pipefail

cd "$(dirname "$0")/.."
OUT=${OUT:-output/m2.15}
WORLD=${WORLD:-earth}
BIN=target/release/fms-climate
HEIGHTMAP=data/ETOPO_2022_v1_60s_N90W180_surface.tif

die() { echo "perf.sh: $*" >&2; exit 1; }
usage() { sed -n '2,24p' "$0" | sed 's/^# \{0,1\}//'; exit "${1:-0}"; }

# The climate grid frequency for a number of cells (10·f² + 2).
frequency() {
    case "$1" in
        642) echo 8 ;;
        2562) echo 16 ;;
        10242) echo 32 ;;
        20252) echo 45 ;;
        *) die "cells must be 642, 2562, 10242 or 20252, not '$1'" ;;
    esac
}

# The world's options.
world_args() {
    case "$WORLD" in
        earth)
            [ -f "$HEIGHTMAP" ] || die "$HEIGHTMAP is missing; run scripts/fetch-etopo.sh"
            echo "--heightmap $HEIGHTMAP"
            ;;
        aqua) echo "--terrain-frequency 32" ;;
        *) die "WORLD must be earth or aqua, not '$WORLD'" ;;
    esac
}

commit() {
    local c
    c=$(git rev-parse --short HEAD)
    [ -z "$(git status --porcelain -- crates Cargo.toml Cargo.lock)" ] || c="$c+changes"
    echo "$c"
}

build() { cargo build --quiet --release -p fms-climate; }

# The machine's state: what BENCHMARKS.md asks every timing to be logged with.
conditions() {
    echo "date:    $(date '+%Y-%m-%d %H:%M:%S %Z')"
    echo "commit:  $(commit)"
    echo "machine: $(sysctl -n machdep.cpu.brand_string), $(sysctl -n hw.ncpu) cores," \
        "$(($(sysctl -n hw.memsize) / 1073741824)) GB"
    echo "system:  macOS $(sw_vers -productVersion) ($(sw_vers -buildVersion)); $(rustc --version)"
    echo "power:   $(pmset -g batt | head -1 | sed "s/Now drawing from //; s/'//g")," \
        "battery $(pmset -g batt | grep -o '[0-9]*%' | head -1)," \
        "powermode $(pmset -g | awk '/powermode/ {print $2; exit}')"
    echo "load:    $(sysctl -n vm.loadavg)"
    # The second sample: the first holds averages since boot.
    top -l 2 -n 6 -o cpu -stats command,cpu | awk '/^Processes:/ {n++} n == 2' | awk '
        /^CPU usage/ {print "cpu:     " substr($0, 12)}
        /^PhysMem/   {print "memory:  " substr($0, 10)}
        busy && NF   {printf "         %s\n", $0}
        /^COMMAND/   {busy = 1; print "busiest:"}'
}

# One value from `/usr/bin/time -l`'s report in FILE: the number before LABEL.
timed() { awk -v label="$2" 'index($0, label) {print $1; exit}' "$1"; }

# run_one KIND CELLS LABEL [options…]: runs `fms-climate run` under
# `/usr/bin/time -l` in $OUT/LABEL, logs the command and adds a row of timings.
run_one() {
    local kind=$1 cells=$2 label=$3
    shift 3
    local dir=$OUT/$label threads=${THREADS:-} prefix=()
    if [ -n "${THROTTLE:-}" ]; then
        prefix=(taskpolicy -b)
        threads=${THREADS:-4}
    fi
    [ -z "$threads" ] || export RAYON_NUM_THREADS=$threads
    # Assigned first: a failure inside a command's arguments wouldn't stop the script.
    local frequency world
    frequency=$(frequency "$cells")
    world=$(world_args)
    # shellcheck disable=SC2086
    local cmd=("$BIN" run --climate-frequency "$frequency" $world "$@" --out-dir "$dir")
    mkdir -p "$dir"
    {
        echo "$(date '+%Y-%m-%d %H:%M') $(commit) $label ($kind, $cells cells, $WORLD," \
            "threads: ${threads:-all}${THROTTLE:+, background priority})"
        echo "   ${RAYON_NUM_THREADS:+RAYON_NUM_THREADS=$RAYON_NUM_THREADS }${cmd[*]}"
    } >>"$OUT/RUNS.txt"
    echo "$label: ${cmd[*]}" >&2
    { echo "== before"; conditions; } >"$dir/conditions.txt"
    local status=0
    if [ "$kind" = profile ]; then
        command -v samply >/dev/null || die "samply is missing; cargo install --locked samply"
        caffeinate -i samply record --save-only --unstable-presymbolicate \
            --rate "${RATE:-1000}" -o "$dir/profile.json.gz" -- "${cmd[@]}" \
            >"$dir/run.log" 2>"$dir/run.err" || status=$?
    else
        caffeinate -i ${prefix[@]+"${prefix[@]}"} /usr/bin/time -l "${cmd[@]}" \
            >"$dir/run.log" 2>"$dir/run.err" || status=$?
    fi
    { echo "== after"; conditions; } >>"$dir/conditions.txt"
    [ "$status" = 0 ] || { tail -20 "$dir/run.err" >&2; die "$label failed; see $dir/run.err"; }

    local speed
    speed=$(sed -n 's/.*ran in \([^ ]*\) (\([0-9.]*\) steps\/s, \([0-9.]*\) CG iterations per solve).*/\1 \2 \3/p' "$dir/run.err")
    [ -n "$speed" ] || die "no speed line in $dir/run.err"
    # shellcheck disable=SC2086
    set -- $speed
    local loop=$1 steps=$2 cg=$3 real=- user=- sys=- switches=- faults=- memory=-
    if [ "$kind" != profile ]; then
        real=$(awk '/ real / && / user / {print $1}' "$dir/run.err")
        user=$(awk '/ real / && / user / {print $3}' "$dir/run.err")
        sys=$(awk '/ real / && / user / {print $5}' "$dir/run.err")
        switches=$(($(timed "$dir/run.err" "voluntary context switches") + $(timed "$dir/run.err" "involuntary context switches")))
        faults=$(($(timed "$dir/run.err" "page reclaims") + $(timed "$dir/run.err" "page faults")))
        memory=$(($(timed "$dir/run.err" "peak memory footprint") / 1048576))
    fi
    [ -f "$OUT/results.tsv" ] || printf 'label\tkind\tcells\tworld\tthreads\tcommit\tsteps_per_s\tcg_per_solve\tloop\tprocess_s\tuser_s\tkernel_s\tcontext_switches\tpage_faults\tpeak_mb\n' >"$OUT/results.tsv"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$label" "$kind" "$cells" "$WORLD" "${threads:-all}" "$(commit)" "$steps" "$cg" "$loop" \
        "$real" "$user" "$sys" "$switches" "$faults" "$memory" >>"$OUT/results.tsv"
    echo "$label: $steps steps/s in the time loop ($loop), $cg CG iterations per solve;" \
        "whole process $real s, user $user s, kernel $sys s -> $dir" >&2
}

restart_file() { echo "$OUT/restart/$WORLD-$1.fmsrst"; }

# The options of a bench or profile run on CELLS cells.
bench_args() {
    if [ -n "${COLD:-}" ]; then
        echo "--ocean-jumps 0 --years ${YEARS:-2}"
    else
        local restart
        restart=$(restart_file "$1")
        [ -f "$restart" ] || die "$restart is missing; run 'scripts/perf.sh spinup $1' first, or set COLD=1"
        echo "--restart $restart --years ${YEARS:-1}"
    fi
}

results() {
    [ -f "$OUT/results.tsv" ] || die "no timings yet in $OUT"
    column -t -s "$(printf '\t')" "$OUT/results.tsv"
}

[ $# -ge 1 ] || usage 1
[ "$WORLD" = earth ] || [ "$WORLD" = aqua ] || die "WORLD must be earth or aqua, not '$WORLD'"
action=$1
shift
# Split the arguments at `--`: ours before, `fms-climate run`'s after.
ours=()
extra=()
while [ $# -gt 0 ]; do
    if [ "$1" = -- ]; then
        shift
        extra=("$@")
        break
    fi
    ours+=("$1")
    shift
done
set -- ${ours[@]+"${ours[@]}"}

case "$action" in
    conditions) conditions ;;
    results) results ;;
    view)
        [ $# -eq 1 ] || usage 1
        samply load "$OUT/$1/profile.json.gz"
        ;;
    spinup)
        [ $# -eq 1 ] || usage 1
        frequency "$1" >/dev/null
        build
        mkdir -p "$OUT/restart"
        run_one spinup "$1" "spinup-$WORLD-$1" --save-restart "$(restart_file "$1")" ${extra[@]+"${extra[@]}"}
        ;;
    bench | profile)
        [ $# -ge 1 ] && [ $# -le 2 ] || usage 1
        frequency "$1" >/dev/null
        build
        args=$(bench_args "$1")
        # shellcheck disable=SC2086
        run_one "$action" "$1" "${2:-$action-$WORLD-$1}" $args ${extra[@]+"${extra[@]}"}
        [ "$action" = bench ] || echo "view it with: scripts/perf.sh view ${2:-$action-$WORLD-$1}" >&2
        ;;
    sweep)
        [ $# -ge 1 ] || usage 1
        cells=$1
        shift
        frequency "$cells" >/dev/null
        build
        args=$(bench_args "$cells")
        for t in "${@:-1 4 8 15}"; do
            for threads in $t; do
                # shellcheck disable=SC2086
                THREADS=$threads run_one bench "$cells" "sweep-$WORLD-$cells-t$threads" $args ${extra[@]+"${extra[@]}"}
            done
        done
        results
        ;;
    -h | --help | help) usage ;;
    *) usage 1 ;;
esac
