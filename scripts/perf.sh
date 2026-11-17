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
#   scripts/perf.sh phases  <label>              # a run's time loop by phase (its timing.csv)
#   scripts/perf.sh compare <before> <after>     # two runs' phases side by side, in ms a step
#   scripts/perf.sh same    <run> <reference>    # whether two runs' results are identical, file by file
#   scripts/perf.sh top     <label> [rows]       # a saved profile as text: CPU time by function
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
usage() { sed -n '2,28p' "$0" | sed 's/^# \{0,1\}//'; exit "${1:-0}"; }

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

# The timing report a run wrote: phase,level,side_by_side,seconds,calls, the
# first row the whole loop (its calls the steps).
timing_csv() {
    local f=$OUT/$1/climate/timing.csv
    [ -f "$f" ] || die "$f is missing (a run from before the timing report?)"
    echo "$f"
}

# A run's phases: ms a step, share of the loop, calls. Indented by nesting;
# || marks phases that run side by side (their times overlap).
phases() {
    local f
    f=$(timing_csv "$1")
    awk -F, -v label="$1" '
        FNR == 1 { next }
        $1 == "loop" {
            loop = $4; steps = $5
            printf "%s: %d steps in %.2f s, %.3f ms a step, %.1f steps/s\n", label, steps, loop,
                1e3 * loop / steps, steps / loop
            printf "  %-44s %9s %7s %8s\n", "phase", "ms/step", "share", "calls"
            next
        }
        $5 == 0 { next }
        {
            name = ($3 == "true") ? "|| " $1 : $1
            for (i = 1; i < $2; i++) name = "  " name
            printf "  %-44s %9.3f %6.1f%% %8d\n", name, 1e3 * $4 / steps, 100 * $4 / loop, $5
        }' "$f"
}

# Two runs' phases side by side: ms a step in each and the change.
compare() {
    local a b
    a=$(timing_csv "$1")
    b=$(timing_csv "$2")
    awk -F, -v A="$1" -v B="$2" '
        FNR == 1 { file++; next }
        $1 == "loop" { loop[file] = $4; steps[file] = $5; next }
        {
            if (!($1 in level)) { order[++n] = $1; level[$1] = $2; side[$1] = $3 }
            ms[file, $1] = 1e3 * $4; calls[file, $1] = $5
        }
        function change(x, y) { return x > 0 ? sprintf("%+.1f%%", 100 * (y - x) / x) : "-" }
        END {
            printf "  %-44s %14s %14s %8s\n", "", A, B, "change"
            x = steps[1] / loop[1]; y = steps[2] / loop[2]
            printf "  %-44s %14.1f %14.1f %8s\n", "steps/s", x, y, change(x, y)
            x = 1e3 * loop[1] / steps[1]; y = 1e3 * loop[2] / steps[2]
            printf "  %-44s %14.3f %14.3f %8s\n", "loop (ms a step)", x, y, change(x, y)
            for (i = 1; i <= n; i++) {
                p = order[i]
                if (!calls[1, p] && !calls[2, p]) continue
                name = (side[p] == "true") ? "|| " p : p
                for (j = 1; j < level[p]; j++) name = "  " name
                x = ms[1, p] / steps[1]; y = ms[2, p] / steps[2]
                printf "  %-44s %14.3f %14.3f %8s\n", name, x, y, change(x, y)
            }
        }' "$a" "$b"
}

# Whether two runs gave the same results: their CSVs and summary, byte for
# byte (the timing report and the images aside).
same() {
    local f a b differ=0
    for f in diagnostics.csv annual_cells.csv monthly_zonal.csv zonal_mean.csv \
        circulation.csv overturning.csv rain_belt.csv ocean_heat_transport.csv summary.txt; do
        a=$OUT/$1/climate/$f
        b=$OUT/$2/climate/$f
        [ -f "$a" ] && [ -f "$b" ] || { echo "missing: $f"; differ=1; continue; }
        cmp -s "$a" "$b" || { echo "differs: $f"; differ=1; }
    done
    [ "$differ" = 0 ] && echo "$1 and $2: identical results"
    return "$differ"
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
    phases)
        [ $# -eq 1 ] || usage 1
        phases "$1"
        ;;
    compare)
        [ $# -eq 2 ] || usage 1
        compare "$1" "$2"
        ;;
    same)
        [ $# -eq 2 ] || usage 1
        same "$1" "$2"
        ;;
    top)
        [ $# -ge 1 ] && [ $# -le 2 ] || usage 1
        python3 scripts/profile_summary.py "$OUT/$1" ${2:+"$2"}
        ;;
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
