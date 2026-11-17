#!/usr/bin/env python3
"""A text summary of a samply profile saved by `scripts/perf.sh profile`.

    scripts/profile_summary.py output/m2.15/<label> [rows]

Reads `profile.json.gz` and its symbols (`profile.json.syms.json`) and prints
CPU time (each sample weighted by the CPU time its thread used since the last
one): by library, by function (self time, then inclusive time), and for the
samples that end in the system's libraries, by the nearest function of ours
above them: where the kernel time comes from.
"""

import bisect
import collections
import gzip
import json
import re
import sys


def load(folder):
    profile = json.load(gzip.open(f"{folder}/profile.json.gz"))
    sidecar = json.load(open(f"{folder}/profile.json.syms.json"))
    strings = sidecar["string_table"]
    symbols = {}
    for lib in sidecar["data"]:
        addresses = [a for a, _ in lib["known_addresses"]]
        entries = [lib["symbol_table"][i] for _, i in lib["known_addresses"]]
        symbols[lib["debug_name"]] = (addresses, entries)
    return profile, strings, symbols


HASH = re.compile(r"::h[0-9a-f]{16}$")
GENERIC = re.compile(r"<[^<>]*>")


def tidy(name):
    """A Rust symbol without its hash and with generics collapsed."""
    name = HASH.sub("", name)
    for _ in range(8):
        shorter = GENERIC.sub("<…>", name)
        if shorter == name:
            break
        name = shorter.replace("<…>", "\0")
    return name.replace("\0", "<…>")


def main():
    folder = sys.argv[1]
    rows = int(sys.argv[2]) if len(sys.argv) > 2 else 40
    profile, strings, symbols = load(folder)
    libs = [lib["debugName"] for lib in profile["libs"]]
    ours = "fms-climate"

    by_lib = collections.Counter()
    self_time = collections.Counter()
    inline_self = collections.Counter()
    inclusive = collections.Counter()
    system_from = collections.Counter()
    system_leaf = collections.Counter()
    total = 0.0

    for thread in profile["threads"]:
        frames, funcs = thread["frameTable"], thread["funcTable"]
        resources, stacks = thread["resourceTable"], thread["stackTable"]

        def describe(frame):
            """(library, symbol, innermost inlined function) of a frame."""
            resource = funcs["resource"][frames["func"][frame]]
            lib = libs[resources["lib"][resource]] if resource >= 0 else "?"
            address = frames["address"][frame]
            if lib in symbols:
                addresses, entries = symbols[lib]
                i = bisect.bisect_left(addresses, address)
                if i < len(addresses) and addresses[i] == address:
                    entry = entries[i]
                    symbol = tidy(strings[entry["symbol"]])
                    inlined = entry.get("frames") or []
                    inner = tidy(strings[inlined[0]["function"]]) if inlined else symbol
                    return lib, symbol, inner
            return lib, f"{lib}+{address:#x}", f"{lib}+{address:#x}"

        described = {}

        def cached(frame):
            if frame not in described:
                described[frame] = describe(frame)
            return described[frame]

        samples = thread["samples"]
        for stack, cpu in zip(samples["stack"], samples["threadCPUDelta"]):
            if stack is None or not cpu:
                continue
            total += cpu
            lib, symbol, inner = cached(stacks["frame"][stack])
            by_lib[lib] += cpu
            self_time[symbol] += cpu
            inline_self[inner] += cpu
            seen = set()
            nearest = None
            s = stack
            while s is not None:
                l, sym, _ = cached(stacks["frame"][s])
                if sym not in seen:
                    seen.add(sym)
                    inclusive[sym] += cpu
                if nearest is None and l == ours and not sym.startswith(
                    ("rayon", "<rayon", "std::", "core::", "alloc::", "crossbeam", "<crossbeam", "__rust", "__rdl")
                ) and "rayon_core" not in sym and "rayon::" not in sym:
                    nearest = sym
                s = stacks["prefix"][s]
            if lib != ours:
                system_from[nearest or "(the thread pool, idle)"] += cpu
                system_leaf[symbol] += cpu

    def table(title, counter, n=rows):
        print(f"\n{title}")
        for name, cpu in counter.most_common(n):
            print(f"  {cpu / 1e6:9.2f} s {100 * cpu / total:6.2f}%  {name[:150]}")

    print(f"{folder}: {total / 1e6:.1f} s of CPU time sampled")
    table("By library (where the sample ended)", by_lib)
    table("Self time by function", self_time)
    table("Self time by innermost inlined function", inline_self)
    table("Inclusive time by function", inclusive, rows + 20)
    table("Samples ending outside our binary, by system function", system_leaf, 25)
    table("The same, by our nearest function above", system_from, 25)


if __name__ == "__main__":
    main()
