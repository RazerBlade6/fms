"""Compare an Earth climate run against the reference in reference/earth.
From the repository root (plain Python, no packages needed):
    python3 scripts/compare_earth.py output/<run>
where <run> is the --out-dir of `fms-climate climate` (or its climate/ folder).

The model's final year is compared with 1991–2020 means. Model months start at
the March equinox, so the reference is interpolated to the model's mid-months.
Model cells are weighted equally (the icosahedral cells are close to equal
area). Sources and caveats: reference/earth/SOURCES.md."""
import csv
import json
import math
import os
import re
import sys

REF = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "reference", "earth")
EARTH_AREA = 4 * math.pi * 6371.0**2 / 1e6  # million km²
YEAR = 365.2422
EQUINOX_DAY = 79.0  # days after 1 January (≈ 20 March)
LAPSE = 6.5e-3  # K/m, to bring a region's temperature to the reference height
MONTHS = "Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec".split()


def read(path):
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def numeric(rows):
    return [{k: float(v) for k, v in r.items()} for r in rows]


def in_box(lat, lon, land_fraction, box):
    if not (box["lat"][0] <= lat <= box["lat"][1] and box["lon"][0] <= lon <= box["lon"][1]):
        return False
    if box["surface"] == "land":
        return land_fraction > 0.5
    if box["surface"] == "ocean":
        return land_fraction <= 0.5
    return True


def mean(values):
    return sum(values) / len(values) if values else float("nan")


def at_model_months(series):
    """Calendar-month means interpolated (cyclically) to the model's 12 mid-months."""
    out = []
    for m in range(12):
        day = (EQUINOX_DAY + (m + 0.5) * YEAR / 12) % YEAR
        x = day / (YEAR / 12) - 0.5  # position on the calendar-month axis
        k = math.floor(x)
        t = x - k
        out.append((1 - t) * series[k % 12] + t * series[(k + 1) % 12])
    return out


def month_label(m):
    start = EQUINOX_DAY + m * YEAR / 12
    a, b = int(start / (YEAR / 12)) % 12, int((start + YEAR / 12) / (YEAR / 12)) % 12
    return f"{MONTHS[a]}–{MONTHS[b]}"


def bands(latitudes):
    """Band edges midway between the model's zonal-mean latitudes."""
    lats = sorted(latitudes)
    edges = [-90.0] + [(a + b) / 2 for a, b in zip(lats, lats[1:])] + [90.0]
    return list(zip(lats, edges[:-1], edges[1:]))


def ref_zonal(rows, key, lo, hi):
    """Reference zonal means in a band, per calendar month, area-weighted."""
    out = []
    for m in range(1, 13):
        sel = [r for r in rows if r["month"] == m and lo <= r["latitude_deg"] < hi]
        w = sum(r["area_share"] for r in sel)
        out.append(sum(r[key] * r["area_share"] for r in sel) / w if w else float("nan"))
    return out


def fmt(x, spec):
    return format(x, spec) if math.isfinite(x) else format("—", spec.split(".")[0].lstrip("+") + "s")


def main():
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    run = sys.argv[1]
    if os.path.isdir(os.path.join(run, "climate")):
        run = os.path.join(run, "climate")
    cells = numeric(read(os.path.join(run, "annual_cells.csv")))
    monthly = numeric(read(os.path.join(run, "monthly_zonal.csv")))
    last = max(r["year"] for r in monthly)
    monthly = [r for r in monthly if r["year"] == last]
    if len({r["month"] for r in monthly}) != 12:
        sys.exit(f"{run}: the final year ({last:.0f}) has fewer than 12 months")
    boxes = {k: v for k, v in json.load(open(os.path.join(REF, "regions.json"))).items()
             if not k.startswith("_")}
    regions = {r["region"]: r for r in read(os.path.join(REF, "regions.csv"))}
    ncep = numeric(read(os.path.join(REF, "zonal_ncep.csv")))
    gpcp = numeric(read(os.path.join(REF, "zonal_gpcp.csv")))
    ice_ref = numeric(read(os.path.join(REF, "sea_ice.csv")))

    print(f"{run}: final year {last:.0f}, {len(cells)} cells, against 1991–2020\n")
    summary = os.path.join(run, "summary.txt")
    if os.path.exists(summary):
        text = open(summary, encoding="utf-8").read()
        found = re.search(r"Final-year global means: (.*)", text)
        if found:
            print(f"Model only (no reference yet): {found.group(1)}\n")

    # Regions, annual means.
    print("Regions, annual mean. T is skin temperature, Ta is near-surface air (model: bottom layer\n"
          "brought to the surface; reference: 2 m). 'at ref z' moves the model's T to the\n"
          f"reference's terrain height at {LAPSE * 1000:.1f} K/km.")
    print(f"{'region':16}{'T':>8}{'ref':>8}{'ΔT':>7}{'Ta':>8}{'ref':>8}{'ΔTa':>7}"
          f"{'rain':>7}{'ref':>7}{'ratio':>7}{'z':>7}{'ref z':>7}{'ΔT at ref z':>12}")
    for name, box in boxes.items():
        sel = [c for c in cells if in_box(c["latitude_deg"], c["longitude_deg"], c["land_fraction"], box)]
        ref = regions.get(name)
        if not sel or ref is None:
            print(f"{name:16} no cells")
            continue
        t = mean([c["rad.surface_temperature"] for c in sel])
        ta = mean([c["atm.surface_air_temperature"] for c in sel]) if "atm.surface_air_temperature" in sel[0] \
            else float("nan")
        rain = mean([c["atm.precipitation"] for c in sel])
        z = mean([c["surface_height_m"] for c in sel]) if "surface_height_m" in sel[0] else float("nan")
        rt, rta, rr, rz = (float(ref[k]) for k in
                           ("skin_temperature_k", "air_2m_k", "precipitation_mm_day", "ncep_height_m"))
        print(f"{name:16}{t:8.2f}{rt:8.2f}{t - rt:+7.2f}{fmt(ta, '8.2f')}{rta:8.2f}{fmt(ta - rta, '+7.2f')}"
              f"{rain:7.2f}{rr:7.2f}{rain / rr:7.2f}{fmt(z, '7.0f')}{rz:7.0f}"
              f"{fmt(t - rt - LAPSE * (rz - z), '+12.2f')}")
    print()

    # Zonal means, annual and seasonal.
    lats = sorted({r["latitude_deg"] for r in monthly})
    rows = []
    for lat, lo, hi in bands(lats):
        sel = sorted((r for r in monthly if r["latitude_deg"] == lat), key=lambda r: r["month"])
        share = sel[0]["area_share"]
        t = [r["rad.surface_temperature"] for r in sel]
        p = [r["atm.precipitation"] for r in sel]
        rt = at_model_months(ref_zonal(ncep, "skin_temperature_k", lo, hi))
        rp = at_model_months(ref_zonal(gpcp, "precipitation_mm_day", lo, hi))
        rows.append((lat, share, t, rt, p, rp))
    print("Zonal means: annual skin temperature (K) and rain (mm/day), and the seasonal range of\n"
          "temperature (warmest minus coldest month).")
    print(f"{'lat':>7}{'T':>8}{'ref':>8}{'ΔT':>7}{'rain':>7}{'ref':>7}{'range':>7}{'ref':>7}")
    for lat, _, t, rt, p, rp in reversed(rows):
        print(f"{lat:7.1f}{mean(t):8.2f}{mean(rt):8.2f}{mean(t) - mean(rt):+7.2f}"
              f"{mean(p):7.2f}{mean(rp):7.2f}{max(t) - min(t):7.1f}{max(rt) - min(rt):7.1f}")
    w = sum(s for _, s, *_ in rows)
    rms_t = math.sqrt(sum(s * sum((a - b) ** 2 for a, b in zip(t, rt)) / 12 for _, s, t, rt, _, _ in rows) / w)
    rms_p = math.sqrt(sum(s * sum((a - b) ** 2 for a, b in zip(p, rp)) / 12 for _, s, _, _, p, rp in rows) / w)
    print(f"Area-weighted RMS over latitude and month: T {rms_t:.2f} K, rain {rms_p:.2f} mm/day\n")

    # Sea ice by month.
    north = [sum(r["ice_cover"] * r["area_share"] for r in monthly if r["month"] == m + 1
                 and r["latitude_deg"] > 0) * EARTH_AREA for m in range(12)]
    south = [sum(r["ice_cover"] * r["area_share"] for r in monthly if r["month"] == m + 1
                 and r["latitude_deg"] < 0) * EARTH_AREA for m in range(12)]
    ref_n = at_model_months([r["north_area"] for r in ice_ref])
    ref_s = at_model_months([r["south_area"] for r in ice_ref])
    print("Sea-ice area (million km²; reference is NSIDC area, not extent). Model months from\n"
          "the March equinox.")
    print(f"{'month':>10}{'north':>8}{'ref':>7}{'south':>8}{'ref':>7}")
    for m in range(12):
        print(f"{month_label(m):>10}{north[m]:8.1f}{ref_n[m]:7.1f}{south[m]:8.1f}{ref_s[m]:7.1f}")
    print(f"{'annual':>10}{mean(north):8.1f}{mean(ref_n):7.1f}{mean(south):8.1f}{mean(ref_s):7.1f}")
    print(f"{'max':>10}{max(north):8.1f}{max(ref_n):7.1f}{max(south):8.1f}{max(ref_s):7.1f}")
    print(f"{'min':>10}{min(north):8.1f}{min(ref_n):7.1f}{min(south):8.1f}{min(ref_s):7.1f}")


if __name__ == "__main__":
    main()
