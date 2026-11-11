"""Compare an Earth climate run against the reference in reference/earth.
From the repository root (plain Python, no packages needed):
    python3 scripts/compare_earth.py output/<run>
where <run> is the --out-dir of `fms-climate climate` (or its climate/ folder).

The model's final year is compared with observations: ERA5 temperature,
terrain height and wind, GPCP rain and NSIDC sea ice for 1991–2020, and CERES
radiation and cloud for 2001–2020. Model months start at the March equinox, so
the reference is interpolated to the model's mid-months. Model cells are
weighted equally (the icosahedral cells are close to equal area). Sources and
caveats: reference/earth/SOURCES.md."""
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


def rms(rows, i, j):
    """Area-weighted RMS over bands and months between columns i and j."""
    w = sum(r[1] for r in rows)
    return math.sqrt(sum(r[1] * sum((a - b) ** 2 for a, b in zip(r[i], r[j])) / 12 for r in rows) / w)


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
    zonal_annual = numeric(read(os.path.join(run, "zonal_mean.csv")))
    boxes = {k: v for k, v in json.load(open(os.path.join(REF, "regions.json"))).items()
             if not k.startswith("_")}
    regions = {r["region"]: r for r in read(os.path.join(REF, "regions.csv"))}
    era5 = numeric(read(os.path.join(REF, "zonal_era5.csv")))
    gpcp = numeric(read(os.path.join(REF, "zonal_gpcp.csv")))
    ceres = numeric(read(os.path.join(REF, "zonal_ceres.csv")))
    ice_ref = numeric(read(os.path.join(REF, "sea_ice.csv")))

    print(f"{run}: final year {last:.0f}, {len(cells)} cells\n")
    albedo = float("nan")
    summary = os.path.join(run, "summary.txt")
    if os.path.exists(summary):
        found = re.search(r"planetary albedo ([\d.]+)", open(summary, encoding="utf-8").read())
        albedo = float(found.group(1)) if found else albedo
    g = regions["Global"]
    print(f"Planetary albedo {fmt(albedo, '.3f')} (CERES {float(g['albedo']):.3f})\n")

    selected = {}
    for name, box in boxes.items():
        selected[name] = [c for c in cells
                          if in_box(c["latitude_deg"], c["longitude_deg"], c["land_fraction"], box)]

    def region_mean(name, key):
        sel = selected[name]
        return mean([c[key] for c in sel]) if sel and key in sel[0] else float("nan")

    # Regions: temperature and rain.
    print("Regions, annual mean. T is skin temperature, Ta near-surface air (model: bottom layer\n"
          "brought to the surface; ERA5: 2 m). 'at ERA5 z' moves the model's T to ERA5's terrain\n"
          f"height at {LAPSE * 1000:.1f} K/km. NCEP is a second reanalysis's T, for scale.")
    print(f"{'region':16}{'T':>8}{'ERA5':>8}{'ΔT':>7}{'NCEP':>8}{'Ta':>8}{'ERA5':>8}{'ΔTa':>7}"
          f"{'rain':>7}{'GPCP':>7}{'ratio':>7}{'z':>7}{'ERA5 z':>7}{'ΔT at ERA5 z':>13}")
    for name in boxes:
        ref = regions.get(name)
        if not selected[name] or ref is None:
            print(f"{name:16} no cells")
            continue
        t, ta = region_mean(name, "rad.surface_temperature"), region_mean(name, "atm.surface_air_temperature")
        rain, z = region_mean(name, "atm.precipitation"), region_mean(name, "surface_height_m")
        rt, rta, rr, rz, nt = (float(ref[k]) for k in ("skin_temperature_k", "air_2m_k",
                                                       "precipitation_mm_day", "height_m",
                                                       "ncep_skin_temperature_k"))
        print(f"{name:16}{t:8.2f}{rt:8.2f}{t - rt:+7.2f}{nt:8.2f}{fmt(ta, '8.2f')}{rta:8.2f}"
              f"{fmt(ta - rta, '+7.2f')}{rain:7.2f}{rr:7.2f}{rain / rr:7.2f}{fmt(z, '7.0f')}{rz:7.0f}"
              f"{fmt(t - rt - LAPSE * (rz - z), '+13.2f')}")
    print()

    # Regions: radiation and cloud.
    print("Regions, annual mean, against CERES (W/m²). ASR absorbed sunlight, OLR outgoing\n"
          "longwave, CRE cloud radiative effect; cloud is cover (0–1).")
    print(f"{'region':16}{'ASR':>7}{'ref':>7}{'OLR':>7}{'ref':>7}{'SW CRE':>8}{'ref':>7}"
          f"{'LW CRE':>8}{'ref':>7}{'cloud':>7}{'ref':>7}")
    for name in boxes:
        ref = regions.get(name)
        if not selected[name] or ref is None:
            continue
        vals = [(region_mean(name, m), float(ref[r])) for m, r in
                (("rad.absorbed_shortwave", "absorbed_shortwave"), ("rad.outgoing_longwave", "outgoing_longwave"),
                 ("rad.shortwave_cloud_effect", "shortwave_cloud_effect"),
                 ("rad.longwave_cloud_effect", "longwave_cloud_effect"))]
        cloud = (region_mean(name, "atm.cloud_cover"), float(ref["cloud_fraction"]))
        print(f"{name:16}" + "".join(f"{a:{7 if i < 2 else 8}.1f}{b:7.1f}" for i, (a, b) in enumerate(vals))
              + f"{cloud[0]:7.2f}{cloud[1]:7.2f}")
    print()

    # Zonal means.
    lats = sorted({r["latitude_deg"] for r in monthly})
    wind = {round(r["latitude_deg"], 2): r["wind_lower_m_s"] for r in zonal_annual}
    rows = []
    for lat, lo, hi in bands(lats):
        sel = sorted((r for r in monthly if r["latitude_deg"] == lat), key=lambda r: r["month"])
        rows.append((
            lat, sel[0]["area_share"],
            [r["rad.surface_temperature"] for r in sel], at_model_months(ref_zonal(era5, "skin_temperature_k", lo, hi)),
            [r["atm.precipitation"] for r in sel], at_model_months(ref_zonal(gpcp, "precipitation_mm_day", lo, hi)),
            [r["rad.absorbed_shortwave"] for r in sel], at_model_months(ref_zonal(ceres, "absorbed_shortwave", lo, hi)),
            [r["rad.outgoing_longwave"] for r in sel], at_model_months(ref_zonal(ceres, "outgoing_longwave", lo, hi)),
            [r["atm.cloud_cover"] for r in sel], at_model_months(ref_zonal(ceres, "cloud_fraction", lo, hi)),
            wind.get(round(lat, 2), float("nan")), mean(ref_zonal(era5, "wind_10m_east_m_s", lo, hi)),
        ))
    print("Zonal means, annual: skin temperature (K) against ERA5 and its seasonal range (warmest\n"
          "minus coldest month), rain (mm/day) against GPCP, absorbed sunlight and outgoing longwave\n"
          "(W/m²) and cloud cover against CERES, and the eastward wind (m/s): the model's bottom\n"
          "layer (~1 km deep) against ERA5 at 10 m, which is weaker.")
    print(f"{'lat':>6}{'T':>7}{'ref':>7}{'ΔT':>6}{'range':>6}{'ref':>5}{'rain':>6}{'ref':>5}"
          f"{'ASR':>6}{'ref':>5}{'OLR':>6}{'ref':>5}{'cloud':>6}{'ref':>5}{'u':>6}{'u10':>5}")
    for r in reversed(rows):
        lat, _, t, rt, p, rp, a, ra, o, ro, c, rc, u, ru = r
        print(f"{lat:6.1f}{mean(t):7.1f}{mean(rt):7.1f}{mean(t) - mean(rt):+6.1f}{max(t) - min(t):6.1f}"
              f"{max(rt) - min(rt):5.1f}{mean(p):6.2f}{mean(rp):5.1f}{mean(a):6.0f}{mean(ra):5.0f}"
              f"{mean(o):6.0f}{mean(ro):5.0f}{mean(c):6.2f}{mean(rc):5.2f}{fmt(u, '6.1f')}{ru:5.1f}")
    print(f"Area-weighted RMS over latitude and month: T {rms(rows, 2, 3):.2f} K, rain "
          f"{rms(rows, 4, 5):.2f} mm/day, ASR {rms(rows, 6, 7):.1f} W/m², OLR {rms(rows, 8, 9):.1f} W/m², "
          f"cloud {rms(rows, 10, 11):.2f}\n")

    # Sea ice by month.
    north = [sum(r["ice_cover"] * r["area_share"] for r in monthly if r["month"] == m + 1
                 and r["latitude_deg"] > 0) * EARTH_AREA for m in range(12)]
    south = [sum(r["ice_cover"] * r["area_share"] for r in monthly if r["month"] == m + 1
                 and r["latitude_deg"] < 0) * EARTH_AREA for m in range(12)]
    ref_n = at_model_months([r["north_area"] for r in ice_ref])
    ref_s = at_model_months([r["south_area"] for r in ice_ref])
    print("Sea-ice area (million km²) against NSIDC area (not extent). Model months from the\n"
          "March equinox.")
    print(f"{'month':>10}{'north':>8}{'ref':>7}{'south':>8}{'ref':>7}")
    for m in range(12):
        print(f"{month_label(m):>10}{north[m]:8.1f}{ref_n[m]:7.1f}{south[m]:8.1f}{ref_s[m]:7.1f}")
    print(f"{'annual':>10}{mean(north):8.1f}{mean(ref_n):7.1f}{mean(south):8.1f}{mean(ref_s):7.1f}")
    print(f"{'max':>10}{max(north):8.1f}{max(ref_n):7.1f}{max(south):8.1f}{max(ref_s):7.1f}")
    print(f"{'min':>10}{min(north):8.1f}{min(ref_n):7.1f}{min(south):8.1f}{min(ref_s):7.1f}")


if __name__ == "__main__":
    main()
