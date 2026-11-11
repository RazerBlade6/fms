"""Build the committed Earth reference tables in reference/earth from the raw
downloads in data/reference/raw (scripts/fetch-reference.sh). From the
repository root, with the project virtualenv (numpy, netCDF4):
    .venv/bin/python scripts/build_reference.py
Sources and caveats: reference/earth/SOURCES.md. Every table covers 1991–2020
and calendar months."""
import csv
import json

import netCDF4
import numpy as np

RAW = "data/reference/raw"
OUT = "reference/earth"
YEARS = range(1991, 2021)
# NSIDC's Arctic pole hole (million km²; user guide v4, Table 8): left out of
# the northern area, assumed fully ice-covered here, as extent assumes it.
POLE_HOLE = {y: (0.31 if y <= 2007 else 0.029) for y in YEARS}


def field(name, var):
    d = netCDF4.Dataset(f"{RAW}/{name}.nc")
    lat = np.asarray(d["lat"][:], dtype=float)
    lon = np.asarray(d["lon"][:], dtype=float)
    data = np.ma.filled(np.ma.asarray(d[var][:], dtype=float), np.nan)
    return lat, lon, data


def band_shares(lat):
    """Area share of each grid row, with row edges midway between rows."""
    order = np.argsort(lat)
    s = np.sort(lat)
    edges = np.concatenate([[-90.0], (s[1:] + s[:-1]) / 2, [90.0]])
    share_sorted = (np.sin(np.radians(edges[1:])) - np.sin(np.radians(edges[:-1]))) / 2
    share = np.empty_like(share_sorted)
    share[order] = share_sorted
    return share


def wrap(lon):
    return np.where(lon > 180, lon - 360, lon)


def region_mask(lat, lon, land, box):
    la, lo = np.meshgrid(lat, wrap(lon), indexing="ij")
    m = (la >= box["lat"][0]) & (la <= box["lat"][1]) & (lo >= box["lon"][0]) & (lo <= box["lon"][1])
    if box["surface"] == "land":
        m &= land > 0.5
    elif box["surface"] == "ocean":
        m &= land <= 0.5
    return m


def weighted(values, weights, mask):
    w = np.where(mask & np.isfinite(values), weights, 0.0)
    return float(np.nansum(np.where(w > 0, values, 0.0) * w) / w.sum()) if w.sum() > 0 else float("nan")


def main():
    regions = {k: v for k, v in json.load(open(f"{OUT}/regions.json")).items() if not k.startswith("_")}

    lat, lon, skt = field("ncep_skt", "skt")
    skt = skt + 273.15  # degC to K
    _, _, air = field("ncep_air2m", "air")
    _, _, land = field("ncep_land", "land")
    _, _, hgt = field("ncep_hgt", "hgt")
    land, hgt = land[0], hgt[0]
    glat, glon, pr = field("gpcp_precip", "precip")

    n_share = band_shares(lat)
    g_share = band_shares(glat)
    n_w = np.repeat(n_share[:, None] / len(lon), len(lon), axis=1)
    g_w = np.repeat(g_share[:, None] / len(glon), len(glon), axis=1)
    # NCEP's land mask on the GPCP grid, nearest point.
    iy = np.abs(lat[None, :] - glat[:, None]).argmin(axis=1)
    ix = np.abs(((lon[None, :] - glon[:, None]) + 180) % 360 - 180).argmin(axis=1)
    g_land = land[np.ix_(iy, ix)]

    with open(f"{OUT}/zonal_ncep.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["month", "latitude_deg", "area_share", "skin_temperature_k", "air_2m_k"])
        for m in range(12):
            for j in np.argsort(lat):
                w.writerow([m + 1, f"{lat[j]:.3f}", f"{n_share[j]:.6f}",
                            f"{np.nanmean(skt[m, j]):.3f}", f"{np.nanmean(air[m, j]):.3f}"])

    with open(f"{OUT}/zonal_gpcp.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["month", "latitude_deg", "area_share", "precipitation_mm_day"])
        for m in range(12):
            for j in np.argsort(glat):
                w.writerow([m + 1, f"{glat[j]:.3f}", f"{g_share[j]:.6f}", f"{np.nanmean(pr[m, j]):.4f}"])

    with open(f"{OUT}/regions.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["region", "skin_temperature_k", "air_2m_k", "precipitation_mm_day",
                    "ncep_height_m", "area_share"])
        for name, box in regions.items():
            nm = region_mask(lat, lon, land, box)
            gm = region_mask(glat, glon, g_land, box)
            w.writerow([name,
                        f"{weighted(skt.mean(axis=0), n_w, nm):.2f}",
                        f"{weighted(air.mean(axis=0), n_w, nm):.2f}",
                        f"{weighted(pr.mean(axis=0), g_w, gm):.3f}",
                        f"{weighted(hgt, n_w, nm):.0f}",
                        f"{float(n_w[nm].sum()):.5f}"])

    ice = {}
    for h in "NS":
        for m in range(1, 13):
            rows = list(csv.reader(open(f"{RAW}/nsidc_{h}_{m:02d}.csv")))[1:]
            rows = [[c.strip() for c in r] for r in rows if r and int(r[0]) in YEARS]
            extent = [float(r[4]) for r in rows]
            area = [float(r[5]) + (POLE_HOLE[int(r[0])] if h == "N" else 0.0) for r in rows]
            if len(rows) != len(YEARS) or min(extent + area) < 0:
                raise SystemExit(f"NSIDC {h} month {m}: missing years in 1991–2020")
            ice[h, m] = (sum(area) / len(area), sum(extent) / len(extent))
    with open(f"{OUT}/sea_ice.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["month", "north_area", "north_extent", "south_area", "south_extent"])
        for m in range(1, 13):
            w.writerow([m] + [f"{v:.2f}" for h in "NS" for v in ice[h, m]])
    print("wrote", ", ".join(f"{OUT}/{n}" for n in ["zonal_ncep.csv", "zonal_gpcp.csv", "regions.csv", "sea_ice.csv"]))


if __name__ == "__main__":
    main()
