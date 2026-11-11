"""Build the committed Earth reference tables in reference/earth from the raw
files in data/reference/raw: those scripts/fetch-reference.sh downloads, plus
two that need an account (reference/earth/SOURCES.md says how to get them):
era5.grib and CERES_EBAF-4.nc. From the repository root, with the project
virtualenv (numpy, netCDF4, eccodes):
    .venv/bin/python scripts/build_reference.py
Reading ERA5 takes a few minutes. All tables use calendar months; ERA5, GPCP,
NCEP and NSIDC cover 1991–2020 and CERES 2001–2020."""
import csv
import json

import eccodes
import netCDF4
import numpy as np

RAW = "data/reference/raw"
OUT = "reference/earth"
YEARS = range(1991, 2021)
CERES_YEARS = range(2001, 2021)
GRAVITY = 9.80665  # ERA5's geopotential to height
# NSIDC's Arctic pole hole (million km²; user guide v4, Table 8): left out of
# the northern area, assumed fully ice-covered here, as extent assumes it.
POLE_HOLE = {y: (0.31 if y <= 2007 else 0.029) for y in YEARS}
ERA5_MONTHLY = ("skt", "2t", "msl", "10u")
ERA5_FIXED = ("z", "lsm")


def field(name, var):
    d = netCDF4.Dataset(f"{RAW}/{name}.nc")
    lat = np.asarray(d["lat"][:], dtype=float)
    lon = np.asarray(d["lon"][:], dtype=float)
    data = np.ma.filled(np.ma.asarray(d[var][:], dtype=float), np.nan)
    return lat, lon, data


def era5():
    """ERA5 monthly climatology, 1991–2020: {name: (12, lat, lon)}, plus the
    fixed fields as (lat, lon), and the grid."""
    sums, counts, fixed = {}, {}, {}
    lat = lon = None
    with open(f"{RAW}/era5.grib", "rb") as f:
        while True:
            h = eccodes.codes_grib_new_from_file(f)
            if h is None:
                break
            try:
                name = eccodes.codes_get(h, "shortName")
                date = eccodes.codes_get(h, "validityDate")
                year, month = date // 10000, date // 100 % 100
                if year not in YEARS or name not in ERA5_MONTHLY + ERA5_FIXED:
                    continue
                if name in ERA5_FIXED and name in fixed:
                    continue
                if lat is None:
                    ni, nj = eccodes.codes_get(h, "Ni"), eccodes.codes_get(h, "Nj")
                    lat0 = eccodes.codes_get(h, "latitudeOfFirstGridPointInDegrees")
                    step = eccodes.codes_get(h, "jDirectionIncrementInDegrees")
                    lat = lat0 - step * np.arange(nj)
                    lon = eccodes.codes_get(h, "longitudeOfFirstGridPointInDegrees") \
                        + eccodes.codes_get(h, "iDirectionIncrementInDegrees") * np.arange(ni)
                values = eccodes.codes_get_values(h).reshape(len(lat), len(lon))
                if name in ERA5_FIXED:
                    fixed[name] = values
                else:
                    if name not in sums:
                        sums[name] = np.zeros((12, len(lat), len(lon)))
                        counts[name] = np.zeros(12)
                    sums[name][month - 1] += values
                    counts[name][month - 1] += 1
            finally:
                eccodes.codes_release(h)
    for name in ERA5_MONTHLY:
        if not (counts[name] == len(YEARS)).all():
            raise SystemExit(f"ERA5 {name}: expected {len(YEARS)} years per month, got {counts[name]}")
    clim = {k: sums[k] / counts[k][:, None, None] for k in ERA5_MONTHLY}
    return lat, lon, clim, fixed


def ceres():
    """CERES EBAF monthly climatology, 2001–2020: {name: (12, lat, lon)}."""
    d = netCDF4.Dataset(f"{RAW}/CERES_EBAF-4.nc")
    dates = netCDF4.num2date(d["time"][:], d["time"].units)
    pick = [[i for i, t in enumerate(dates) if t.year in CERES_YEARS and t.month == m]
            for m in range(1, 13)]
    if any(len(p) != len(CERES_YEARS) for p in pick):
        raise SystemExit("CERES: missing months in 2001–2020")

    def clim(var):
        a = np.ma.filled(np.ma.asarray(d[var][:], dtype=float), np.nan)
        return np.stack([a[p].mean(axis=0) for p in pick])

    solar, sw, lw = clim("solar_mon"), clim("toa_sw_all_mon"), clim("toa_lw_all_mon")
    sw_clr, lw_clr = clim("toa_sw_clr_c_mon"), clim("toa_lw_clr_c_mon")
    fields = {
        "solar": solar,
        "absorbed_shortwave": solar - sw,
        "outgoing_longwave": lw,
        "shortwave_cloud_effect": sw_clr - sw,
        "longwave_cloud_effect": lw_clr - lw,
        "cloud_fraction": clim("cldarea_total_daynight_mon") / 100,
        "reflected_shortwave": sw,
    }
    return np.asarray(d["lat"][:], dtype=float), np.asarray(d["lon"][:], dtype=float), fields


def band_shares(lat):
    """Area share of each grid row, with row edges midway between rows."""
    order = np.argsort(lat)
    s = np.sort(lat)
    edges = np.concatenate([[-90.0], (s[1:] + s[:-1]) / 2, [90.0]])
    share_sorted = (np.sin(np.radians(edges[1:])) - np.sin(np.radians(edges[:-1]))) / 2
    share = np.empty_like(share_sorted)
    share[order] = share_sorted
    return share


def weights(lat, lon):
    return np.repeat(band_shares(lat)[:, None] / len(lon), len(lon), axis=1)


def wrap(lon):
    return np.where(lon > 180, lon - 360, lon)


def nearest(src_lat, src_lon, src, lat, lon):
    """A source field sampled at the nearest point of another grid."""
    iy = np.abs(src_lat[None, :] - lat[:, None]).argmin(axis=1)
    ix = np.abs(((src_lon[None, :] - lon[:, None]) + 180) % 360 - 180).argmin(axis=1)
    return src[np.ix_(iy, ix)]


def region_mask(lat, lon, land, box):
    la, lo = np.meshgrid(lat, wrap(lon), indexing="ij")
    m = (la >= box["lat"][0]) & (la <= box["lat"][1]) & (lo >= box["lon"][0]) & (lo <= box["lon"][1])
    if box["surface"] == "land":
        m &= land > 0.5
    elif box["surface"] == "ocean":
        m &= land <= 0.5
    return m


def weighted(values, w, mask):
    w = np.where(mask & np.isfinite(values), w, 0.0)
    return float(np.sum(np.where(w > 0, values, 0.0) * w) / w.sum()) if w.sum() > 0 else float("nan")


def zonal_rows(lat, monthly, coarsen=None):
    """Monthly zonal means per grid row, or per `coarsen`-degree band:
    yields (month, latitude, area share, [values])."""
    share = band_shares(lat)
    if coarsen is None:
        groups = [(lat[j], [j]) for j in np.argsort(lat)]
    else:
        key = np.clip(np.floor((lat + 90) / coarsen), 0, 180 / coarsen - 1)
        groups = [(-90 + (k + 0.5) * coarsen, list(np.nonzero(key == k)[0])) for k in np.unique(key)]
    for m in range(12):
        for centre, rows in groups:
            w = share[rows]
            yield m + 1, centre, w.sum(), [float((np.nanmean(f[m, rows], axis=1) * w).sum() / w.sum())
                                           for f in monthly]


def write_zonal(path, header, rows, digits):
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["month", "latitude_deg", "area_share"] + header)
        for m, lat, share, values in rows:
            w.writerow([m, f"{lat:.3f}", f"{share:.6f}"] + [f"{v:.{d}f}" for v, d in zip(values, digits)])


def main():
    regions = {k: v for k, v in json.load(open(f"{OUT}/regions.json")).items() if not k.startswith("_")}

    print("reading ERA5 ...", flush=True)
    elat, elon, e, efixed = era5()
    eland, eheight = efixed["lsm"], efixed["z"] / GRAVITY
    clat, clon, c = ceres()
    nlat, nlon, nskt = field("ncep_skt", "skt")
    nskt = nskt + 273.15  # degC to K
    _, _, nland = field("ncep_land", "land")
    nland = nland[0]
    glat, glon, pr = field("gpcp_precip", "precip")
    g_land = nearest(elat, elon, eland, glat, glon)
    c_land = nearest(elat, elon, eland, clat, clon)

    write_zonal(f"{OUT}/zonal_era5.csv",
                ["skin_temperature_k", "air_2m_k", "sea_level_pressure_hpa", "wind_10m_east_m_s"],
                zonal_rows(elat, [e["skt"], e["2t"], e["msl"] / 100, e["10u"]], coarsen=1.0), [3, 3, 2, 3])
    write_zonal(f"{OUT}/zonal_gpcp.csv", ["precipitation_mm_day"], zonal_rows(glat, [pr]), [4])
    ckeys = ["solar", "absorbed_shortwave", "outgoing_longwave", "shortwave_cloud_effect",
             "longwave_cloud_effect", "cloud_fraction"]
    write_zonal(f"{OUT}/zonal_ceres.csv", ckeys, zonal_rows(clat, [c[k] for k in ckeys]), [2, 2, 2, 2, 2, 4])

    ew, nw, gw, cw = weights(elat, elon), weights(nlat, nlon), weights(glat, glon), weights(clat, clon)
    annual = {k: v.mean(axis=0) for k, v in c.items()}
    with open(f"{OUT}/regions.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["region", "skin_temperature_k", "air_2m_k", "height_m", "ncep_skin_temperature_k",
                    "precipitation_mm_day", "solar", "absorbed_shortwave", "outgoing_longwave",
                    "shortwave_cloud_effect", "longwave_cloud_effect", "cloud_fraction", "albedo",
                    "area_share"])
        for name, box in regions.items():
            em = region_mask(elat, elon, eland, box)
            cm = region_mask(clat, clon, c_land, box)
            rad = {k: weighted(annual[k], cw, cm) for k in annual}
            w.writerow([name,
                        f"{weighted(e['skt'].mean(axis=0), ew, em):.2f}",
                        f"{weighted(e['2t'].mean(axis=0), ew, em):.2f}",
                        f"{weighted(eheight, ew, em):.0f}",
                        f"{weighted(nskt.mean(axis=0), nw, region_mask(nlat, nlon, nland, box)):.2f}",
                        f"{weighted(pr.mean(axis=0), gw, region_mask(glat, glon, g_land, box)):.3f}"]
                       + [f"{rad[k]:.2f}" for k in ckeys[:-1]]
                       + [f"{rad['cloud_fraction']:.3f}", f"{rad['reflected_shortwave'] / rad['solar']:.4f}",
                          f"{float(ew[em].sum()):.5f}"])

    ice = {}
    for h in "NS":
        for m in range(1, 13):
            rows = list(csv.reader(open(f"{RAW}/nsidc_{h}_{m:02d}.csv")))[1:]
            rows = [[x.strip() for x in r] for r in rows if r and int(r[0]) in YEARS]
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
    print("wrote", ", ".join(f"{OUT}/{n}" for n in
                             ["zonal_era5.csv", "zonal_gpcp.csv", "zonal_ceres.csv", "regions.csv", "sea_ice.csv"]))


if __name__ == "__main__":
    main()
