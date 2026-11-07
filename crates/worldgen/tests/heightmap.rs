//! Heightmap import: file readers, resampling onto a mesh, and sub-grid terrain
//! statistics, checked against synthetic heightmaps with known answers.

use std::f64::consts::PI;
use std::fs::File;
use std::io::BufWriter;
use std::path::{Path, PathBuf};

use fms_geo::coords::{from_lat_lon, to_lat_lon};
use fms_geo::spherical::polygon_moment;
use fms_geo::{LloydOptions, MeshGeometry, SphereMesh, TwoLevelMap};
use fms_worldgen::subgrid::{slope, surface_height};
use fms_worldgen::{HeightmapError, LoadOptions, Raster, TerrainStats, load, resample};
use glam::DVec3;
use tiff::encoder::{TiffEncoder, colortype};
use tiff::tags::Tag;

const RADIUS: f64 = 6.371e6;

fn mesh(frequency: u32) -> SphereMesh {
    SphereMesh::centroidal(frequency, LloydOptions::default()).0
}

/// A linear function of position, so its exact mean over a cell is known.
const A: DVec3 = DVec3::new(1000.0, -500.0, 2000.0);

fn linear(lat: f64, lon: f64) -> f64 {
    A.dot(from_lat_lon(lat, lon))
}

/// Exact area-weighted mean of `A · p` over each cell: `A · ∫p dA / area`.
fn exact_linear_means(mesh: &SphereMesh) -> Vec<f64> {
    let areas = mesh.cell_areas();
    (0..mesh.n_cells() as u32)
        .map(|c| {
            let ring = mesh
                .cell_vertices(c)
                .iter()
                .map(|&v| mesh.vertices()[v as usize]);
            A.dot(polygon_moment(ring)) / areas[c as usize]
        })
        .collect()
}

fn max_error(a: &[f64], b: &[f64]) -> f64 {
    a.iter()
        .zip(b)
        .map(|(x, y)| (x - y).abs())
        .fold(0.0, f64::max)
}

fn scratch(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join("fms-worldgen-tests");
    std::fs::create_dir_all(&dir).unwrap();
    dir.join(format!("{}-{name}", std::process::id()))
}

// ---------------------------------------------------------------------------
// Resampling
// ---------------------------------------------------------------------------

#[test]
fn constant_heightmap_resamples_exactly() {
    let m = mesh(16);
    let raster = Raster::from_fn(720, 360, |_, _| -1234.5);
    for z in resample(&raster, &m) {
        assert!((z + 1234.5).abs() < 1e-9, "{z}");
    }
}

#[test]
fn cell_means_match_the_exact_mean_and_converge() {
    let m = mesh(16);
    let exact = exact_linear_means(&m);
    let errors: Vec<f64> = [(360, 180), (1440, 720), (5760, 2880)]
        .iter()
        .map(|&(w, h)| max_error(&resample(&Raster::from_fn(w, h, linear), &m), &exact))
        .collect();
    // The field varies by ~|A| ≈ 2300 m per radian; a cell is ~0.08 rad across.
    // Error comes from which pixels fall inside each cell, so it shrinks with
    // the pixel size. Measured: 12.9, 2.8, 0.57 m.
    eprintln!("max error vs exact cell mean: {errors:?} m");
    assert!(errors[2] < 1.5, "{errors:?}");
    assert!(errors[1] < 0.5 * errors[0] && errors[2] < 0.5 * errors[1]);
}

#[test]
fn poles_and_the_antimeridian_have_no_seams() {
    // Cells containing a pole, or straddling ±180°, get pixels from both sides.
    let m = mesh(16);
    let exact = exact_linear_means(&m);
    let got = resample(&Raster::from_fn(1440, 720, linear), &m);
    for (c, &p) in m.cell_centers().iter().enumerate() {
        let (lat, lon) = to_lat_lon(p);
        if lat.abs() > 1.4 || lon.abs() > PI - 0.1 {
            assert!((got[c] - exact[c]).abs() < 3.0, "cell {c} at {lat},{lon}");
        }
    }
}

#[test]
fn coarse_heightmaps_fall_back_to_bilinear_samples() {
    // 16×8 pixels on 5762 cells: most cells contain no pixel centre.
    let m = mesh(24);
    let raster = Raster::from_fn(16, 8, linear);
    let got = resample(&raster, &m);
    let (lo, hi) = raster
        .data()
        .iter()
        .fold((f64::MAX, f64::MIN), |(a, b), &v| {
            (a.min(v as f64), b.max(v as f64))
        });
    let mut sampled = 0;
    for (c, &p) in m.cell_centers().iter().enumerate() {
        assert!(got[c].is_finite() && got[c] >= lo - 1e-9 && got[c] <= hi + 1e-9);
        if (got[c] - raster.sample(p)).abs() < 1e-9 {
            sampled += 1;
        }
    }
    assert!(sampled > m.n_cells() / 2, "only {sampled} cells sampled");
}

// ---------------------------------------------------------------------------
// GeoTIFF
// ---------------------------------------------------------------------------

/// Write a single-band float GeoTIFF with lat/lon georeferencing.
fn write_geotiff(path: &Path, w: u32, h: u32, data: &[f32], geo: Option<(f64, f64, bool)>) {
    let mut tiff = TiffEncoder::new(BufWriter::new(File::create(path).unwrap())).unwrap();
    let mut image = tiff.new_image::<colortype::Gray32Float>(w, h).unwrap();
    if let Some((spacing_deg, origin_lon, pixel_is_point)) = geo {
        let enc = image.encoder();
        enc.write_tag(
            Tag::ModelPixelScaleTag,
            &[spacing_deg, spacing_deg, 0.0][..],
        )
        .unwrap();
        // Pixel (0, 0)'s corner (is-area) or centre (is-point) is at 90° N.
        enc.write_tag(
            Tag::ModelTiepointTag,
            &[0.0, 0.0, 0.0, origin_lon, 90.0, 0.0][..],
        )
        .unwrap();
        let raster_type = if pixel_is_point { 2 } else { 1 };
        let keys: [u16; 12] = [1, 1, 0, 2, 1024, 0, 1, 2, 1025, 0, 1, raster_type];
        enc.write_tag(Tag::GeoKeyDirectoryTag, &keys[..]).unwrap();
        enc.write_tag(Tag::GdalNodata, "-99999").unwrap();
    }
    image.write_data(data).unwrap();
}

fn assert_same_raster(a: &Raster, b: &Raster) {
    assert_eq!((a.width(), a.height()), (b.width(), b.height()));
    for row in [0, a.height() / 2, a.height() - 1] {
        assert!((a.lat(row) - b.lat(row)).abs() < 1e-12);
    }
    for col in [0, a.width() / 2, a.width() - 1] {
        let d = (a.lon(col) - b.lon(col)).rem_euclid(2.0 * PI);
        assert!(d.min(2.0 * PI - d) < 1e-12);
    }
    assert!(
        max_error(
            &a.data().iter().map(|&v| v as f64).collect::<Vec<_>>(),
            &b.data().iter().map(|&v| v as f64).collect::<Vec<_>>(),
        ) < 1e-3
    );
}

#[test]
fn geotiff_pixel_is_area_round_trips() {
    let expected = Raster::from_fn(72, 36, linear);
    let path = scratch("area.tif");
    write_geotiff(&path, 72, 36, expected.data(), Some((5.0, -180.0, false)));
    assert_same_raster(&load(&path, &LoadOptions::default()).unwrap(), &expected);
}

#[test]
fn plain_tiff_is_read_as_a_global_image() {
    let expected = Raster::from_fn(40, 20, linear);
    let path = scratch("plain.tiff");
    write_geotiff(&path, 40, 20, expected.data(), None);
    assert_same_raster(&load(&path, &LoadOptions::default()).unwrap(), &expected);
}

#[test]
fn geotiff_grid_registered_drops_the_duplicate_column() {
    // ETOPO1 style: 73×37 samples at the grid nodes, with ±180° both present
    // and rows exactly at the poles.
    let (w, h) = (73usize, 37usize);
    let mut data = Vec::with_capacity(w * h);
    for row in 0..h {
        for col in 0..w {
            let lat = (90.0 - 5.0 * row as f64).to_radians();
            let lon = (-180.0 + 5.0 * col as f64).to_radians();
            data.push(linear(lat, lon) as f32);
        }
    }
    let path = scratch("point.tif");
    write_geotiff(&path, w as u32, h as u32, &data, Some((5.0, -180.0, true)));
    let r = load(&path, &LoadOptions::default()).unwrap();
    assert_eq!((r.width(), r.height()), (72, 37));
    assert!((r.lat(0) - PI / 2.0).abs() < 1e-12 && (r.lon(0) + PI).abs() < 1e-12);
    assert!((r.get(18, 36) as f64 - linear(0.0, 0.0)).abs() < 1e-3);
}

#[test]
fn geotiff_decimates_while_reading() {
    let full = Raster::from_fn(360, 180, linear);
    let path = scratch("decimate.tif");
    write_geotiff(&path, 360, 180, full.data(), Some((1.0, -180.0, false)));
    let options = LoadOptions {
        target_spacing: Some(4.5f64.to_radians()),
        ..LoadOptions::default()
    };
    let r = load(&path, &options).unwrap();
    // Largest block ≤ 4 that divides 360 and 180.
    assert_eq!((r.width(), r.height()), (90, 45));
    // A block mean of a linear field is close to its value at the block centre.
    // Near the poles the area-weighted mean sits measurably equatorward of the
    // centre (the unit tests check the exact block mean), so compare below 60°.
    let mut worst = 0.0f64;
    for row in (0..r.height()).filter(|&row| r.lat(row).abs() < 60f64.to_radians()) {
        for col in 0..r.width() {
            worst = worst.max((r.get(row, col) - linear(r.lat(row), r.lon(col))).abs());
        }
    }
    assert!(
        worst < 5.0,
        "block mean differs from the centre value by {worst} m"
    );
}

#[test]
fn geotiff_nodata_is_rejected() {
    let mut data = Raster::from_fn(72, 36, linear).data().to_vec();
    data[100] = -99999.0;
    let path = scratch("nodata.tif");
    write_geotiff(&path, 72, 36, &data, Some((5.0, -180.0, false)));
    let err = load(&path, &LoadOptions::default()).unwrap_err();
    assert!(err.to_string().contains("1 pixels have no data"), "{err}");
}

#[test]
fn regional_geotiff_is_rejected() {
    let path = scratch("regional.tif");
    write_geotiff(
        &path,
        72,
        36,
        &vec![0.0; 72 * 36],
        Some((1.0, -10.0, false)),
    );
    let err = load(&path, &LoadOptions::default()).unwrap_err();
    assert!(matches!(err, HeightmapError::Raster(..)), "{err}");
}

// ---------------------------------------------------------------------------
// PNG
// ---------------------------------------------------------------------------

fn write_png(
    path: &Path,
    w: u32,
    h: u32,
    color: png::ColorType,
    depth: png::BitDepth,
    bytes: &[u8],
) {
    let mut enc = png::Encoder::new(BufWriter::new(File::create(path).unwrap()), w, h);
    enc.set_color(color);
    enc.set_depth(depth);
    enc.write_header().unwrap().write_image_data(bytes).unwrap();
}

#[test]
fn png_16_bit_maps_grey_levels_to_the_range() {
    let (w, h) = (64u32, 32u32);
    let mut bytes = Vec::new();
    for row in 0..h {
        for col in 0..w {
            let v = ((row * w + col) * 32) as u16; // 0 … 65504
            bytes.extend_from_slice(&v.to_be_bytes());
        }
    }
    let path = scratch("height16.png");
    write_png(
        &path,
        w,
        h,
        png::ColorType::Grayscale,
        png::BitDepth::Sixteen,
        &bytes,
    );
    let options = LoadOptions {
        png_range: Some((-8000.0, 6000.0)),
        ..LoadOptions::default()
    };
    let r = load(&path, &options).unwrap();
    assert_eq!((r.width(), r.height()), (64, 32));
    assert!((r.get(0, 0) + 8000.0).abs() < 1e-3);
    let last = -8000.0 + 14000.0 * (2047.0 * 32.0) / 65535.0;
    assert!((r.get(31, 63) - last).abs() < 1e-2);
}

#[test]
fn png_8_bit_and_rejections() {
    let path = scratch("height8.png");
    write_png(
        &path,
        8,
        4,
        png::ColorType::Grayscale,
        png::BitDepth::Eight,
        &[255; 32],
    );
    let options = LoadOptions {
        png_range: Some((0.0, 1000.0)),
        ..LoadOptions::default()
    };
    assert!((load(&path, &options).unwrap().get(2, 3) - 1000.0).abs() < 1e-3);

    let err = load(&path, &LoadOptions::default()).unwrap_err();
    assert!(err.to_string().contains("black and white"), "{err}");

    let rgb = scratch("colour.png");
    write_png(
        &rgb,
        8,
        4,
        png::ColorType::Rgb,
        png::BitDepth::Eight,
        &[0; 96],
    );
    let err = load(&rgb, &options).unwrap_err();
    assert!(err.to_string().contains("greyscale"), "{err}");

    let err = load(Path::new("heights.jpg"), &options).unwrap_err();
    assert!(err.to_string().contains("unknown format"), "{err}");
}

// ---------------------------------------------------------------------------
// Sub-grid statistics
// ---------------------------------------------------------------------------

struct Grids {
    fine: SphereMesh,
    geometry: MeshGeometry,
    map: TwoLevelMap,
}

fn grids() -> Grids {
    let coarse = mesh(8);
    let fine = mesh(48);
    let geometry = MeshGeometry::new(&fine);
    let map = TwoLevelMap::new(&coarse, &fine, &geometry.cell_areas);
    Grids {
        fine,
        geometry,
        map,
    }
}

fn stats(g: &Grids, elevation: &[f64]) -> TerrainStats {
    let s = slope(&surface_height(elevation), &g.fine, &g.geometry, RADIUS);
    TerrainStats::compute(elevation, &s, &g.map)
}

#[test]
fn flat_land_and_flat_ocean() {
    let g = grids();
    let land = stats(&g, &vec![500.0; g.fine.n_cells()]);
    for c in 0..land.land_fraction.len() {
        assert!((land.land_fraction[c] - 1.0).abs() < 1e-12);
        assert!((land.surface_height[c] - 500.0).abs() < 1e-9);
        assert!(land.surface_height_stddev[c] < 1e-6);
        assert!(land.mean_slope[c] < 1e-12);
    }
    let ocean = stats(&g, &vec![-3000.0; g.fine.n_cells()]);
    assert!(ocean.land_fraction.iter().all(|&f| f == 0.0));
    assert!(ocean.surface_height.iter().all(|&h| h == 0.0));
}

#[test]
fn slope_of_a_tilted_surface() {
    // z = 6000 + H·sin(lat), all above sea level: gradient magnitude H·cos(lat) / r.
    let g = grids();
    let h = 5000.0;
    let elevation: Vec<f64> = g
        .fine
        .cell_centers()
        .iter()
        .map(|p| 6000.0 + h * p.z)
        .collect();
    let s = slope(&surface_height(&elevation), &g.fine, &g.geometry, RADIUS);
    let rel: Vec<f64> = g
        .fine
        .cell_centers()
        .iter()
        .zip(&s)
        .filter(|(p, _)| p.z.abs() < 0.9)
        .map(|(p, &got)| {
            let expected = h * (1.0 - p.z * p.z).sqrt() / RADIUS;
            (got - expected).abs() / expected
        })
        .collect();
    let worst = rel.iter().fold(0.0f64, |a, &b| a.max(b));
    let rms = (rel.iter().map(|e| e * e).sum::<f64>() / rel.len() as f64).sqrt();
    // Perot reconstruction errors peak at the pentagons (OBSERVATIONS.md O-1).
    eprintln!("slope relative error: max {worst:.4}, rms {rms:.5}");
    assert!(worst < 0.05 && rms < 0.005, "max {worst}, rms {rms}");
}

#[test]
fn hemispheric_continent() {
    // Northern hemisphere land with a ridge; southern ocean.
    let g = grids();
    let elevation: Vec<f64> = g
        .fine
        .cell_centers()
        .iter()
        .map(|p| {
            if p.z > 0.0 {
                200.0 + 2000.0 * p.x.max(0.0)
            } else {
                -4000.0
            }
        })
        .collect();
    let st = stats(&g, &elevation);
    let areas = mesh(8).cell_areas();
    let total: f64 = areas.iter().sum();
    let land: f64 = st
        .land_fraction
        .iter()
        .zip(&areas)
        .map(|(f, a)| f * a)
        .sum::<f64>()
        / total;
    assert!((land - 0.5).abs() < 0.01, "land fraction {land}");
    for c in 0..areas.len() {
        let f = st.land_fraction[c];
        assert!((0.0..=1.0).contains(&f));
        // Coastal cells mix sea level and land, so they have spread; the mean
        // surface height stays between the lowest and highest surface.
        assert!(st.surface_height[c] >= 0.0 && st.surface_height[c] <= 2200.0);
        if f > 1e-9 && f < 1.0 - 1e-9 {
            assert!(st.surface_height_stddev[c] > 10.0);
        }
    }
}
