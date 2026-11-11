//! End-to-end: build a small world, step it, check the insolation fields.

use fms_sim::{World, WorldConfig};

fn small_world() -> World {
    World::new(WorldConfig {
        terrain_frequency: 20,
        climate_frequency: 8,
        ..WorldConfig::default()
    })
}

#[test]
fn builds_registers_core_fields_and_steps() {
    let mut world = small_world();
    assert!(world.adjustments.is_empty());
    assert_eq!(world.fields.len(), 13);
    for name in [
        "geo.terrain_latitude",
        "geo.terrain_cell_area",
        "geo.climate_latitude",
        "geo.climate_cell_area",
        "planet.toa_insolation",
        "planet.daily_insolation",
        "terrain.elevation",
        "terrain.slope",
        "terrain.sea",
        "terrain.land_fraction",
        "terrain.surface_height",
        "terrain.surface_height_stddev",
        "terrain.mean_slope",
    ] {
        assert!(world.fields.key(name).is_ok(), "{name}");
    }
    let report = world.step();
    assert_eq!(report.date.step, 0);
    assert_eq!(world.time(), world.scheduler.calendar().dt());
}

#[test]
fn cell_areas_cover_the_planet() {
    let world = small_world();
    let r = world.planet.params().radius;
    let sphere = 4.0 * std::f64::consts::PI * r * r;
    for id in [world.core.terrain_cell_area, world.core.climate_cell_area] {
        let total: f64 = world.fields.read(id).iter().sum();
        assert!((total / sphere - 1.0).abs() < 1e-10);
    }
}

#[test]
fn global_mean_insolation_is_a_quarter_of_the_flux() {
    let mut world = small_world();
    world.run_steps(17); // arbitrary time of day
    let areas = world.fields.read(world.core.climate_cell_area).to_vec();
    let total: f64 = areas.iter().sum();
    let t = 16.0 * world.scheduler.calendar().dt(); // values are from the last step
    let r = world.planet.orbit_state(t).distance_ratio;
    let expected = world.planet.mean_flux() / 4.0 / (r * r);
    for id in [world.core.toa_insolation, world.core.daily_insolation] {
        let mean: f64 = world
            .fields
            .read(id)
            .iter()
            .zip(&areas)
            .map(|(v, a)| v * a)
            .sum::<f64>()
            / total;
        // Coarse grid (642 cells) quadrature error on a non-smooth terminator.
        assert!((mean / expected - 1.0).abs() < 0.01, "{mean} vs {expected}");
    }
}

#[test]
fn identical_runs_are_bit_identical() {
    let run = || {
        let mut w = small_world();
        w.run_steps(100);
        w.fields.read(w.core.toa_insolation).to_vec()
    };
    let (a, b) = (run(), run());
    assert!(a.iter().zip(&b).all(|(x, y)| x.to_bits() == y.to_bits()));
}

/// Hash of the optimized climate grid and the insolation field after 100 steps.
/// CI runs this on macOS, Windows and Linux: a mismatch on any of them means
/// some computation is not bit-reproducible across platforms.
#[test]
fn results_are_bit_identical_across_platforms() {
    let mut w = small_world();
    w.run_steps(100);
    let hash = |values: &mut dyn Iterator<Item = f64>| {
        values.fold(0xcbf2_9ce4_8422_2325u64, |h, v| {
            (h ^ v.to_bits()).wrapping_mul(0x0100_0000_01b3)
        })
    };
    let grid = hash(
        &mut w
            .climate
            .mesh
            .cell_centers()
            .iter()
            .flat_map(|p| p.to_array()),
    );
    let insolation = hash(&mut w.fields.read(w.core.toa_insolation).iter().copied());
    let daily = hash(&mut w.fields.read(w.core.daily_insolation).iter().copied());
    assert_eq!(
        [grid, insolation, daily],
        GOLDEN,
        "results differ from the reference platform: {grid:#018x}, {insolation:#018x}, {daily:#018x}"
    );
}

const GOLDEN: [u64; 3] = [
    0x3d60_e023_2d57_f745,
    0xc52d_c4e9_50d4_e23e,
    0xba12_f7fe_2830_1500,
];

#[test]
fn terrain_starts_as_a_blank_ocean() {
    let world = small_world();
    let t = world.core.terrain;
    assert!(
        world
            .fields
            .read(t.elevation)
            .iter()
            .all(|&z| z == fms_sim::BLANK_OCEAN_ELEVATION)
    );
    for id in [t.land_fraction, t.surface_height, t.mean_slope] {
        assert!(world.fields.read(id).iter().all(|&v| v == 0.0));
    }
}

#[test]
fn setting_elevation_updates_the_derived_fields() {
    let mut world = small_world();
    // Land north of the equator, rising toward the pole; ocean to the south.
    let elevation: Vec<f64> = world
        .terrain
        .mesh
        .cell_centers()
        .iter()
        .map(|p| if p.z > 0.0 { 3000.0 * p.z } else { -2000.0 })
        .collect();
    world.set_elevation(elevation);
    let t = world.core.terrain;
    let areas = world.fields.read(world.core.climate_cell_area).to_vec();
    let land = world.fields.read(t.land_fraction);
    let total: f64 = areas.iter().sum();
    let fraction: f64 = land.iter().zip(&areas).map(|(f, a)| f * a).sum::<f64>() / total;
    assert!((fraction - 0.5).abs() < 0.03, "land fraction {fraction}");
    let slope = world.fields.read(t.slope);
    assert!(slope.iter().any(|&s| s > 0.0));
    assert!(world.fields.first_invalid(t.surface_height.key()).is_none());
}

#[test]
fn loads_a_png_heightmap() {
    use std::io::BufWriter;
    let path = std::env::temp_dir().join(format!("fms-sim-{}-heights.png", std::process::id()));
    // Left half white (land at 1000 m), right half black (sea floor at -1000 m).
    let (w, h) = (64u32, 32u32);
    let bytes: Vec<u8> = (0..w * h)
        .map(|i| if i % w < w / 2 { 255 } else { 0 })
        .collect();
    let mut enc = png::Encoder::new(BufWriter::new(std::fs::File::create(&path).unwrap()), w, h);
    enc.set_color(png::ColorType::Grayscale);
    enc.set_depth(png::BitDepth::Eight);
    enc.write_header()
        .unwrap()
        .write_image_data(&bytes)
        .unwrap();

    let mut world = small_world();
    let import = world
        .load_heightmap(&path, Some((-1000.0, 1000.0)))
        .unwrap();
    assert_eq!(import.size, (64, 32));
    let elevation = world.fields.read(world.core.terrain.elevation).to_vec();
    for (p, &z) in world
        .terrain
        .mesh
        .cell_centers()
        .iter()
        .zip(elevation.iter())
    {
        let lon = fms_grid::coords::to_lat_lon(*p).1.to_degrees();
        // Cells near 0° and ±180° straddle the boundary and mix both halves.
        if lon < -10.0 && lon > -170.0 {
            assert!((z - 1000.0).abs() < 1e-3, "west {lon}: {z}");
        } else if lon > 10.0 && lon < 170.0 {
            assert!((z + 1000.0).abs() < 1e-3, "east {lon}: {z}");
        }
    }
    assert!(world.load_heightmap(&path, None).is_err());
}
