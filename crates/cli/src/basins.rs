//! `fms-climate basins`: the world's sea mask and the ocean basins of the
//! wind-driven gyres, as maps and tables, without running anything.

use std::fmt::Write as _;
use std::time::Instant;

use clap::Args;
use fms_fields::{ColorScale, FieldSpec, Location};
use fms_grid::coords::to_lat_lon;
use fms_grid::{Operators, SphereMesh};
use fms_ocean::{Basin, Basins};
use fms_worldgen::SeaMask;
use glam::DVec3;

use crate::{ImageArgs, WorldArgs, build_world, ensure_dir, write_one};

#[derive(Args)]
pub struct BasinsArgs {
    #[command(flatten)]
    world: WorldArgs,
    /// Bodies below sea level to print besides the world ocean, largest
    /// first (all of them go to `sea_bodies.csv`).
    #[arg(long, default_value_t = 40)]
    bodies: usize,
    #[command(flatten)]
    image: ImageArgs,
}

/// Latitude and longitude of a unit vector, in degrees.
fn lat_lon(p: DVec3) -> (f64, f64) {
    let (lat, lon) = to_lat_lon(p);
    (lat.to_degrees(), lon.to_degrees())
}

/// Area-weighted mean position of a basin's vertices, and the latitude range
/// they span (degrees).
fn extent(mesh: &SphereMesh, areas: &[f64], basin: &Basin) -> ((f64, f64), (f64, f64)) {
    let centre = basin
        .vertices
        .iter()
        .map(|&v| areas[v as usize] * mesh.vertices()[v as usize])
        .sum::<DVec3>()
        .normalize_or_zero();
    let (low, high) = basin
        .vertices
        .iter()
        .map(|&v| lat_lon(mesh.vertices()[v as usize]).0)
        .fold((f64::MAX, f64::MIN), |(a, b), x| (a.min(x), b.max(x)));
    (lat_lon(centre), (low, high))
}

pub fn run(args: BasinsArgs) -> Result<(), String> {
    let world = build_world(&args.world)?;
    let r = world.planet.params().radius;
    let sphere_km2 = 4.0 * std::f64::consts::PI * r * r / 1e6;
    let dir = &args.image.out_dir;
    ensure_dir(dir)?;

    // The sea mask again, for its bodies and its cost.
    let elevation = world.fields.read(world.core.terrain.elevation).to_vec();
    let (mesh, areas) = (&world.terrain.mesh, &world.terrain.geometry.cell_areas);
    let mut seconds = f64::INFINITY;
    let mut mask = None;
    for _ in 0..3 {
        let t = Instant::now();
        mask = Some(SeaMask::new(&elevation, mesh, areas));
        seconds = seconds.min(t.elapsed().as_secs_f64());
    }
    let mask = mask.expect("built three times");
    println!(
        "sea mask: {} terrain cells in {:.1} ms; {} bodies at or below sea level",
        mesh.n_cells(),
        seconds * 1e3,
        mask.bodies.len()
    );
    let mut bodies: Vec<_> = mask.bodies.iter().collect();
    bodies.sort_by(|a, b| b.area_fraction.total_cmp(&a.area_fraction));
    let mut csv = String::from(
        "area_km2,area_fraction,cells,latitude_deg,longitude_deg,lowest_m,world_ocean,sea\n",
    );
    for b in &bodies {
        let (lat, lon) = lat_lon(b.centre);
        let _ = writeln!(
            csv,
            "{:.0},{:e},{},{lat:.2},{lon:.2},{:.0},{},{}",
            b.area_fraction * sphere_km2,
            b.area_fraction,
            b.cells,
            b.lowest,
            b.world_ocean,
            b.sea
        );
    }
    let path = dir.join("sea_bodies.csv");
    std::fs::write(&path, csv).map_err(|e| format!("cannot write {}: {e}", path.display()))?;
    println!(
        "  {:>12} {:>7} {:>8} {:>8} {:>7}  verdict",
        "km2", "cells", "lat", "lon", "lowest"
    );
    for b in bodies.iter().take(args.bodies + 1) {
        let (lat, lon) = lat_lon(b.centre);
        println!(
            "  {:>12.0} {:>7} {:>8.2} {:>8.2} {:>7.0}  {}",
            b.area_fraction * sphere_km2,
            b.cells,
            lat,
            lon,
            b.lowest,
            if b.world_ocean {
                "world ocean"
            } else if b.sea {
                "sea"
            } else {
                "land"
            }
        );
    }
    let (dry, dry_area) = mask
        .bodies
        .iter()
        .filter(|b| !b.sea)
        .fold((0, 0.0), |(n, a), b| (n + 1, a + b.area_fraction));
    println!(
        "  {dry} bodies ({:.0} km2) become land; {} enclosed bodies stay sea; -> {}",
        dry_area * sphere_km2,
        mask.bodies
            .iter()
            .filter(|b| b.sea && !b.world_ocean)
            .count(),
        path.display()
    );

    // The land fraction with land as `elevation > 0`, as before the mask.
    let land = world.fields.read(world.core.terrain.land_fraction).to_vec();
    let above: Vec<f64> = elevation
        .iter()
        .map(|&z| if z > 0.0 { 1.0 } else { 0.0 })
        .collect();
    let mut old = vec![0.0; land.len()];
    world.two_level.restrict_mean(&above, &mut old);
    // As `TerrainStats` does, clamp away round-off.
    old.iter_mut().for_each(|f| *f = f.clamp(0.0, 1.0));
    let cell_areas = &world.climate.geometry.cell_areas;
    let mean = |f: &[f64]| {
        f.iter().zip(cell_areas).map(|(f, a)| f * a).sum::<f64>() / cell_areas.iter().sum::<f64>()
    };
    let (worst, change) =
        (0..land.len())
            .map(|c| (c, land[c] - old[c]))
            .fold((0, 0.0), |best, (c, d)| {
                if d.abs() > f64::abs(best.1) {
                    (c, d)
                } else {
                    best
                }
            });
    let (lat, lon) = lat_lon(world.climate.mesh.cell_centers()[worst]);
    println!(
        "land fraction: global {:.5} with elevation > 0, {:.5} with the sea mask; \
         largest change in a climate cell {change:+.4} (cell {worst} at {lat:.1}, {lon:.1}: \
         {:.4} -> {:.4}); {} cells changed",
        mean(&old),
        mean(&land),
        old[worst],
        land[worst],
        (0..land.len()).filter(|&c| land[c] != old[c]).count()
    );

    // Basins.
    let ops = Operators::new(&world.climate.mesh, &world.climate.geometry, r);
    let t = Instant::now();
    let basins = Basins::new(&ops, &land);
    println!(
        "basins: {} kept, {} dropped as too small ({:.1} ms)",
        basins.basins().len(),
        basins.dropped().len(),
        t.elapsed().as_secs_f64() * 1e3
    );
    let vertex_areas = &world.climate.geometry.vertex_areas;
    let mut csv = String::from(
        "basin,kept,area_km2,interior_vertices,landmasses,touches_coast,\
         latitude_deg,longitude_deg,lowest_latitude_deg,highest_latitude_deg\n",
    );
    println!(
        "  {:>5} {:>12} {:>8} {:>10} {:>5} {:>8} {:>8}  latitudes",
        "basin", "km2", "vertices", "landmasses", "coast", "lat", "lon"
    );
    let listed = basins
        .basins()
        .iter()
        .enumerate()
        .map(|(k, b)| (format!("{}", k + 1), true, b))
        .chain(basins.dropped().iter().map(|b| ("-".to_string(), false, b)));
    for (label, kept, b) in listed {
        let ((lat, lon), (low, high)) = extent(&world.climate.mesh, vertex_areas, b);
        let _ = writeln!(
            csv,
            "{label},{kept},{:.0},{},{},{},{lat:.2},{lon:.2},{low:.2},{high:.2}",
            b.area / 1e6,
            b.vertices.len(),
            b.landmasses,
            b.touches_coast
        );
        println!(
            "  {label:>5} {:>12.0} {:>8} {:>10} {:>5} {lat:>8.2} {lon:>8.2}  {low:.1} to {high:.1}",
            b.area / 1e6,
            b.vertices.len(),
            b.landmasses,
            if b.touches_coast { "yes" } else { "no" },
        );
    }
    let path = dir.join("basins.csv");
    std::fs::write(&path, csv).map_err(|e| format!("cannot write {}: {e}", path.display()))?;
    let mut csv = String::from(
        "cell,latitude_deg,longitude_deg,land_fraction_above_sea_level,land_fraction,wet,basin\n",
    );
    for (c, &p) in world.climate.mesh.cell_centers().iter().enumerate() {
        let (lat, lon) = lat_lon(p);
        let _ = writeln!(
            csv,
            "{c},{lat:.3},{lon:.3},{},{},{},{}",
            old[c],
            land[c],
            basins.is_wet(c as u32),
            basins.cell_basin(c as u32).map_or(0, |b| b + 1)
        );
    }
    let path = dir.join("basin_cells.csv");
    std::fs::write(&path, csv).map_err(|e| format!("cannot write {}: {e}", path.display()))?;
    println!("  -> {}, basins.csv", path.display());

    // Maps: the mask, the land fraction and the basins.
    let mut world = world;
    let basin = world
        .fields
        .register_with(
            FieldSpec {
                description: "Ocean basin of the wind-driven gyres, numbered from 1 by area \
                              (0 where there is none)",
                color_scale: ColorScale::Categorical,
                persistent: false,
                ..FieldSpec::new("ocean.basin", "", Location::CLIMATE_CELLS)
            },
            (0..land.len())
                .map(|c| basins.cell_basin(c as u32).map_or(0.0, |b| b as f64 + 1.0))
                .collect(),
        )
        .map_err(|e| e.to_string())?;
    for key in [
        world.core.terrain.sea.key(),
        world.core.terrain.land_fraction.key(),
        basin.key(),
    ] {
        write_one(&world, key, None, &args.image)?;
    }
    Ok(())
}
