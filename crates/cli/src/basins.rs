//! `fms-climate basins`: the world's sea mask and the ocean basins of the
//! wind-driven gyres, as maps and tables, without running anything; with
//! `--gyres`, the gyres themselves for an idealized wind.

use std::fmt::Write as _;
use std::path::Path;
use std::time::Instant;

use clap::Args;
use fms_fields::{ColorScale, FieldSpec, Location};
use fms_grid::CellId;
use fms_grid::coords::to_lat_lon;
use fms_grid::{Operators, SphereMesh};
use fms_ocean::{Basin, Basins, GyreParams, Gyres};
use fms_worldgen::SeaMask;
use glam::DVec3;

use crate::colormap::{Colormap, MISSING};
use crate::render::{CellRaster, Legend, format_value, image_height, render, write_png};
use crate::{ImageArgs, WorldArgs, build_world, ensure_dir, write_one};

#[derive(Args)]
pub struct BasinsArgs {
    #[command(flatten)]
    world: WorldArgs,
    /// Bodies below sea level to print besides the world ocean, largest
    /// first (all of them go to `sea_bodies.csv`).
    #[arg(long, default_value_t = 40)]
    bodies: usize,
    /// Also solve the wind-driven gyres for an idealized zonal wind stress,
    /// the same in both hemispheres: `τ_x = cos φ · Σ a (g(φ) + g(−φ))`,
    /// `g(φ) = exp(−((φ − φ₀)/σ)²)`, with trades (a = −0.062 N/m², φ₀ = 15°,
    /// σ = 13°), westerlies (0.149, 48°, 12°) and polar easterlies (−0.097,
    /// 72°, 10°): −0.03 N/m² at the equator, −0.06 at 15°, +0.10 at 47° and
    /// −0.03 at 71°, taken along each mesh edge's normal at its midpoint.
    /// Prints each component's solve, the largest coast groups' streamfunction
    /// and the flow through Drake Passage, and writes
    /// `ocean.gyre_streamfunction.png` (Sv), `gyre_vertices.csv` and
    /// `gyre_groups.csv`.
    #[arg(long)]
    gyres: bool,
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
    if args.gyres {
        gyres(&world, &ops, &basins, &args.image)?;
    }

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

/// The trades, westerlies and polar easterlies of `--gyres`: amplitude
/// (N m⁻²), centre and width (degrees) of each Gaussian.
const WINDS: [(f64, f64, f64); 3] = [
    (-0.062, 15.0, 13.0),
    (0.149, 48.0, 12.0),
    (-0.097, 72.0, 10.0),
];

/// `G(φ)` and `∂G/∂φ` (N m⁻² rad⁻¹) of the idealized zonal stress
/// `τ_x = cos φ · G(φ)` of `--gyres`, at latitude `lat` (radians).
fn wind_profile(lat: f64) -> (f64, f64) {
    let (mut g, mut slope) = (0.0, 0.0);
    for (a, centre, width) in WINDS {
        let width = width.to_radians();
        for centre in [centre.to_radians(), -centre.to_radians()] {
            let x = (lat - centre) / width;
            let e = a * fms_math::exp(-x * x);
            g += e;
            slope -= 2.0 * x / width * e;
        }
    }
    (g, slope)
}

/// The curl (N m⁻³) at latitude `lat` (radians) of the idealized zonal stress
/// `τ_x = cos φ · G(φ)` of `--gyres` on a planet of radius `r` (m):
/// `−(1/(r cos φ)) ∂(τ_x cos φ)/∂φ = (2 sin φ · G − cos φ · ∂G/∂φ)/r`, finite
/// at the poles. Only written out, for comparison: the solve takes the
/// stress.
fn wind_curl(lat: f64, r: f64) -> f64 {
    let (g, slope) = wind_profile(lat);
    (2.0 * fms_math::sin(lat) * g - fms_math::cos(lat) * slope) / r
}

/// The idealized stress of `--gyres` along each mesh edge's normal (N m⁻²):
/// `τ_x` at the edge point times its eastward unit vector's component along
/// the normal.
fn wind_stress(ops: &Operators) -> Vec<f64> {
    let geometry = ops.geometry();
    (0..ops.mesh().n_edges())
        .map(|e| {
            let p = geometry.edge_points[e];
            let lat = to_lat_lon(p).0;
            let east = DVec3::Z.cross(p).normalize_or_zero();
            fms_math::cos(lat) * wind_profile(lat).0 * east.dot(geometry.edge_normals[e])
        })
        .collect()
}

/// The cell of `mesh` nearest `(lat, lon)` (degrees).
fn nearest_cell(mesh: &SphereMesh, lat: f64, lon: f64) -> CellId {
    let p = fms_grid::coords::from_lat_lon(lat.to_radians(), lon.to_radians());
    let centres = mesh.cell_centers();
    (0..centres.len())
        .max_by(|&a, &b| centres[a].dot(p).total_cmp(&centres[b].dot(p)))
        .expect("a mesh has cells") as CellId
}

/// `--gyres`: solve, print each component and the largest groups, and draw
/// the streamfunction.
fn gyres(
    world: &fms_sim::World,
    ops: &Operators,
    basins: &Basins,
    image: &ImageArgs,
) -> Result<(), String> {
    let (mesh, r) = (ops.mesh(), ops.radius());
    let rotation = world.planet.rotation_rate();
    let params = GyreParams::default();
    let t = Instant::now();
    let mut gyres = Gyres::new(ops, basins, rotation, params);
    let built = t.elapsed().as_secs_f64();
    let stress = wind_stress(ops);
    let t = Instant::now();
    gyres.solve(&stress);
    let solved = t.elapsed().as_secs_f64();
    println!(
        "gyres: friction {:.1e} s-1, rotation {rotation:.4e} rad/s; {} coast groups; \
         matrices built in {:.1} ms, solved in {:.1} ms",
        params.friction,
        basins.groups().len(),
        built * 1e3,
        solved * 1e3
    );
    const SV: f64 = 1e9;
    let psi = gyres.streamfunction();
    let group_psi = gyres.group_streamfunction();
    println!(
        "  {:>4} {:>6} {:>6} {:>8} {:>7} {:>10} {:>9} {:>8} {:>16} {:>8} {:>16}",
        "comp",
        "basins",
        "groups",
        "vertices",
        "unknown",
        "iterations",
        "residual",
        "min Sv",
        "at",
        "max Sv",
        "at"
    );
    for (k, (component, stats)) in gyres.components().iter().zip(gyres.stats()).enumerate() {
        let vertices: usize = component
            .basins
            .iter()
            .map(|&b| basins.basins()[b].vertices.len())
            .sum();
        let head = format!(
            "  {:>4} {:>6} {:>6} {:>8}",
            k + 1,
            component.basins.len(),
            component.groups.len(),
            vertices
        );
        let Some(stats) = stats else {
            println!("{head}  not solved: no coast");
            continue;
        };
        let at = |&v: &u32| lat_lon(mesh.vertices()[v as usize]);
        let by_psi = |a: &&u32, b: &&u32| psi[**a as usize].total_cmp(&psi[**b as usize]);
        let all = || {
            component
                .basins
                .iter()
                .flat_map(|&b| basins.basins()[b].vertices.iter())
        };
        let extremes = all().min_by(by_psi).zip(all().max_by(by_psi));
        let extremes = extremes.map_or(String::new(), |(low, high)| {
            let ((low_lat, low_lon), (high_lat, high_lon)) = (at(low), at(high));
            format!(
                " {:>8.1} {:>16} {:>8.1} {:>16}",
                psi[*low as usize] / SV,
                format!("{low_lat:.1}, {low_lon:.1}"),
                psi[*high as usize] / SV,
                format!("{high_lat:.1}, {high_lon:.1}"),
            )
        });
        println!(
            "{head} {:>7} {:>10} {:>9.2e}{extremes}{}",
            format!("{}+{}", component.vertex_unknowns, component.group_unknowns),
            stats.iterations,
            stats.residual,
            if stats.converged {
                ""
            } else {
                "  NOT CONVERGED"
            }
        );
    }

    // The coast groups: Ψ, land area and the land's area-weighted centre.
    let cell_areas = &world.climate.geometry.cell_areas;
    let mut centres = vec![DVec3::ZERO; basins.groups().len()];
    for (c, (&area, &p)) in cell_areas.iter().zip(mesh.cell_centers()).enumerate() {
        if let Some(g) = basins.cell_group(c as CellId) {
            centres[g] += area * p;
        }
    }
    let mut csv = String::from(
        "group,component,reference,psi_sv,land_km2,landmasses,vertices,latitude_deg,longitude_deg\n",
    );
    let shown = 10;
    println!(
        "  {:>5} {:>4} {:>9} {:>12} {:>10} {:>8} {:>8} {:>8}  (largest {shown} by land area)",
        "group", "comp", "psi Sv", "land km2", "landmasses", "vertices", "lat", "lon"
    );
    for (k, group) in basins.groups().iter().enumerate() {
        let (lat, lon) = lat_lon(centres[k].normalize_or_zero());
        let component = gyres.group_component(k);
        let reference = gyres.components()[component].reference == Some(k);
        let _ = writeln!(
            csv,
            "{},{},{reference},{},{:.0},{},{},{lat:.2},{lon:.2}",
            k + 1,
            component + 1,
            group_psi[k] / SV,
            group.area / 1e6,
            group.landmasses.len(),
            group.vertices.len()
        );
        if k < shown {
            println!(
                "  {:>5} {:>4} {:>9.2} {:>12.0} {:>10} {:>8} {lat:>8.2} {lon:>8.2}{}",
                k + 1,
                component + 1,
                group_psi[k] / SV,
                group.area / 1e6,
                group.landmasses.len(),
                group.vertices.len(),
                if reference { "  reference" } else { "" }
            );
        }
    }
    let path = image.out_dir.join("gyre_groups.csv");
    std::fs::write(&path, csv).map_err(|e| format!("cannot write {}: {e}", path.display()))?;

    // Drake Passage, for Earth-like terrain: Ψ of Antarctica's group less
    // South America's is the eastward flow between them.
    let group_at = |lat: f64, lon: f64| basins.cell_group(nearest_cell(mesh, lat, lon));
    match (group_at(-80.0, 0.0), group_at(-15.0, -60.0)) {
        (Some(south), Some(north)) if south != north => println!(
            "  Drake Passage (diagnostic, Earth-like terrain only): group {} (land at 80S 0E) \
             minus group {} (land at 15S 60W) = {:+.1} Sv, positive eastward",
            south + 1,
            north + 1,
            (group_psi[south] - group_psi[north]) / SV
        ),
        (Some(_), Some(_)) => println!(
            "  Drake Passage (diagnostic): the land at 80S 0E and at 15S 60W is one group: closed"
        ),
        _ => println!("  Drake Passage (diagnostic): no land at 80S 0E or at 15S 60W"),
    }

    let mut csv = String::from("vertex,latitude_deg,longitude_deg,basin,group,curl_n_m3,psi_sv\n");
    for (v, &p) in mesh.vertices().iter().enumerate() {
        let (lat, lon) = lat_lon(p);
        let _ = writeln!(
            csv,
            "{v},{lat:.3},{lon:.3},{},{},{:e},{}",
            basins.vertex_basin(v as u32).map_or(0, |b| b + 1),
            basins.vertex_group(v as u32).map_or(0, |g| g + 1),
            wind_curl(lat.to_radians(), r),
            psi[v] / SV
        );
    }
    let path = image.out_dir.join("gyre_vertices.csv");
    std::fs::write(&path, csv).map_err(|e| format!("cannot write {}: {e}", path.display()))?;
    println!("  -> {}, gyre_groups.csv", path.display());

    // The map: each wet cell shows the mean of its vertices, land is grey.
    let values: Vec<f64> = (0..mesh.n_cells())
        .map(|c| {
            let c = c as u32;
            if !basins.is_wet(c) {
                return f64::NAN;
            }
            let ring = mesh.cell_vertices(c);
            ring.iter().map(|&v| psi[v as usize]).sum::<f64>() / ring.len() as f64 / SV
        })
        .collect();
    write_gyre_map(
        &world.climate,
        &values,
        "ocean.gyre_streamfunction",
        "Sv",
        "",
        image,
        &image.out_dir,
    )
}

/// The share of wet cells whose size a gyre map's colours span; the rest
/// saturate. Not higher: on Earth the circumpolar current is a tenth of the
/// wet cells and four times the subtropical gyres.
const GYRE_MAP_QUANTILE: f64 = 0.9;

/// A diverging map of a gyre field, `<dir>/<name>.png`. `values` is NaN
/// where a cell isn't wet ([`Basins::is_wet`]), drawn grey. Zero is white
/// and the range symmetric: `image.range` if given, else ± the
/// [`GYRE_MAP_QUANTILE`] quantile of the wet cells' |value|, beyond which the
/// colours saturate, so one strong current doesn't wash out the rest.
pub fn write_gyre_map(
    grid: &fms_sim::Grid,
    values: &[f64],
    name: &str,
    units: &str,
    stamp: &str,
    image: &ImageArgs,
    dir: &Path,
) -> Result<(), String> {
    let mut sizes: Vec<f64> = values
        .iter()
        .filter(|v| v.is_finite())
        .map(|v| v.abs())
        .collect();
    sizes.sort_by(f64::total_cmp);
    let largest = sizes.last().copied().unwrap_or(0.0);
    let (min, max, title) = match image.range {
        Some((min, max)) => (min, max, format!("{name} ({units}), land grey")),
        None => {
            let k = ((GYRE_MAP_QUANTILE * sizes.len() as f64).ceil() as usize).max(1) - 1;
            let cap = [sizes.get(k).copied().unwrap_or(0.0), largest, 1.0]
                .into_iter()
                .find(|&c| c > 0.0)
                .unwrap_or(1.0);
            let title = format!(
                "{name} ({units}), land grey, range: {:.0}th percentile of magnitude (largest {})",
                GYRE_MAP_QUANTILE * 100.0,
                format_value(largest)
            );
            (-cap, cap, title)
        }
    };
    let colormap = Colormap::new(fms_fields::ColorScale::Diverging, min, max);
    let raster = CellRaster::new(grid, image.width, false);
    let legend = Legend {
        title: &title,
        subtitle: stamp,
    };
    let pixels = render(&raster, values, &colormap, &legend, None);
    let mut palette = *colormap.palette();
    palette[MISSING as usize] = [150, 150, 150];
    let path = dir.join(format!("{name}.png"));
    write_png(
        &path,
        raster.width,
        image_height(&raster),
        &pixels,
        &palette,
    )
    .map_err(|e| format!("cannot write {}: {e}", path.display()))?;
    println!(
        "{name} [{} .. {}], largest {}, land grey -> {}",
        format_value(min),
        format_value(max),
        format_value(largest),
        path.display()
    );
    Ok(())
}
