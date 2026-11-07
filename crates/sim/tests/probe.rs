//! Point probes read the nearest cell, resample non-cell fields like the map
//! renderer, and report vector fields as east/north components.

use fms_fields::{FieldSpec, GridLevel, Location, Quantity};
use fms_geo::coords::from_lat_lon;
use fms_sim::probe::{Component, local_frame};
use fms_sim::{Probe, ProbeError, World, WorldConfig};
use glam::DVec3;

fn small_world() -> World {
    World::new(WorldConfig {
        terrain_frequency: 20,
        climate_frequency: 12,
        ..WorldConfig::default()
    })
}

#[test]
fn picks_the_nearest_cell_on_each_grid() {
    let world = small_world();
    let keys = [
        world.core.climate_latitude.key(),
        world.core.terrain_latitude.key(),
    ];
    let probe = Probe::new(&world, 33.9, 18.4, &keys).unwrap();
    let p = from_lat_lon(33.9f64.to_radians(), 18.4f64.to_radians());
    for (column, grid) in probe
        .columns()
        .iter()
        .zip([GridLevel::Climate, GridLevel::Terrain])
    {
        assert_eq!(column.site.grid, grid);
        assert_eq!(column.component, Component::Value);
        // No other cell centre is closer than the chosen one.
        let centers = world.grid(grid).mesh.cell_centers();
        let chosen = centers[column.site.cell as usize].dot(p);
        assert!(centers.iter().all(|c| c.dot(p) <= chosen + 1e-15));
    }
    // A latitude field read at a cell is that cell's latitude.
    let values = probe.sample(&world);
    for (v, column) in values.iter().zip(probe.columns()) {
        assert!(
            (v - column.site.lat_deg).abs() < 1e-9,
            "{v} vs {:?}",
            column.site
        );
    }
}

#[test]
fn insolation_matches_the_planet_at_the_cell_centre() {
    let mut world = small_world();
    let probe = Probe::new(&world, -10.0, 75.0, &[world.core.toa_insolation.key()]).unwrap();
    let cell = probe.columns()[0].site.cell as usize;
    let center = world.climate.mesh.cell_centers()[cell];
    let dt = world.scheduler.calendar().dt();
    for step in 0..48 {
        world.step();
        let expected = world.planet.insolation(center, step as f64 * dt);
        let got = probe.sample(&world)[0];
        assert_eq!(got.to_bits(), expected.to_bits(), "step {step}");
    }
}

#[test]
fn vector_fields_give_east_and_north_components() {
    let mut world = small_world();
    // Solid-body rotation about the polar axis: u = U cos(lat) eastward, no
    // northward component.
    let speed = 20.0;
    let geometry = &world.climate.geometry;
    let normal: Vec<f64> = geometry
        .edge_points
        .iter()
        .zip(&geometry.edge_normals)
        .map(|(&x, &n)| speed * DVec3::Z.cross(x).dot(n))
        .collect();
    let wind = world
        .fields
        .register_with(
            FieldSpec {
                quantity: Quantity::NormalComponent,
                ..FieldSpec::new("test.wind", "m s-1", Location::CLIMATE_EDGES)
            },
            normal,
        )
        .unwrap();
    for (lat, lon) in [(0.0, 0.0), (45.0, -120.0), (-70.0, 30.0)] {
        let probe = Probe::new(&world, lat, lon, &[wind.key()]).unwrap();
        let labels: Vec<_> = probe.columns().iter().map(|c| c.label.as_str()).collect();
        assert_eq!(labels, ["test.wind.east", "test.wind.north"]);
        let [east, north] = probe.sample(&world)[..] else {
            unreachable!()
        };
        // cos(lat) of the cell centre, without platform trigonometry.
        let c = world.climate.mesh.cell_centers()[probe.columns()[0].site.cell as usize];
        let expected = speed * c.truncate().length();
        // Perot reconstruction is accurate to ~2% at worst (OBSERVATIONS.md O-1).
        assert!(
            (east - expected).abs() < 0.03 * speed,
            "{east} vs {expected}"
        );
        assert!(north.abs() < 0.03 * speed, "north {north}");
    }
}

#[test]
fn vertex_and_edge_scalars_are_averaged_onto_the_cell() {
    let mut world = small_world();
    let vertex = world
        .fields
        .register(
            FieldSpec::new("test.vertex", "1", Location::CLIMATE_VERTICES),
            3.5,
        )
        .unwrap();
    let edge = world
        .fields
        .register(
            FieldSpec::new("test.edge", "1", Location::CLIMATE_EDGES),
            -2.0,
        )
        .unwrap();
    let probe = Probe::new(&world, 60.0, 10.0, &[vertex.key(), edge.key()]).unwrap();
    let values = probe.sample(&world);
    assert!((values[0] - 3.5).abs() < 1e-12);
    assert!((values[1] + 2.0).abs() < 1e-12);

    // Same numbers as the whole-field resampling the renderer uses.
    world.fields.write(vertex)[7] = 100.0;
    let grid = &world.climate;
    let at_cells = grid
        .geometry
        .vertices_to_cells(&grid.mesh, &world.fields.read(vertex));
    let cell = probe.columns()[0].site.cell as usize;
    assert_eq!(probe.sample(&world)[0].to_bits(), at_cells[cell].to_bits());
}

#[test]
fn rejects_invalid_locations() {
    let world = small_world();
    let keys = [world.core.climate_latitude.key()];
    for (lat, lon) in [
        (91.0, 0.0),
        (-90.5, 0.0),
        (f64::NAN, 0.0),
        (0.0, f64::INFINITY),
    ] {
        assert!(matches!(
            Probe::new(&world, lat, lon, &keys),
            Err(ProbeError::InvalidLocation { .. })
        ));
    }
    // Poles and any longitude are fine.
    assert!(Probe::new(&world, 90.0, 720.0, &keys).is_ok());
}

#[test]
fn local_frame_is_orthonormal_including_the_poles() {
    for p in [
        from_lat_lon(0.3, 1.0),
        DVec3::Z,
        -DVec3::Z,
        from_lat_lon(-1.2, -2.5),
    ] {
        let (e, n) = local_frame(p);
        assert!((e.length() - 1.0).abs() < 1e-12 && (n.length() - 1.0).abs() < 1e-12);
        assert!(e.dot(n).abs() < 1e-12 && e.dot(p).abs() < 1e-12 && n.dot(p).abs() < 1e-12);
    }
    // Away from the poles, north points towards +z and east is right-handed.
    let p = from_lat_lon(0.5, 2.0);
    let (e, n) = local_frame(p);
    assert!(n.z > 0.0 && e.cross(n).dot(p) > 0.0);
}
