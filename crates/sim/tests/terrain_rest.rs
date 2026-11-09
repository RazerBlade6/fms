//! Terrain in the climate run (M2 step 10.2.1): a stratified atmosphere at
//! rest stays nearly at rest over terrain, measuring the spurious winds of
//! the σ coordinates.

use std::path::Path;

use fms_atmos::column_physics::{CLIMATE_LEVELS, initial_state};
use fms_atmos::{
    Forcing, Gas, HyperdiffusionOptions, LayeredModel, LayeredSemiImplicit,
    LayeredSemiImplicitOptions, LayeredState, Levels, PhysicsParams,
};
use fms_grid::Operators;
use fms_sim::{World, WorldConfig};

const DAY: f64 = 86_400.0;
/// Largest persistent wind over the idealized continent (m/s).
const CONTINENT_LIMIT: f64 = 1.0;

fn world(terrain_frequency: u32, climate_frequency: u32) -> World {
    World::new(WorldConfig {
        terrain_frequency,
        climate_frequency,
        ..WorldConfig::default()
    })
}

fn model(world: &World) -> LayeredModel<'_> {
    let params = world.planet.params();
    let omega = world.planet.rotation_rate();
    let (mesh, geometry) = (&world.climate.mesh, &world.climate.geometry);
    LayeredModel::new(
        Operators::new(mesh, geometry, params.radius),
        params.gravity,
        Gas::EARTH_AIR,
        Levels::new(&CLIMATE_LEVELS),
        mesh.vertices().iter().map(|p| 2.0 * omega * p.z).collect(),
        world.surface_geopotential(),
    )
    .with_humidity_profile(PhysicsParams::default().radiation)
}

/// Rayleigh damping of every layer's wind.
struct Damping(f64);

impl Forcing for Damping {
    fn add_tendency(&mut self, state: &LayeredState, out: &mut LayeredState) {
        for (o, u) in out.u.iter_mut().zip(&state.u) {
            for (o, u) in o.iter_mut().zip(u) {
                *o -= self.0 * u;
            }
        }
    }
}

/// The largest |u| (m/s) after `days` of damping the resting start with a
/// one-day timescale, and the relative change in air mass. The damping
/// removes the start's gravity waves, leaving the wind the terrain keeps
/// forcing: a persistent pressure-gradient error of `a` m/s per day holds a
/// wind of `a` m/s.
fn persistent_wind(world: &World, days: f64) -> (f64, f64) {
    let model = model(world);
    let params = PhysicsParams::default();
    let (state, _) = initial_state(&model, &world.planet, &params);
    let dt = LayeredSemiImplicit::stable_dt(&model, 120.0, 0.6);
    let dt = DAY / (DAY / dt).ceil();
    let diffusion = HyperdiffusionOptions::strongest_for_step(model.operators(), 4, dt);
    let model = model
        .with_hyperdiffusion(diffusion)
        .with_forcing(Box::new(Damping(1.0 / DAY)));
    let mass = model.diagnostics(&state).mass;
    let mut si =
        LayeredSemiImplicit::new(model, state, dt, LayeredSemiImplicitOptions::default(), 0.5);
    while (si.steps() as f64) * dt < days * DAY {
        si.step();
        si.check().expect("the resting atmosphere blew up");
    }
    let wind = si
        .state()
        .u
        .iter()
        .flatten()
        .fold(0.0f64, |a, &u| a.max(u.abs()));
    let end = si.model().diagnostics(si.state()).mass;
    (wind, (end - mass).abs() / mass)
}

#[test]
fn a_flat_world_starts_at_the_planets_surface_pressure() {
    let w = world(16, 8);
    let model = model(&w);
    let (state, _) = initial_state(&model, &w.planet, &PhysicsParams::default());
    let p0 = w.planet.params().surface_pressure;
    assert!(state.ps.iter().all(|&p| p == p0));
}

#[test]
fn surface_pressure_falls_with_height_over_the_continent() {
    let mut w = world(32, 16);
    w.set_idealized_continent();
    let model = model(&w);
    let (state, _) = initial_state(&model, &w.planet, &PhysicsParams::default());
    let p0 = w.planet.params().surface_pressure;
    let (mut highest, mut p_highest) = (0.0, p0);
    for (&phi, &p) in model.surface_geopotential().iter().zip(&state.ps) {
        assert!(phi >= 0.0);
        if phi == 0.0 {
            assert_eq!(p, p0);
        } else {
            assert!(p < p0);
        }
        if phi > highest {
            (highest, p_highest) = (phi, p);
        }
    }
    let z = highest / w.planet.params().gravity;
    // Against a 280 K scale height, within 3%.
    let expected = p0 * fms_math::exp(-highest / (287.0 * 280.0));
    eprintln!(
        "highest cell {z:.0} m: {:.1} hPa (280 K: {:.1})",
        p_highest / 100.0,
        expected / 100.0
    );
    assert!(z > 1000.0, "the ridge should reach the climate grid: {z} m");
    assert!((p_highest / expected - 1.0).abs() < 0.03);
}

#[test]
fn a_stratified_atmosphere_stays_nearly_at_rest_over_the_continent() {
    let flat = world(32, 16);
    let mut continent = world(32, 16);
    continent.set_idealized_continent();
    let (flat_wind, _) = persistent_wind(&flat, 15.0);
    let (wind, mass) = persistent_wind(&continent, 15.0);
    eprintln!("persistent wind: flat {flat_wind:.3} m/s, continent {wind:.3} m/s; mass {mass:.1e}");
    assert!(flat_wind < 0.05, "flat: {flat_wind} m/s");
    assert!(wind < CONTINENT_LIMIT, "continent: {wind} m/s");
    assert!(mass < 1e-12, "mass changed by {mass}");
}

/// The same on Earth's terrain; needs `scripts/fetch-etopo.sh`.
#[test]
#[ignore = "needs data/ETOPO_2022_v1_60s_N90W180_surface.tif; run with --ignored (~1 minute)"]
fn a_stratified_atmosphere_stays_nearly_at_rest_over_earth() {
    let path = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../data/ETOPO_2022_v1_60s_N90W180_surface.tif");
    let mut w = world(64, 32);
    w.load_heightmap(&path, None).expect("read ETOPO");
    let phi = w.surface_geopotential();
    let land = w.fields.read(w.core.terrain.land_fraction);
    let g = w.planet.params().gravity;
    let ocean: Vec<f64> = phi
        .iter()
        .zip(land.iter())
        .filter(|&(_, &f)| f == 0.0)
        .map(|(p, _)| p / g)
        .collect();
    let lifted = |z: f64| ocean.iter().filter(|&&h| h > z).count();
    eprintln!(
        "all-ocean cells: {}, above 100 m {}, above 500 m {}",
        ocean.len(),
        lifted(100.0),
        lifted(500.0)
    );
    let (wind, mass) = persistent_wind(&w, 15.0);
    eprintln!("Earth, 10k cells: persistent wind {wind:.3} m/s; mass {mass:.1e}");
    assert!(wind < 1.5, "{wind} m/s");
    assert!(mass < 1e-12, "mass changed by {mass}");
}
