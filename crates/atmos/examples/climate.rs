//! The two-layer model driven by radiation over a slab aquaplanet (Earth's
//! orbit and CO₂): global means, energy budget and zonal surface
//! temperature.
//!
//! ```sh
//! cargo run --release -p fms-atmos --example climate [frequency] [years]
//! ```
//!
//! Environment: `FMS_CO2` (ppm), `FMS_CLOUD_ALBEDO`.

use std::time::Instant;

use fms_atmos::williamson::{GRAVITY, RADIUS};
use fms_atmos::{
    ColumnPhysics, Gas, HyperdiffusionOptions, PhysicsParams, SingleColumn, TwoLayer,
    TwoLayerSemiImplicit, TwoLayerSemiImplicitOptions, TwoLayerState,
};
use fms_geo::{LatitudeBands, LloydOptions, MeshGeometry, Operators, SphereMesh};
use fms_planet::Planet;

const DAY: f64 = 86_400.0;

fn env(name: &str, default: f64) -> f64 {
    std::env::var(name)
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(default)
}

fn main() {
    let mut args = std::env::args().skip(1);
    let k: u32 = args.next().map_or(16, |s| s.parse().unwrap());
    let years: f64 = args.next().map_or(2.0, |s| s.parse().unwrap());
    let planet = Planet::earth();
    let (mesh, _) = SphereMesh::centroidal(k, LloydOptions::default());
    let geometry = MeshGeometry::new(&mesh);
    let mut params = PhysicsParams::default();
    params.co2 = env("FMS_CO2", params.co2);
    params.radiation.cloud_albedo = env("FMS_CLOUD_ALBEDO", params.radiation.cloud_albedo);

    // Start from the single-column equilibrium, at rest.
    let mut column = SingleColumn::earth(planet.mean_flux());
    column.params = PhysicsParams {
        slab_depth: 1.0,
        ..params
    };
    let eq = column.equilibrium();
    println!(
        "single column: T_s {:.1} K, layers {:.1}/{:.1} K",
        eq.surface_temperature, eq.temperature[0], eq.temperature[1]
    );
    let model = TwoLayer::new(
        Operators::new(&mesh, &geometry, RADIUS),
        GRAVITY,
        Gas::EARTH_AIR,
        0.5,
        mesh.vertices()
            .iter()
            .map(|p| 2.0 * planet.rotation_rate() * p.z)
            .collect(),
        vec![0.0; mesh.n_cells()],
    );
    let (nc, ne) = (mesh.n_cells(), mesh.n_edges());
    let ps = planet.params().surface_pressure;
    let exner = fms_atmos::column_physics::layer_exner(Gas::EARTH_AIR, 0.5, ps);
    let mut x = 0x2545_F491_4F6C_DD1Du64;
    let mut noise = || {
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        0.1 * ((x >> 11) as f64 / (1u64 << 52) as f64 - 1.0)
    };
    let theta: [Vec<f64>; 2] = std::array::from_fn(|k| {
        (0..nc)
            .map(|_| eq.temperature[k] / exner[k] + noise())
            .collect()
    });
    let state = TwoLayerState::new(vec![ps; nc], [vec![0.0; ne], vec![0.0; ne]], theta);
    let dt = TwoLayerSemiImplicit::stable_dt(&model, env("FMS_MAX_WIND", 120.0), 0.6);
    let dt = DAY / (DAY / dt).ceil();
    let diffusion = HyperdiffusionOptions::strongest_for_step(model.operators(), 4, dt);
    let physics = ColumnPhysics::new(
        &model,
        &planet,
        params,
        0.0,
        &state,
        vec![eq.surface_temperature; nc],
    );
    let model = model
        .with_hyperdiffusion(diffusion)
        .with_forcing(Box::new(physics));
    let start = model.diagnostics(&state);
    let budget0 = model.forcing().unwrap().budget().unwrap();
    let mut si = TwoLayerSemiImplicit::new(
        model,
        state,
        dt,
        TwoLayerSemiImplicitOptions::default(),
        0.5,
    );
    let year = planet.year_seconds();
    let per_year = (year / dt).round() as u64;
    let total = (years * per_year as f64) as u64;
    let bands = LatitudeBands::new(&mesh, &geometry, 18);
    let area = 4.0 * std::f64::consts::PI * RADIUS * RADIUS;
    let mut last_year: Option<(f64, f64, f64)> = None;
    let (mut ts_sum, mut ts_zonal, mut samples) = (0.0, vec![0.0; 18], 0.0);
    let t = Instant::now();
    while si.steps() < total {
        si.step();
        if let Err(b) = si.check() {
            println!("BLOW-UP day {:.1}: {b}", si.steps() as f64 * dt / DAY);
            return;
        }
        let forcing = si.model().forcing().unwrap();
        let sampling = si.steps() + per_year >= total && si.steps() % 8 == 0;
        let reporting = si.steps() % (per_year / 8) == 0;
        if !(sampling || reporting || si.steps() + per_year == total) {
            continue;
        }
        let fields = forcing.cell_fields();
        let ts = &fields[0].values;
        if sampling {
            let mean: f64 = ts
                .iter()
                .zip(&geometry.cell_areas)
                .map(|(t, a)| t * a)
                .sum::<f64>()
                / (4.0 * std::f64::consts::PI);
            ts_sum += mean;
            for (z, m) in ts_zonal.iter_mut().zip(bands.mean(&geometry, ts)) {
                *z += m;
            }
            samples += 1.0;
        }
        if si.steps() + per_year == total {
            let b = forcing.budget().unwrap();
            last_year = Some((b.absorbed, b.outgoing, b.stored));
        }
        if reporting {
            let b = forcing.budget().unwrap();
            let d = si.model().diagnostics(si.state());
            let e = si.energy_budget();
            let closure = (d.energy - start.energy) + (b.stored - budget0.stored)
                - (b.absorbed - b.outgoing)
                - e.forcing.kinetic
                - e.diffusion.total();
            let mean_ts: f64 = ts
                .iter()
                .zip(&geometry.cell_areas)
                .map(|(t, a)| t * a)
                .sum::<f64>()
                / (4.0 * std::f64::consts::PI);
            let h = si.check().unwrap();
            println!(
                "day {:>6.1}: T_s {:.2} K  max |u| {:.0}  closure residual {:+.2e} J ({:+.4} W/m² over the run)",
                si.steps() as f64 * dt / DAY,
                mean_ts,
                h.max_wind,
                closure,
                closure / (si.steps() as f64 * dt) / area
            );
        }
    }
    println!("{years} years in {:.0} s", t.elapsed().as_secs_f64());
    if let Some((a0, o0, s0)) = last_year {
        let b = si.model().forcing().unwrap().budget().unwrap();
        let wm2 = |j: f64| j / year / area;
        let asr = wm2(b.absorbed - a0);
        let olr = wm2(b.outgoing - o0);
        println!(
            "last year: T_s {:.2} K, ASR {asr:.1}, OLR {olr:.1}, imbalance {:+.2} W/m², slab uptake {:+.2} W/m², planetary albedo {:.3}",
            ts_sum / samples,
            asr - olr,
            wm2(b.stored - s0),
            1.0 - asr / (planet.mean_flux() / 4.0)
        );
        let lat = bands.latitudes();
        for (l, z) in lat.iter().zip(&ts_zonal) {
            println!("  {:>6.1}°  {:.1} K", l.to_degrees(), z / samples);
        }
    }
}
