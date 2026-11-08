//! Radiation (M2 step 8) in 3-D: the two-layer model over a slab aquaplanet
//! with Earth's orbit. The planet's energy budget closes, and (slow, run with
//! `--ignored`) the equilibrium climate is Earth-like.

use fms_atmos::williamson::{GRAVITY, RADIUS};
use fms_atmos::{
    ColumnPhysics, Gas, HyperdiffusionOptions, PhysicsParams, SingleColumn, TwoLayer,
    TwoLayerSemiImplicit, TwoLayerSemiImplicitOptions, TwoLayerState, column_physics::layer_exner,
};
use fms_geo::{LloydOptions, MeshGeometry, Operators, SphereMesh};
use fms_planet::Planet;

const DAY: f64 = 86_400.0;

/// Run `days` on frequency `k` from the single-column equilibrium, calling
/// `each` after every step.
fn run(k: u32, days: f64, mut each: impl FnMut(&TwoLayerSemiImplicit, &MeshGeometry, f64)) {
    let planet = Planet::earth();
    let (mesh, _) = SphereMesh::centroidal(k, LloydOptions::default());
    let geometry = MeshGeometry::new(&mesh);
    let params = PhysicsParams::default();
    let mut column = SingleColumn::earth(planet.mean_flux());
    column.params.slab_depth = 1.0;
    let eq = column.equilibrium();
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
    let exner = layer_exner(Gas::EARTH_AIR, 0.5, ps);
    let mut x = 0x2545_F491_4F6C_DD1Du64;
    let theta: [Vec<f64>; 2] = std::array::from_fn(|k| {
        (0..nc)
            .map(|_| {
                x ^= x << 13;
                x ^= x >> 7;
                x ^= x << 17;
                eq.temperature[k] / exner[k] + 0.1 * ((x >> 11) as f64 / (1u64 << 52) as f64 - 1.0)
            })
            .collect()
    });
    let state = TwoLayerState::new(vec![ps; nc], [vec![0.0; ne], vec![0.0; ne]], theta);
    // Radiatively driven jets and storms reach ~100 m/s.
    let dt = TwoLayerSemiImplicit::stable_dt(&model, 120.0, 0.6);
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
    let mut si = TwoLayerSemiImplicit::new(
        model,
        state,
        dt,
        TwoLayerSemiImplicitOptions::default(),
        0.5,
    );
    while (si.steps() as f64) * dt < days * DAY {
        si.step();
        si.check().expect("the radiative climate blew up");
        each(&si, &geometry, dt);
    }
}

#[test]
fn the_planets_energy_budget_closes() {
    // What the atmosphere and the slab gain equals sunlight absorbed minus
    // longwave emitted, minus the kinetic energy friction and diffusion
    // remove. Only the dynamics' time-stepping error is left.
    let mut first = None;
    let mut last = None;
    run(16, 90.0, |si, _, dt| {
        let budget = si.model().forcing().unwrap().budget().unwrap();
        let row = (
            si.model().diagnostics(si.state()),
            si.energy_budget(),
            budget,
            si.steps() as f64 * dt,
        );
        if first.is_none() {
            first = Some(row);
        }
        last = Some(row);
    });
    let ((d0, e0, b0, t0), (d1, e1, b1, t1)) = (first.unwrap(), last.unwrap());
    let area = 4.0 * std::f64::consts::PI * RADIUS * RADIUS;
    let wm2 = |j: f64| j / (t1 - t0) / area;
    let toa = (b1.absorbed - b0.absorbed) - (b1.outgoing - b0.outgoing);
    let uptake = (d1.energy - d0.energy) + (b1.stored - b0.stored);
    let removed =
        (e1.forcing.kinetic - e0.forcing.kinetic) + (e1.diffusion.total() - e0.diffusion.total());
    let residual = toa - uptake + removed;
    eprintln!(
        "90 days: net in at the top {:+.3} W/m², uptake {:+.3}, kinetic energy removed {:.3}, residual {:+.5}",
        wm2(toa),
        wm2(uptake),
        -wm2(removed),
        wm2(residual)
    );
    assert!(
        wm2(residual).abs() < 0.01,
        "residual {} W/m²",
        wm2(residual)
    );
    assert!((d1.mass / d0.mass - 1.0).abs() < 1e-13);
}

#[test]
#[ignore = "slow: five model years (~4 minutes); run with --ignored"]
fn the_equilibrium_climate_is_earthlike() {
    let planet = Planet::earth();
    let year = planet.year_seconds();
    let years = 5.0;
    let mut start = None;
    let (mut ts, mut n) = (0.0, 0.0);
    let mut end = None;
    run(16, years * year / DAY, |si, geometry, dt| {
        let t = si.steps() as f64 * dt;
        if t < (years - 1.0) * year {
            return;
        }
        let forcing = si.model().forcing().unwrap();
        let row = (
            forcing.budget().unwrap(),
            si.model().diagnostics(si.state()),
            si.energy_budget(),
            t,
        );
        if start.is_none() {
            start = Some(row);
        }
        if si.steps() % 8 == 0 {
            let fields = forcing.cell_fields();
            let total: f64 = geometry.cell_areas.iter().sum();
            ts += fields[0]
                .values
                .iter()
                .zip(&geometry.cell_areas)
                .map(|(v, a)| v * a)
                .sum::<f64>()
                / total;
            n += 1.0;
        }
        end = Some(row);
    });
    let ((b0, d0, e0, t0), (b1, d1, e1, t1)) = (start.unwrap(), end.unwrap());
    let area = 4.0 * std::f64::consts::PI * RADIUS * RADIUS;
    let wm2 = |j: f64| j / (t1 - t0) / area;
    let (asr, olr) = (
        wm2(b1.absorbed - b0.absorbed),
        wm2(b1.outgoing - b0.outgoing),
    );
    let albedo = 1.0 - asr / (planet.mean_flux() / 4.0);
    let uptake = wm2((d1.energy - d0.energy) + (b1.stored - b0.stored));
    let removed =
        -wm2((e1.forcing.kinetic - e0.forcing.kinetic)
            + (e1.diffusion.total() - e0.diffusion.total()));
    let ts = ts / n;
    eprintln!(
        "final year: T_s {ts:.2} K, ASR {asr:.1}, OLR {olr:.1}, albedo {albedo:.3}; \
         ASR − OLR {:+.2} W/m² = uptake {uptake:+.2} + kinetic energy removed {removed:.2}",
        asr - olr
    );
    // Earth: 288 K, 240 W/m² each way, albedo 0.30.
    assert!((286.0..290.0).contains(&ts), "T_s {ts}");
    assert!((230.0..250.0).contains(&olr), "OLR {olr}");
    assert!((0.28..0.32).contains(&albedo), "albedo {albedo}");
    // Settled: little heat still going into the atmosphere and slab. The
    // top of the atmosphere doesn't balance exactly: friction and diffusion
    // remove kinetic energy without returning it as heat (O-21).
    assert!(uptake.abs() < 0.5, "still drifting: {uptake} W/m²");
    assert!((asr - olr - uptake - removed).abs() < 0.05);
}
