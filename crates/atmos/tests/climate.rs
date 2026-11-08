//! Radiation and moisture (M2 steps 8–9) in 3-D: the two-layer model over a
//! slab aquaplanet with Earth's orbit. The planet's energy and water budgets
//! close, and (slow, run with `--ignored`) the equilibrium climate is that of
//! an ice-free, Earth-like aquaplanet.

use fms_atmos::column_physics::{LATENT_HEAT, initial_state};
use fms_atmos::williamson::{GRAVITY, RADIUS};
use fms_atmos::{
    ColumnPhysics, Gas, HyperdiffusionOptions, PhysicsParams, TwoLayer, TwoLayerSemiImplicit,
    TwoLayerSemiImplicitOptions,
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
    let (state, ts) = initial_state(&model, &planet, &params);
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
        vec![ts; mesh.n_cells()],
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
        si.check().expect("the moist climate blew up");
        each(&si, &geometry, dt);
    }
}

#[test]
fn the_planets_energy_and_water_budgets_close() {
    // What the atmosphere (enthalpy and latent energy) and the slab gain
    // equals sunlight absorbed minus longwave emitted, minus the kinetic
    // energy friction and diffusion remove; only the dynamics'
    // time-stepping error is left. The atmosphere's water changes by
    // evaporation minus precipitation (and the time filter's trace), to
    // round-off.
    let mut first = None;
    let mut last = None;
    run(16, 90.0, |si, _, dt| {
        let budget = si.model().forcing().unwrap().budget().unwrap();
        let row = (
            si.model().diagnostics(si.state()),
            si.energy_budget(),
            budget,
            si.steps() as f64 * dt,
            si.water(),
            si.water_budget(),
        );
        if first.is_none() {
            first = Some(row);
        }
        last = Some(row);
    });
    let ((d0, e0, b0, t0, w0, wb0), (d1, e1, b1, t1, w1, wb1)) = (first.unwrap(), last.unwrap());
    let area = 4.0 * std::f64::consts::PI * RADIUS * RADIUS;
    let wm2 = |j: f64| j / (t1 - t0) / area;
    let toa = (b1.absorbed - b0.absorbed) - (b1.outgoing - b0.outgoing);
    let uptake = (d1.energy - d0.energy) + (b1.stored - b0.stored) + LATENT_HEAT * (w1 - w0);
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

    let exchanged = (b1.evaporated - b0.evaporated) - (b1.precipitated - b0.precipitated);
    let (filter, fixer) = (wb1.filter - wb0.filter, wb1.fixer - wb0.fixer);
    let water_residual = (w1 - w0) - exchanged - filter - fixer;
    let mm_day = |kg: f64| kg / area / (t1 - t0) * DAY;
    eprintln!(
        "water: {:.1} kg/m², evaporation − precipitation {:+.4} mm/day, filter {:+.1e}, fixer {:+.1e}, residual {:+.1e} of the water",
        w1 / area,
        mm_day(exchanged),
        filter / w1,
        fixer / w1,
        water_residual / w1
    );
    assert!(
        water_residual.abs() < 1e-12 * w1,
        "water residual {water_residual:e} kg"
    );
    assert!(filter.abs() < 1e-3 * w1 && fixer.abs() < 1e-9 * w1);
}

#[test]
#[ignore = "slow: five model years (~6 minutes); run with --ignored"]
fn the_equilibrium_climate_is_an_earthlike_aquaplanet() {
    let planet = Planet::earth();
    let year = planet.year_seconds();
    let years = 5.0;
    let mut start = None;
    let mut sums: Vec<f64> = Vec::new();
    let mut names: Vec<&str> = Vec::new();
    let mut n = 0.0;
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
            si.water(),
        );
        if start.is_none() {
            start = Some(row);
        }
        if si.steps() % 8 == 0 {
            let fields = forcing.cell_fields();
            let total: f64 = geometry.cell_areas.iter().sum();
            if sums.is_empty() {
                sums = vec![0.0; fields.len()];
                names = fields.iter().map(|f| f.name).collect();
            }
            for (s, f) in sums.iter_mut().zip(&fields) {
                *s += f
                    .values
                    .iter()
                    .zip(&geometry.cell_areas)
                    .map(|(v, a)| v * a)
                    .sum::<f64>()
                    / total;
            }
            n += 1.0;
        }
        end = Some(row);
    });
    let ((b0, d0, e0, t0, w0), (b1, d1, e1, t1, w1)) = (start.unwrap(), end.unwrap());
    let area = 4.0 * std::f64::consts::PI * RADIUS * RADIUS;
    let wm2 = |j: f64| j / (t1 - t0) / area;
    let (asr, olr) = (
        wm2(b1.absorbed - b0.absorbed),
        wm2(b1.outgoing - b0.outgoing),
    );
    let albedo = 1.0 - asr / (planet.mean_flux() / 4.0);
    let uptake = wm2((d1.energy - d0.energy) + (b1.stored - b0.stored) + LATENT_HEAT * (w1 - w0));
    let removed =
        -wm2((e1.forcing.kinetic - e0.forcing.kinetic)
            + (e1.diffusion.total() - e0.diffusion.total()));
    let mean = |name: &str| sums[names.iter().position(|n| *n == name).unwrap()] / n;
    let ts = mean("rad.surface_temperature");
    let (precipitation, cloud) = (mean("atm.precipitation"), mean("atm.cloud_cover"));
    let rh = [
        mean("atm.relative_humidity_upper"),
        mean("atm.relative_humidity_lower"),
    ];
    eprintln!(
        "final year: T_s {ts:.2} K, ASR {asr:.1}, OLR {olr:.1}, albedo {albedo:.3}, cloud {cloud:.2}, \
         RH {:.2}/{:.2}, P {precipitation:.2} mm/day, water {:.1} kg/m²; \
         ASR − OLR {:+.2} W/m² = uptake {uptake:+.2} + kinetic energy removed {removed:.2}",
        rh[0],
        rh[1],
        mean("atm.water_vapour"),
        asr - olr
    );
    // An ice-free aquaplanet with Earth's clouds runs warmer than Earth
    // (288 K), whose ice, snow and land brighten the clear sky; measured
    // 297.1 K (O-23). Earth: 240 W/m² each way, albedo 0.29–0.30, cloud
    // 0.6–0.7; the humidity targets are the calibrated radiation's.
    // Measured: albedo 0.286, cloud 0.63, RH 0.31/0.63, P 4.8 mm/day.
    assert!((292.0..300.0).contains(&ts), "T_s {ts}");
    assert!((225.0..250.0).contains(&olr), "OLR {olr}");
    assert!((0.27..0.33).contains(&albedo), "albedo {albedo}");
    assert!((0.55..0.72).contains(&cloud), "cloud {cloud}");
    assert!(
        (rh[0] - 0.31).abs() < 0.04 && (rh[1] - 0.63).abs() < 0.04,
        "RH {rh:?}"
    );
    // Too wet (O-24), but bounded.
    assert!((3.0..6.0).contains(&precipitation), "P {precipitation}");
    // Nearly settled: little heat still going into the atmosphere and slab
    // (measured 0.5 W/m²; the moist aquaplanet settles more slowly than the
    // fixed-humidity one did). The top of the atmosphere doesn't balance
    // exactly: friction and diffusion remove kinetic energy without
    // returning it as heat (O-21). The budget closes to the leapfrog's
    // time-stepping error, ~0.1 W/m² with moist physics (measured 0.12).
    assert!(uptake.abs() < 1.0, "still drifting: {uptake} W/m²");
    assert!((asr - olr - uptake - removed).abs() < 0.2);
}
