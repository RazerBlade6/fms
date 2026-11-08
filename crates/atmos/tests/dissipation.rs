//! Dissipation and safety (M2 step 6): hyperdiffusion keeps budgets closed
//! and leaves balanced rest alone, carries a storm life cycle through a month
//! that the inviscid core can't survive cleanly, and blow-ups are caught
//! within a step.

use std::sync::OnceLock;

use fms_atmos::williamson::{Case, GRAVITY, OMEGA, RADIUS};
use fms_atmos::{
    Fault, Gas, HyperdiffusionOptions, Limits, ShearedFlow, TwoLayer, TwoLayerSemiImplicit,
    TwoLayerSemiImplicitOptions, TwoLayerState, Variable, safety,
};
use fms_geo::coords::to_lat_lon;
use fms_geo::{LloydOptions, MeshGeometry, Operators, SphereMesh};

const GAS: Gas = Gas::EARTH_AIR;
const DAY: f64 = 86_400.0;
/// Design wind for the step: storms on the 30 m/s-shear jet peak near 100 m/s.
const MAX_WIND: f64 = 120.0;

fn mesh() -> &'static (SphereMesh, MeshGeometry) {
    static M: OnceLock<(SphereMesh, MeshGeometry)> = OnceLock::new();
    M.get_or_init(|| {
        let (m, _) = SphereMesh::centroidal(16, LloydOptions::default());
        let g = MeshGeometry::new(&m);
        (m, g)
    })
}

fn bare_model(mountain: bool) -> TwoLayer<'static> {
    let (m, g) = mesh();
    let phi_s = if mountain {
        m.cell_centers()
            .iter()
            .map(|&p| GRAVITY * Case::Mountain.topography(p))
            .collect()
    } else {
        vec![0.0; m.n_cells()]
    };
    TwoLayer::new(
        Operators::new(m, g, RADIUS),
        GRAVITY,
        GAS,
        0.5,
        m.vertices().iter().map(|p| 2.0 * OMEGA * p.z).collect(),
        phi_s,
    )
}

/// The step for [`MAX_WIND`], a whole number per day.
fn step_for(model: &TwoLayer) -> f64 {
    let dt = TwoLayerSemiImplicit::stable_dt(model, MAX_WIND, 0.6);
    DAY / (DAY / dt).ceil()
}

/// A model with the default (strongest stable) ∇⁴ hyperdiffusion, and its
/// step.
fn diffused_model(mountain: bool) -> (TwoLayer<'static>, f64) {
    let model = bare_model(mountain);
    let dt = step_for(&model);
    let options = HyperdiffusionOptions::strongest_for_step(model.operators(), 4, dt);
    (model.with_hyperdiffusion(options), dt)
}

/// The storm setup of `fms baroclinic`: a 30 m/s-shear jet seeded with
/// 0.5 K waves 4–8 at 45° N and S.
fn storm_state(model: &TwoLayer) -> TwoLayerState {
    let flow = ShearedFlow {
        speed: [40.0, 10.0],
        stability: 40.0,
        theta_equator: 315.0,
        ps_equator: 1.0e5,
        omega: OMEGA,
    };
    let mut state = flow.state(model);
    for k in 0..2 {
        for (i, &p) in mesh().0.cell_centers().iter().enumerate() {
            let (lat, lon) = to_lat_lon(p);
            let band =
                |c: f64| fms_math::exp(-((lat - c.to_radians()) / 15f64.to_radians()).powi(2));
            let waves: f64 = (4..=8)
                .map(|m| fms_math::cos(f64::from(m) * lon + 0.7 * f64::from(m)))
                .sum();
            state.theta_ps[k][i] += state.ps[i] * 0.5 * (band(45.0) + band(-45.0)) * waves;
        }
    }
    state
}

fn stepper(
    model: TwoLayer<'static>,
    state: TwoLayerState,
    dt: f64,
) -> TwoLayerSemiImplicit<'static> {
    TwoLayerSemiImplicit::new(
        model,
        state,
        dt,
        TwoLayerSemiImplicitOptions::default(),
        0.5,
    )
}

fn max_abs(x: &[f64]) -> f64 {
    x.iter().fold(0.0f64, |a, v| a.max(v.abs()))
}

#[test]
fn diffusion_conserves_theta_and_removes_kinetic_energy() {
    let (mut model, _) = diffused_model(true);
    let mut state = storm_state(&model);
    // Add grid-scale noise so every part of the operator is exercised.
    let mut s = 12345u64;
    let mut noise = || {
        s ^= s << 13;
        s ^= s >> 7;
        s ^= s << 17;
        (s >> 11) as f64 / (1u64 << 52) as f64 - 1.0
    };
    for k in 0..2 {
        for u in &mut state.u[k] {
            *u += 5.0 * noise();
        }
        for (t, p) in state.theta_ps[k].iter_mut().zip(&state.ps) {
            *t += p * 2.0 * noise();
        }
    }
    let mut d = state.clone();
    model.diffusion_tendency(&state, &mut d);
    assert!(d.ps.iter().all(|&x| x == 0.0), "surface pressure diffused");
    let areas = &mesh().1.cell_areas;
    for k in 0..2 {
        let net: f64 = d.theta_ps[k].iter().zip(areas).map(|(x, a)| x * a).sum();
        let scale: f64 = d.theta_ps[k]
            .iter()
            .zip(areas)
            .map(|(x, a)| (x * a).abs())
            .sum();
        assert!(net.abs() < 1e-13 * scale, "layer {k}: net Θ source {net:e}");
    }
    // The wind part alone only removes energy.
    let winds_only = TwoLayerState {
        theta_ps: [vec![0.0; d.ps.len()], vec![0.0; d.ps.len()]],
        ..d.clone()
    };
    let rate = model.energy_tendency(&state, &winds_only).total();
    eprintln!("kinetic energy dissipation {rate:.3e} W");
    assert!(rate < 0.0);
}

#[test]
fn isentropic_rest_over_a_mountain_is_untouched() {
    // Uniform θ and no wind: diffusion acts on θ (not Θ = p_s θ), so the
    // varying surface pressure over the mountain doesn't stir anything.
    let (model, dt) = diffused_model(true);
    let theta0 = 300.0;
    let exner_ref = GAS.exner(1.0e5);
    let ps: Vec<f64> = model
        .surface_geopotential()
        .iter()
        .map(|phi| {
            1.0e5
                * fms_math::pow(
                    exner_ref - phi / (GAS.heat_capacity * theta0),
                    1.0 / GAS.kappa(),
                )
        })
        .collect();
    let (nc, ne) = (ps.len(), mesh().0.n_edges());
    let state = TwoLayerState::new(
        ps,
        [vec![0.0; ne], vec![0.0; ne]],
        [vec![theta0; nc], vec![theta0; nc]],
    );
    let mut si = stepper(model, state, dt);
    while si.steps() < (5.0 * DAY / dt) as u64 {
        si.step();
    }
    let wind = si.state().u.iter().map(|u| max_abs(u)).fold(0.0, f64::max);
    eprintln!("isentropic rest with diffusion, day 5: max |u| {wind:.1e} m/s");
    assert!(wind < 1e-3, "spurious wind {wind} m/s");
}

/// Run the storm for `days`, checking every step. Returns the stepper, the
/// peak wind, the largest lower- and upper-layer |ζ| seen, and the initial
/// diagnostics; or the day it blew up.
fn storm(
    model: TwoLayer<'static>,
    dt: f64,
    days: f64,
) -> Result<
    (
        TwoLayerSemiImplicit<'static>,
        f64,
        f64,
        fms_atmos::TwoLayerDiagnostics,
    ),
    f64,
> {
    let state = storm_state(&model);
    let start = model.diagnostics(&state);
    let mut si = stepper(model, state, dt);
    let per_day = (DAY / dt).round() as u64;
    let (mut peak, mut zeta) = (0.0f64, 0.0f64);
    while si.steps() < (days * per_day as f64) as u64 {
        si.step();
        match si.check() {
            Ok(h) => peak = peak.max(h.max_wind),
            Err(_) => return Err(si.steps() as f64 / per_day as f64),
        }
        if si.steps() % (per_day / 2) == 0 {
            for k in 0..2 {
                zeta = zeta.max(max_abs(&si.model().vorticity(si.state(), k)));
            }
        }
    }
    Ok((si, peak, zeta, start))
}

#[test]
fn a_storm_life_cycle_survives_a_month_with_closed_budgets() {
    let (model, dt) = diffused_model(false);
    let tau = model.hyperdiffusion().unwrap().options().e_folding;
    let (si, peak, zeta, start) = storm(model, dt, 30.0).expect("the diffused storm blew up");
    let end = si.model().diagnostics(si.state());
    let rel = |a: f64, b: f64| b / a - 1.0;
    let dissipated = si.dissipated_energy();
    let residual = end.energy - start.energy - dissipated;
    eprintln!(
        "∇⁴, τ = {:.1} h: peak |u| {peak:.1} m/s, max |ζ| {zeta:.1e} s⁻¹; mass {:+.1e} Θ {:+.1e}; \
         energy {:+.2e} J = dissipated {:+.2e} + residual {:+.2e}",
        tau / 3600.0,
        rel(start.mass, end.mass),
        rel(start.theta, end.theta),
        end.energy - start.energy,
        dissipated,
        residual,
    );
    // Measured: 90 m/s, 2.4e-4 s⁻¹ (inviscid: 129 m/s, 6e-4 s⁻¹ and rising).
    assert!(peak < 100.0, "peak wind {peak} m/s");
    assert!(zeta < 3e-4, "vorticity reached {zeta} s⁻¹");
    assert!(rel(start.mass, end.mass).abs() < 1e-13);
    assert!(rel(start.theta, end.theta).abs() < 1e-13);
    // Dissipation removed ~4e-4 of the total energy; time-stepping error
    // accounts for the rest (measured: 0.5% of the dissipation).
    assert!(dissipated < -1e-5 * start.energy, "too little dissipation");
    assert!(
        residual.abs() < 0.02 * dissipated.abs(),
        "budget residual {residual:e} J"
    );
}

#[test]
fn without_diffusion_the_storm_collapses_to_the_grid_scale() {
    // The control for the test above, at the same step: the inviscid core
    // survives the month only because the step was sized for fast winds,
    // and its vorticity keeps sharpening.
    let model = bare_model(false);
    let dt = step_for(&model);
    let (_, peak, zeta, _) = match storm(model, dt, 30.0) {
        Ok(run) => run,
        Err(day) => {
            eprintln!("inviscid storm blew up on day {day:.1}");
            return;
        }
    };
    eprintln!("inviscid: peak |u| {peak:.1} m/s, max |ζ| {zeta:.1e} s⁻¹");
    assert!(zeta > 4e-4, "vorticity only reached {zeta} s⁻¹");
}

#[test]
fn an_injected_nan_is_caught_within_a_step() {
    let (model, dt) = diffused_model(false);
    let state = storm_state(&model);
    let mut si = stepper(model, state, dt);
    si.step();
    assert!(si.check().is_ok());
    si.state_mut().theta_ps[1][100] = f64::NAN;
    si.step();
    let blow_up = si.check().unwrap_err();
    eprintln!("{blow_up}");
    assert_eq!(blow_up.fault, Fault::NotFinite);
}

#[test]
fn a_wind_beyond_the_courant_limit_is_reported_where_it_is() {
    let (model, dt) = diffused_model(false);
    let state = storm_state(&model);
    let mut si = stepper(model, state, dt);
    si.step();
    let limit = si.courant_limit();
    si.state_mut().u[0][500] = -1.5 * limit;
    let blow_up = si.check().unwrap_err();
    eprintln!("{blow_up}");
    assert_eq!(
        (blow_up.variable, blow_up.fault, blow_up.index),
        (Variable::Wind(0), Fault::TooFast, 500)
    );
    assert_eq!((blow_up.value, blow_up.limit), (-1.5 * limit, limit));
}

#[test]
fn a_runaway_is_caught_before_it_turns_into_nans() {
    // Inviscid, with the step sized for 45 m/s: the storms outrun it.
    let model = bare_model(false);
    let dt = TwoLayerSemiImplicit::stable_dt(&model, 45.0, 0.6);
    let dt = DAY / (DAY / dt).ceil();
    let state = storm_state(&model);
    let mut si = stepper(model, state, dt);
    let fault = loop {
        si.step();
        if let Err(b) = si.check() {
            eprintln!("day {:.2}: {b}", si.steps() as f64 * dt / DAY);
            break b.fault;
        }
        assert!(si.steps() as f64 * dt < 30.0 * DAY, "no blow-up in 30 days");
    };
    assert_eq!(fault, Fault::TooFast);
}

#[test]
fn out_of_range_values_are_reported() {
    let model = bare_model(false);
    let mut state = storm_state(&model);
    let limits = Limits::around(200.0, 1.0e5, [330.0, 290.0]);
    assert!(safety::check(&state, &limits).is_ok());
    state.ps[7] = 2.0e4;
    let blow_up = safety::check(&state, &limits).unwrap_err();
    assert_eq!(
        (blow_up.variable, blow_up.fault, blow_up.index),
        (Variable::SurfacePressure, Fault::OutOfRange, 7)
    );
}
