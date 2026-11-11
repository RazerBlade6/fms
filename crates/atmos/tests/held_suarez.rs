//! Held–Suarez benchmark (M2 step 7): the forcing's layer averages, and the
//! climate it drives on the two-layer core — jets in the right place, the
//! two hemispheres agreeing, and a steady, closed energy budget.
//!
//! The climate tests share one run (k = 16, 200 days of spin-up, then 250
//! days of statistics). The full benchmark (1200 days at k = 32) is
//! `fms-climate held-suarez`; see docs/climate/CLI.md.

use fms_atmos::column_physics::CLIMATE_LEVELS;
use std::sync::OnceLock;

use fms_atmos::williamson::{GRAVITY, OMEGA, RADIUS};
use fms_atmos::{
    EnergyBudget, Forcing, Gas, HeldSuarez, HeldSuarezParams, HyperdiffusionOptions,
    LayeredDiagnostics, LayeredModel, LayeredSemiImplicit, LayeredSemiImplicitOptions, Levels,
    ZonalClimatology, ZonalMeans,
};
use fms_grid::{LloydOptions, MeshGeometry, Operators, SphereMesh};

const DAY: f64 = 86_400.0;
const SPINUP: u64 = 200;
const DAYS: u64 = 450;
const BANDS: usize = 36;

fn mesh() -> &'static (SphereMesh, MeshGeometry) {
    static M: OnceLock<(SphereMesh, MeshGeometry)> = OnceLock::new();
    M.get_or_init(|| {
        let (m, _) = SphereMesh::centroidal(16, LloydOptions::default());
        let g = MeshGeometry::new(&m);
        (m, g)
    })
}

fn model() -> LayeredModel<'static> {
    model_with(Levels::two(0.5))
}

fn model_with(levels: Levels) -> LayeredModel<'static> {
    let (m, g) = mesh();
    LayeredModel::new(
        Operators::new(m, g, RADIUS),
        GRAVITY,
        Gas::EARTH_AIR,
        levels,
        m.vertices().iter().map(|p| 2.0 * OMEGA * p.z).collect(),
        vec![0.0; m.n_cells()],
    )
}

#[test]
fn layer_averages_of_the_forcing() {
    let model = model();
    let hs = HeldSuarez::new(&model, HeldSuarezParams::default());
    let p = HeldSuarezParams::default();
    // ⟨max(0, (σ − 0.7)/0.3)⟩ over σ ∈ [0.5, 1] is 0.3.
    let &[upper, lower] = hs.friction() else {
        panic!("two layers")
    };
    assert_eq!(upper, 0.0);
    assert!((lower * DAY - 0.3).abs() < 1e-12, "{}", lower * DAY);
    // Brute-force midpoint averages of k_T and k_T θ_eq for comparison.
    let brute = |sin2: f64, ps: f64, (a, b): (f64, f64)| {
        let n = 200_000;
        let (mut rate, mut forcing) = (0.0, 0.0);
        for j in 0..n {
            let sigma = a + (b - a) * (j as f64 + 0.5) / n as f64;
            let f = ((sigma - p.sigma_boundary) / (1.0 - p.sigma_boundary)).max(0.0);
            let cos2 = 1.0 - sin2;
            let k =
                p.relaxation_free + (p.relaxation_surface - p.relaxation_free) * f * cos2 * cos2;
            rate += k / n as f64;
            forcing += k * hs.theta_eq(sin2, sigma * ps) / n as f64;
        }
        (rate, forcing / rate)
    };
    let centers = mesh().0.cell_centers();
    for cell in [0, 400, 1281, 2000] {
        let sin2 = centers[cell].z * centers[cell].z;
        for ps in [0.9e5, 1.0e5, 1.03e5] {
            let got = hs.equilibrium(cell, ps);
            for (k, range) in [(0.0, 0.5), (0.5, 1.0)].into_iter().enumerate() {
                let (rate, theta) = brute(sin2, ps, range);
                assert!(
                    (got[k].0 / rate - 1.0).abs() < 1e-6,
                    "k_T layer {k}: {} vs {rate}",
                    got[k].0
                );
                // Limited by the kink where the stratospheric floor takes over.
                assert!(
                    (got[k].1 - theta).abs() < 0.02,
                    "θ_eq layer {k} at sin²φ {sin2:.2}, p_s {ps}: {} vs {theta}",
                    got[k].1
                );
            }
        }
    }
}

#[test]
fn the_forcing_vanishes_at_equilibrium_and_at_rest() {
    let mut model = model();
    let mut hs = HeldSuarez::new(&model, HeldSuarezParams::default());
    let mut state = hs.initial_state(&model, 0.0);
    for i in 0..state.ps.len() {
        let eq = hs.equilibrium(i, state.ps[i]);
        for (theta, (_, eq)) in state.theta_ps.iter_mut().zip(eq) {
            theta[i] = state.ps[i] * eq;
        }
    }
    let mut out = state.clone();
    model.forcing_tendency(&state, &mut out);
    let worst = out
        .theta_ps
        .iter()
        .flatten()
        .zip(state.theta_ps.iter().flatten())
        .fold(0.0f64, |a, (d, t)| a.max((d / t).abs()));
    assert!(worst < 1e-15, "relative Θ tendency {worst:e} s⁻¹");
    assert!(out.u.iter().flatten().all(|&d| d == 0.0));
    // And friction is a pure linear drag on the lower layer.
    state.u[1].iter_mut().for_each(|u| *u = 10.0);
    out.u[1].fill(0.0);
    hs.add_tendency(&state, &mut out);
    assert!(out.u[1].iter().all(|&d| (d + 3.0 / DAY).abs() < 1e-18));
}

struct Climate {
    means: ZonalMeans,
    /// Diagnostics and budget at the end of spin-up and at the end.
    spun_up: (LayeredDiagnostics, EnergyBudget),
    end: (LayeredDiagnostics, EnergyBudget),
    peak_wind: f64,
}

fn climate() -> &'static Climate {
    static C: OnceLock<Climate> = OnceLock::new();
    C.get_or_init(|| {
        let model = model();
        let dt = LayeredSemiImplicit::stable_dt(&model, 80.0, 0.6);
        let per_day = (DAY / dt).ceil() as u64;
        let dt = DAY / per_day as f64;
        let diffusion = HyperdiffusionOptions::strongest_for_step(model.operators(), 4, dt);
        let hs = HeldSuarez::new(&model, HeldSuarezParams::default());
        let state = hs.initial_state(&model, 0.1);
        let model = model
            .with_hyperdiffusion(diffusion)
            .with_forcing(Box::new(hs));
        let mut zonal = ZonalClimatology::new(&model, BANDS);
        let mut si =
            LayeredSemiImplicit::new(model, state, dt, LayeredSemiImplicitOptions::default(), 0.5);
        let mut spun_up = None;
        let mut peak_wind = 0.0f64;
        while si.steps() < DAYS * per_day {
            si.step();
            let health = si.check().expect("the Held–Suarez run blew up");
            peak_wind = peak_wind.max(health.max_wind);
            if si.steps() == SPINUP * per_day {
                spun_up = Some((si.model().diagnostics(si.state()), si.energy_budget()));
            }
            if si.steps() > SPINUP * per_day && si.steps() % (per_day / 4) == 0 {
                zonal.sample(si.model(), si.state());
            }
        }
        Climate {
            means: zonal.means(),
            spun_up: spun_up.unwrap(),
            end: (si.model().diagnostics(si.state()), si.energy_budget()),
            peak_wind,
        }
    })
}

/// Index and value of the maximum of `v` over bands whose |latitude| is in
/// `[lo, hi]` degrees, in one hemisphere.
fn max_in(m: &ZonalMeans, v: &[f64], north: bool, (lo, hi): (f64, f64)) -> (usize, f64) {
    (0..v.len())
        .filter(|&i| {
            let lat = m.latitude[i];
            (lat > 0.0) == north && (lo..=hi).contains(&lat.abs())
        })
        .map(|i| (i, v[i]))
        .fold((0, f64::NEG_INFINITY), |a, b| if b.1 > a.1 { b } else { a })
}

fn band_mean(m: &ZonalMeans, v: &[f64], (lo, hi): (f64, f64)) -> f64 {
    let picked: Vec<f64> = (0..v.len())
        .filter(|&i| (lo..=hi).contains(&m.latitude[i].abs()))
        .map(|i| v[i])
        .collect();
    picked.iter().sum::<f64>() / picked.len() as f64
}

#[test]
fn jets_and_temperatures_are_in_the_right_place() {
    let c = climate();
    let m = &c.means;
    eprintln!("peak wind {:.1} m/s", c.peak_wind);
    for north in [true, false] {
        let (i, jet) = max_in(m, &m.wind[0], north, (20.0, 70.0));
        let (j, low) = max_in(m, &m.wind[1], north, (20.0, 70.0));
        eprintln!(
            "{}: upper jet {jet:.1} m/s at {:.1}°, lower westerlies {low:.1} m/s at {:.1}°",
            if north { "north" } else { "south" },
            m.latitude[i],
            m.latitude[j]
        );
        // Measured: 25–28 m/s at 47.5–52.5°; 5–6 m/s at 52.5°.
        assert!((20.0..40.0).contains(&jet), "upper jet {jet} m/s");
        assert!(
            (35.0..=55.0).contains(&m.latitude[i].abs()),
            "jet at {}°",
            m.latitude[i]
        );
        assert!(low > 3.0, "lower-layer westerlies only {low} m/s");
        assert!(
            (35.0..=60.0).contains(&m.latitude[j].abs()),
            "westerlies at {}°",
            m.latitude[j]
        );
    }
    // Easterlies in the lower layer's tropics and subtropics (measured −1.5).
    let tropics = band_mean(m, &m.wind[1], (0.0, 30.0));
    assert!(tropics < -0.5, "tropical lower-layer wind {tropics} m/s");
    // Equator-to-pole contrast of the lower layer: HS's T_eq has 60 K, the
    // eddies carry heat poleward and cut it (measured 43 K).
    let contrast = band_mean(m, &m.temperature[1], (0.0, 10.0))
        - band_mean(m, &m.temperature[1], (80.0, 90.0));
    eprintln!("tropical lower-layer wind {tropics:.2} m/s; temperature contrast {contrast:.1} K");
    assert!((30.0..55.0).contains(&contrast), "contrast {contrast} K");
}

#[test]
fn the_hemispheres_agree() {
    // Symmetric forcing on a nearly symmetric mesh: the time means of the two
    // hemispheres differ only by sampling noise (measured: 1.5 m/s, 0.5 K),
    // and the broad jets peak within a band of each other (47.5° and
    // 52.5° N against 47.5° S, at 5° bands).
    let m = &climate().means;
    let mirror = m.mirrored();
    let worst = |a: &[f64], b: &[f64]| {
        a.iter()
            .zip(b)
            .fold(0.0f64, |w, (x, y)| w.max((x - y).abs()))
    };
    let du = worst(&m.wind[0], &mirror.wind[0]);
    let dl = worst(&m.wind[1], &mirror.wind[1]);
    let dt = worst(&m.temperature[1], &mirror.temperature[1]);
    eprintln!(
        "hemispheric differences: upper wind {du:.2} m/s, lower {dl:.2} m/s, lower T {dt:.2} K"
    );
    assert!(du < 4.0 && dl < 2.0 && dt < 1.5);
    let (n, _) = max_in(m, &m.wind[0], true, (20.0, 70.0));
    let (s, _) = max_in(m, &m.wind[0], false, (20.0, 70.0));
    let mirror_band = m.latitude.len() - 1 - s;
    assert!(n.abs_diff(mirror_band) <= 1, "jets at different latitudes");
}

#[test]
fn the_energy_budget_is_steady_and_closed() {
    let c = climate();
    let ((d0, b0), (d1, b1)) = (c.spun_up, c.end);
    let seconds = (DAYS - SPINUP) as f64 * DAY;
    let area = 4.0 * std::f64::consts::PI * RADIUS * RADIUS;
    let wm2 = |j: f64| j / seconds / area;
    let change = d1.energy - d0.energy;
    let heating = b1.forcing.thermal - b0.forcing.thermal;
    let friction = b1.forcing.kinetic - b0.forcing.kinetic;
    let diffusion = b1.diffusion.total() - b0.diffusion.total();
    let residual = change - heating - friction - diffusion;
    eprintln!(
        "after spin-up (W/m²): dE/dt {:+.3} = Newtonian {:+.3} + friction {:+.3} + diffusion {:+.3} + residual {:+.4}",
        wm2(change),
        wm2(heating),
        wm2(friction),
        wm2(diffusion),
        wm2(residual)
    );
    // Steady: what Newtonian heating puts in, friction and diffusion take
    // out (measured: the net is 6% of the input, sampling noise).
    assert!(heating > 0.0 && friction < 0.0 && diffusion < 0.0);
    assert!(change.abs() < 0.15 * heating, "energy drifting");
    // Closed: the residual is time-stepping error (measured 0.1%).
    assert!(residual.abs() < 0.01 * heating, "residual {residual:e} J");
    assert!((d1.mass / d0.mass - 1.0).abs() < 1e-13, "mass drifted");
}

#[test]
fn more_layers_spin_up_westerly_jets() {
    // Held–Suarez on three layers (a thin bottom one, σ 0.85–1), the
    // climate's four and five: 200 days at the long step without blowing
    // up, mass conserved, and westerly jets in both hemispheres at the layer
    // holding 250 hPa. On the climate's layers (O-18) the lowest layer's
    // winds under the jets are westerly too.
    const DAYS: u64 = 200;
    for levels in [
        Levels::new(&[0.4, 0.85]),
        Levels::new(&CLIMATE_LEVELS),
        Levels::equal_mass(5),
    ] {
        let climate_layers = levels.interior() == CLIMATE_LEVELS;
        let n = levels.count();
        let jet_layer = levels.layer_of(0.25);
        let model = model_with(levels);
        let dt = LayeredSemiImplicit::stable_dt(&model, 80.0, 0.6);
        let per_day = (DAY / dt).ceil() as u64;
        let dt = DAY / per_day as f64;
        let diffusion = HyperdiffusionOptions::strongest_for_step(model.operators(), 4, dt);
        let hs = HeldSuarez::new(&model, HeldSuarezParams::default());
        let state = hs.initial_state(&model, 0.1);
        let model = model
            .with_hyperdiffusion(diffusion)
            .with_forcing(Box::new(hs));
        let mut zonal = ZonalClimatology::new(&model, BANDS);
        let start = model.diagnostics(&state);
        let mut si =
            LayeredSemiImplicit::new(model, state, dt, LayeredSemiImplicitOptions::default(), 0.5);
        let mut peak = 0.0f64;
        while si.steps() < DAYS * per_day {
            si.step();
            peak = peak.max(
                si.check()
                    .expect("the layered Held–Suarez run blew up")
                    .max_wind,
            );
            if si.steps() > DAYS / 2 * per_day && si.steps() % (per_day / 4) == 0 {
                zonal.sample(si.model(), si.state());
            }
        }
        let end = si.model().diagnostics(si.state());
        let m = zonal.means();
        let jets = [true, false].map(|north| max_in(&m, &m.wind[jet_layer], north, (20.0, 70.0)));
        let surface = [true, false].map(|north| {
            let (i, _) = max_in(&m, &m.wind[jet_layer], north, (20.0, 70.0));
            m.wind[n - 1][i]
        });
        eprintln!(
            "{n} layers: jets {:.1} m/s at {:.0}° N, {:.1} m/s at {:.0}° S (layer {jet_layer}); \
             lowest layer under them {:+.1}, {:+.1} m/s; peak wind {peak:.0} m/s; mass {:+.1e}",
            jets[0].1,
            m.latitude[jets[0].0],
            jets[1].1,
            m.latitude[jets[1].0],
            surface[0],
            surface[1],
            end.mass / start.mass - 1.0
        );
        assert!((end.mass / start.mass - 1.0).abs() < 1e-12);
        for (_, jet) in jets {
            assert!(jet > 15.0, "{n} layers: jet only {jet} m/s");
        }
        if climate_layers {
            for wind in surface {
                assert!(wind > 0.0, "{n} layers: {wind} m/s under the jet");
            }
        }
    }
}
