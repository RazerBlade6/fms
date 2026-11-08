//! The two-layer model stepped semi-implicitly against RK4: balance,
//! conservation, rest over a mountain, baroclinic growth, cost.
//!
//! ```sh
//! cargo run --release -p fms-atmos --example two_layer_si [frequency]
//! ```

use std::time::Instant;

use fms_atmos::williamson::{Case, GRAVITY, OMEGA, RADIUS};
use fms_atmos::{
    Gas, ShearedFlow, TwoLayer, TwoLayerSemiImplicit, TwoLayerSemiImplicitOptions, TwoLayerState,
};
use fms_geo::coords::to_lat_lon;
use fms_geo::{LloydOptions, MeshGeometry, Operators, SphereMesh};

const DAY: f64 = 86_400.0;
/// Design maximum wind for choosing the semi-implicit step (m/s).
const MAX_WIND: f64 = 60.0;
const COURANT: f64 = 0.6;

struct Setup {
    mesh: SphereMesh,
    geometry: MeshGeometry,
}

impl Setup {
    fn model(&self, mountain: bool) -> TwoLayer<'_> {
        let phi_s = if mountain {
            self.mesh
                .cell_centers()
                .iter()
                .map(|&p| GRAVITY * Case::Mountain.topography(p))
                .collect()
        } else {
            vec![0.0; self.mesh.n_cells()]
        };
        TwoLayer::new(
            Operators::new(&self.mesh, &self.geometry, RADIUS),
            GRAVITY,
            Gas::EARTH_AIR,
            0.5,
            self.mesh
                .vertices()
                .iter()
                .map(|p| 2.0 * OMEGA * p.z)
                .collect(),
            phi_s,
        )
    }
}

/// Either scheme behind one interface.
enum Stepper<'a> {
    /// Model, state, step, steps taken.
    Rk4(Box<TwoLayer<'a>>, TwoLayerState, f64, u64),
    Si(Box<TwoLayerSemiImplicit<'a>>),
}

impl Stepper<'_> {
    fn new(model: TwoLayer<'_>, state: TwoLayerState, semi_implicit: bool) -> Stepper<'_> {
        if semi_implicit {
            let dt = TwoLayerSemiImplicit::stable_dt(&model, MAX_WIND, COURANT);
            let dt = DAY / (DAY / dt).ceil();
            let si = TwoLayerSemiImplicit::new(
                model,
                state,
                dt,
                TwoLayerSemiImplicitOptions::default(),
                0.5,
            );
            Stepper::Si(Box::new(si))
        } else {
            let dt = model.stable_dt(&state, 0.5);
            let dt = DAY / (DAY / dt).ceil();
            Stepper::Rk4(Box::new(model), state, dt, 0)
        }
    }
    fn dt(&self) -> f64 {
        match self {
            Self::Rk4(_, _, dt, _) => *dt,
            Self::Si(si) => si.dt(),
        }
    }
    /// Advance to the end of `days` from the start.
    fn run_days(&mut self, days: f64) {
        let target = (days * DAY / self.dt()).round() as u64;
        match self {
            Self::Rk4(m, s, dt, taken) => {
                while *taken < target {
                    m.step(s, *dt);
                    *taken += 1;
                }
            }
            Self::Si(si) => {
                while si.steps() < target {
                    si.step();
                }
            }
        }
    }
    fn state(&self) -> &TwoLayerState {
        match self {
            Self::Rk4(_, s, _, _) => s,
            Self::Si(si) => si.state(),
        }
    }
    fn model(&self) -> &TwoLayer<'_> {
        match self {
            Self::Rk4(m, _, _, _) => m,
            Self::Si(si) => si.model(),
        }
    }
}

fn main() {
    let k: u32 = std::env::args().nth(1).map_or(16, |s| s.parse().unwrap());
    let (mesh, _) = SphereMesh::centroidal(k, LloydOptions::default());
    let setup = Setup {
        geometry: MeshGeometry::new(&mesh),
        mesh,
    };
    println!("frequency {k}: {} cells", setup.mesh.n_cells());
    let flow = |shear: f64| ShearedFlow {
        speed: [10.0 + shear, 10.0],
        stability: 40.0,
        theta_equator: 315.0,
        ps_equator: 1.0e5,
        omega: OMEGA,
    };

    for si in [false, true] {
        let name = if si { "semi-implicit" } else { "RK4" };
        // Sheared flow over the mountain, 10 days: conservation, cost.
        let model = setup.model(true);
        let state = flow(30.0).state(&model);
        let start = model.diagnostics(&state);
        let explicit_dt = model.stable_dt(&state, 0.5);
        let mut stepper = Stepper::new(model, state, si);
        let t = Instant::now();
        stepper.run_days(10.0);
        let elapsed = t.elapsed().as_secs_f64();
        let end = stepper.model().diagnostics(stepper.state());
        let max_u = stepper
            .state()
            .u
            .iter()
            .flatten()
            .fold(0.0f64, |a, &u| a.max(u.abs()));
        print!(
            "{name:<14} dt {:>5.0} s ({:.1}× explicit)  10 days in {elapsed:.2}s  mass {:+.1e} θ {:+.1e} energy {:+.1e}  max |u| {max_u:.0}",
            stepper.dt(),
            stepper.dt() / explicit_dt,
            end.mass / start.mass - 1.0,
            end.theta / start.theta - 1.0,
            end.energy / start.energy - 1.0,
        );
        if let Stepper::Si(s) = &stepper {
            let st = s.solver_stats();
            print!(
                "  CG {:.1}/solve  c = {:.0?} m/s",
                st.iterations as f64 / st.solves as f64,
                s.wave_speeds()
            );
        }
        println!();

        // Isothermal rest over the mountain, 5 days.
        let model = setup.model(true);
        let t0 = 260.0;
        let ps: Vec<f64> = model
            .surface_geopotential()
            .iter()
            .map(|phi| 1.0e5 * fms_math::exp(-phi / (Gas::EARTH_AIR.gas_constant * t0)))
            .collect();
        let theta = {
            let mut th = [vec![0.0; ps.len()], vec![0.0; ps.len()]];
            for (i, &p) in ps.iter().enumerate() {
                let e = model.layer_exner(p);
                th[0][i] = t0 / e[0];
                th[1][i] = t0 / e[1];
            }
            th
        };
        let ne = setup.mesh.n_edges();
        let state = TwoLayerState::new(ps, [vec![0.0; ne], vec![0.0; ne]], theta);
        let mut stepper = Stepper::new(model, state, si);
        stepper.run_days(5.0);
        let spurious = stepper
            .state()
            .u
            .iter()
            .flatten()
            .fold(0.0f64, |a, &u| a.max(u.abs()));
        println!("{name:<14} isothermal rest over the cone, day 5: max |u| {spurious:.3} m/s");

        // Baroclinic growth of m = 6 at 45 m/s shear, against a control.
        let f = flow(45.0);
        let reference = setup.model(false);
        let (theory, _) = f.phillips_growth_rate(&reference, 6, 45f64.to_radians());
        let make = |seed: f64| {
            let model = setup.model(false);
            let mut state = f.state(&model);
            for k in 0..2 {
                for (i, &p) in setup.mesh.cell_centers().iter().enumerate() {
                    let (lat, lon) = to_lat_lon(p);
                    let bump =
                        fms_math::exp(-((lat - 45f64.to_radians()) / 15f64.to_radians()).powi(2));
                    state.theta_ps[k][i] += state.ps[i] * seed * bump * fms_math::cos(6.0 * lon);
                }
            }
            Stepper::new(model, state, si)
        };
        let (mut control, mut perturbed) = (make(0.0), make(0.01));
        let mut samples = Vec::new();
        for day in 0..=8 {
            if day > 0 {
                control.run_days(f64::from(day));
                perturbed.run_days(f64::from(day));
            }
            let a = wave6(&setup, perturbed.state(), control.state());
            samples.push((f64::from(day), fms_math::ln(a)));
        }
        let fit = &samples[2..];
        let n = fit.len() as f64;
        let (mt, my) = fit
            .iter()
            .fold((0.0, 0.0), |(a, b), (t, y)| (a + t / n, b + y / n));
        let (num, den) = fit.iter().fold((0.0, 0.0), |(a, b), (t, y)| {
            (a + (t - mt) * (y - my), b + (t - mt).powi(2))
        });
        let growth = num / den;
        println!(
            "{name:<14} baroclinic m=6 at 45 m/s: growth {growth:.3} /day (theory {:.3}, ratio {:.2})",
            theory * DAY,
            growth / (theory * DAY)
        );
    }
}

fn wave6(setup: &Setup, a: &TwoLayerState, b: &TwoLayerState) -> f64 {
    let (mut re, mut im, mut area) = (0.0, 0.0, 0.0);
    for (i, &p) in setup.mesh.cell_centers().iter().enumerate() {
        let (lat, lon) = to_lat_lon(p);
        if (30f64.to_radians()..60f64.to_radians()).contains(&lat) {
            let th = |s: &TwoLayerState| 0.5 * (s.theta_ps[0][i] + s.theta_ps[1][i]) / s.ps[i];
            let (s, c) = fms_math::sin_cos(6.0 * lon);
            let d = setup.geometry.cell_areas[i] * (th(a) - th(b));
            re += d * c;
            im += d * s;
            area += setup.geometry.cell_areas[i];
        }
    }
    2.0 * (re * re + im * im).sqrt() / area
}
