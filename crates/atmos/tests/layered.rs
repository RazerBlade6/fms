//! Layered dry core: conservation and rest over mountains for several layer
//! sets, the vertical modes, and (two layers) baroclinic instability
//! against linear theory.

// The energy-budget test indexes cells, edges and layers explicitly so each
// sum reads like the formula it checks.
#![allow(clippy::needless_range_loop)]

use std::sync::OnceLock;

use fms_atmos::williamson::{Case, GRAVITY, OMEGA, RADIUS};
use fms_atmos::{Gas, LayeredModel, LayeredState, Levels, ShearedFlow};
use fms_geo::coords::to_lat_lon;
use fms_geo::{LloydOptions, MeshGeometry, Operators, SphereMesh};

const GAS: Gas = Gas::EARTH_AIR;

fn mesh(k: u32) -> &'static (SphereMesh, MeshGeometry) {
    static M: [OnceLock<(SphereMesh, MeshGeometry)>; 2] = [OnceLock::new(), OnceLock::new()];
    let slot = match k {
        16 => &M[0],
        32 => &M[1],
        _ => panic!("unsupported test frequency"),
    };
    slot.get_or_init(|| {
        let (m, _) = SphereMesh::centroidal(k, LloydOptions::default());
        let g = MeshGeometry::new(&m);
        (m, g)
    })
}

fn noise(n: usize, seed: u64) -> Vec<f64> {
    let mut x = seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1;
    (0..n)
        .map(|_| {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            (x >> 11) as f64 / (1u64 << 52) as f64 - 1.0
        })
        .collect()
}

fn earth_coriolis(m: &SphereMesh) -> Vec<f64> {
    m.vertices().iter().map(|p| 2.0 * OMEGA * p.z).collect()
}

/// Surface geopotential of the Williamson case-5 cone (2 km high).
fn mountain(m: &SphereMesh) -> Vec<f64> {
    m.cell_centers()
        .iter()
        .map(|&p| GRAVITY * Case::Mountain.topography(p))
        .collect()
}

fn model<'a>(m: &'a SphereMesh, g: &'a MeshGeometry, phi_s: Vec<f64>) -> LayeredModel<'a> {
    model_with(m, g, phi_s, Levels::two(0.5))
}

fn model_with<'a>(
    m: &'a SphereMesh,
    g: &'a MeshGeometry,
    phi_s: Vec<f64>,
    levels: Levels,
) -> LayeredModel<'a> {
    LayeredModel::new(
        Operators::new(m, g, RADIUS),
        GRAVITY,
        GAS,
        levels,
        earth_coriolis(m),
        phi_s,
    )
}

/// Two equal layers, three uneven ones with a thin bottom layer, and five.
fn layer_sets() -> [Levels; 3] {
    [
        Levels::two(0.5),
        Levels::new(&[0.4, 0.85]),
        Levels::equal_mass(5),
    ]
}

/// θ decreasing linearly with σ from `top` in the top layer's centre to
/// `bottom` in the lowest's, plus `spread` K of noise (seeded by `seed`).
fn stratified(
    levels: &Levels,
    n: usize,
    top: f64,
    bottom: f64,
    spread: f64,
    seed: u64,
) -> Vec<Vec<f64>> {
    let centre = |k: usize| {
        let (a, b) = levels.bounds(k);
        0.5 * (a + b)
    };
    let (c0, c1) = (centre(0), centre(levels.count() - 1));
    (0..levels.count())
        .map(|k| {
            let x = if c1 > c0 {
                (centre(k) - c0) / (c1 - c0)
            } else {
                0.0
            };
            let mean = top + (bottom - top) * x;
            noise(n, seed + k as u64)
                .iter()
                .map(|x| mean + spread * x)
                .collect()
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Conservation
// ---------------------------------------------------------------------------

#[test]
fn semi_discrete_energy_is_conserved() {
    // For a random state over random terrain, the spatial tendencies change
    // total energy E = Σ Aᵢ [Σₖ δpₖ (Kₖ + c_p θₖ Pₖ) + p_s Φ_s] / g by zero:
    //   g dE/dt = Σₑ lₑdₑ Σₖ δpₖ,ₑ uₖ ∂uₖ
    //           + Σᵢ Aᵢ [Σₖ (Kₖ δσₖ ∂p_s + c_p Pₖ δσₖ ∂Θₖ + c_p κ δpₖ θₖ Pₖ ∂p_s / p_s)
    //                    + Φ_s ∂p_s]
    let (m, g) = mesh(16);
    for levels in layer_sets() {
        let n = levels.count();
        let phi_s: Vec<f64> = noise(m.n_cells(), 1)
            .iter()
            .map(|x| 1.0e4 * (1.0 + x))
            .collect();
        let mut model = model_with(m, g, phi_s.clone(), levels);
        let ps: Vec<f64> = noise(m.n_cells(), 2)
            .iter()
            .map(|x| 9.5e4 + 5.0e3 * x)
            .collect();
        let theta = stratified(model.levels(), m.n_cells(), 330.0, 290.0, 10.0, 3);
        let u = (0..n)
            .map(|k| {
                noise(m.n_edges(), 20 + k as u64)
                    .iter()
                    .map(|x| (30.0 - 15.0 * k as f64 / n as f64) * x)
                    .collect()
            })
            .collect();
        let state = LayeredState::new(ps, u, theta);
        let mut d = state.clone();
        model.tendency(&state, &mut d);

        let ops = Operators::new(m, g, RADIUS);
        let r2 = RADIUS * RADIUS;
        let ds = model.dsigma().to_vec();
        let cp = GAS.heat_capacity;
        let kappa = GAS.kappa();
        let mut terms = Vec::new();
        let mut ps_edge = vec![0.0; m.n_edges()];
        ops.cell_to_edge(&state.ps, &mut ps_edge);
        for k in 0..n {
            for e in 0..m.n_edges() {
                let dp = ds[k] * ps_edge[e];
                terms.push(
                    g.edge_lengths[e]
                        * g.dual_edge_lengths[e]
                        * r2
                        * dp
                        * state.u[k][e]
                        * d.u[k][e],
                );
            }
        }
        let theta: Vec<Vec<f64>> = (0..n).map(|k| state.theta(k)).collect();
        let mut kinetic = vec![vec![0.0; m.n_cells()]; n];
        for k in 0..n {
            ops.kinetic_energy(&state.u[k], &mut kinetic[k]);
        }
        for i in 0..m.n_cells() {
            let a = g.cell_areas[i] * r2;
            let ps = state.ps[i];
            let exner = model.layer_exner(ps);
            terms.push(a * phi_s[i] * d.ps[i]);
            for k in 0..n {
                terms.push(a * kinetic[k][i] * ds[k] * d.ps[i]);
                terms.push(a * cp * exner[k] * ds[k] * d.theta_ps[k][i]);
                terms.push(a * cp * kappa * ds[k] * theta[k][i] * exner[k] * d.ps[i]);
            }
        }
        let total: f64 = terms.iter().sum();
        let scale: f64 = terms.iter().map(|t| t.abs()).sum();
        assert!(
            total.abs() < 1e-12 * scale,
            "{n} layers: dE/dt = {total:e} of {scale:e}"
        );
    }
}

#[test]
fn mass_and_theta_are_conserved_and_energy_drifts_little() {
    // A sheared, baroclinic flow over a mountain, run for 2 days with RK4.
    let (m, g) = mesh(16);
    let mut model = model(m, g, mountain(m));
    let state0 = baroclinic_state(&model, 1.0);
    let start = model.diagnostics(&state0);
    let mut state = state0;
    let dt = model.stable_dt(&state, 0.5);
    let steps = (2.0 * 86_400.0 / dt).ceil() as u32;
    let dt = 2.0 * 86_400.0 / f64::from(steps);
    for _ in 0..steps {
        model.step(&mut state, dt);
    }
    let end = model.diagnostics(&state);
    let rel = |a: f64, b: f64| b / a - 1.0;
    eprintln!(
        "dt {dt:.0} s, {steps} steps: mass {:+.1e}, θ {:+.1e}, energy {:+.1e}",
        rel(start.mass, end.mass),
        rel(start.theta, end.theta),
        rel(start.energy, end.energy)
    );
    assert!(rel(start.mass, end.mass).abs() < 1e-13);
    assert!(rel(start.theta, end.theta).abs() < 1e-13);
    // The spatial scheme conserves energy exactly (test above); what remains
    // is RK4's dissipation, which converges like Δt⁵. Measured: −8.9e-8 at
    // Courant 0.5, −3.4e-9 at 0.25, −1.1e-10 at 0.125.
    assert!(rel(start.energy, end.energy).abs() < 2e-7);
}

#[test]
fn more_layers_conserve_mass_theta_and_water() {
    // A sheared, stratified, moist flow over a mountain on three layers (a
    // thin bottom one) and five, 2 days of RK4: mass, θ and water to
    // round-off, energy to RK4's dissipation. That converges like Δt⁵: at
    // Courant 0.5 it is −3.3e-7 (three layers), and at 0.25, run here,
    // −1.5e-8 (three) and −9.9e-9 (five).
    let (m, g) = mesh(16);
    for levels in [Levels::new(&[0.4, 0.85]), Levels::equal_mass(5)] {
        let n = levels.count();
        let mut model = model_with(m, g, mountain(m), levels);
        let theta = stratified(model.levels(), m.n_cells(), 330.0, 290.0, 2.0, 7);
        let u = (0..n)
            .map(|k| {
                let speed = 10.0 + 30.0 * (n - 1 - k) as f64 / (n - 1) as f64;
                m.edges()
                    .iter()
                    .zip(&g.edge_lengths)
                    .map(|(e, &l)| {
                        let psi = |v: u32| -RADIUS * speed * m.vertices()[v as usize].z;
                        -(psi(e.vertices[1]) - psi(e.vertices[0])) / (l * RADIUS)
                    })
                    .collect()
            })
            .collect();
        let q = (0..n)
            .map(|k| vec![1e-2 * (k + 1) as f64 / n as f64; m.n_cells()])
            .collect();
        let mut state = LayeredState::new(vec![1.0e5; m.n_cells()], u, theta).with_humidity(q);
        let start = model.diagnostics(&state);
        let dt = model.stable_dt(&state, 0.25);
        let steps = (2.0 * 86_400.0 / dt).ceil() as u32;
        let dt = 2.0 * 86_400.0 / f64::from(steps);
        for _ in 0..steps {
            model.step(&mut state, dt);
        }
        let end = model.diagnostics(&state);
        let rel = |a: f64, b: f64| b / a - 1.0;
        eprintln!(
            "{n} layers, {steps} steps: mass {:+.1e}, θ {:+.1e}, water {:+.1e}, energy {:+.1e}",
            rel(start.mass, end.mass),
            rel(start.theta, end.theta),
            rel(start.water, end.water),
            rel(start.energy, end.energy)
        );
        assert!(rel(start.mass, end.mass).abs() < 1e-13);
        assert!(rel(start.theta, end.theta).abs() < 1e-13);
        assert!(rel(start.water, end.water).abs() < 1e-13);
        assert!(rel(start.energy, end.energy).abs() < 5e-8);
    }
}

// ---------------------------------------------------------------------------
// Vertical modes
// ---------------------------------------------------------------------------

#[test]
fn two_layers_reproduce_the_closed_form_linearization() {
    // The general continuity and thermodynamic rows of B against the
    // two-layer closed forms, at an interface where the arithmetic differs.
    let (m, g) = mesh(16);
    let model = model_with(m, g, vec![0.0; m.n_cells()], Levels::two(0.4));
    let (ps, theta) = (9.8e4, [325.0, 288.0]);
    let b = model.linearization(ps, &theta).b;
    let s = 0.4;
    let th_hat = 0.5 * (theta[0] + theta[1]);
    let closed = [
        [-ps * s, -ps * (1.0 - s)],
        [
            -(1.0 - s) * (theta[0] - th_hat),
            (1.0 - s) * (theta[0] - th_hat),
        ],
        [-s * (th_hat - theta[1]), s * (th_hat - theta[1])],
    ];
    for r in 0..3 {
        for c in 0..2 {
            let scale = closed[r][c].abs().max(1.0);
            assert!(
                (b[r][c] - closed[r][c]).abs() < 1e-14 * scale,
                "B[{r}][{c}] {} vs {}",
                b[r][c],
                closed[r][c]
            );
        }
    }
}

#[test]
fn gravity_wave_speeds_with_more_layers() {
    // An isothermal 260 K atmosphere: the external (Lamb) wave's speed
    // barely depends on the layers, √(γ R T) ≈ 323 m/s in the continuum
    // limit (Lamb waves in an isothermal atmosphere travel at the speed of
    // sound); each extra layer adds a slower internal mode.
    let (m, g) = mesh(16);
    let t0 = 260.0;
    let gamma = 1.0 / (1.0 - GAS.kappa());
    let sound = (gamma * GAS.gas_constant * t0).sqrt();
    let mut external = Vec::new();
    for levels in [
        Levels::two(0.5),
        Levels::equal_mass(3),
        Levels::equal_mass(5),
        Levels::equal_mass(10),
    ] {
        let n = levels.count();
        let model = model_with(m, g, vec![0.0; m.n_cells()], levels);
        let exner = model.layer_exner(1.0e5);
        let theta: Vec<f64> = exner.iter().map(|e| t0 / e).collect();
        let speeds = model.gravity_wave_speeds(1.0e5, &theta);
        eprintln!("{n} layers: {speeds:.1?} m/s");
        assert_eq!(speeds.len(), n);
        assert!(speeds.windows(2).all(|w| w[0] > w[1]), "{speeds:?}");
        assert!(speeds[n - 1] > 0.0);
        external.push(speeds[0]);
    }
    for c in &external {
        assert!(
            (c / sound - 1.0).abs() < 0.1,
            "external {c} vs sound {sound}"
        );
    }
}

// ---------------------------------------------------------------------------
// Rest over mountains
// ---------------------------------------------------------------------------

#[test]
fn isentropic_atmosphere_at_rest_stays_at_rest() {
    // Constant θ: balance means Φ_s + c_p θ Π_s is uniform, and the scheme
    // keeps that exactly, for any layers.
    let (m, g) = mesh(16);
    let phi_s = mountain(m);
    let theta0 = 300.0;
    let exner_ref = GAS.exner(1.0e5);
    let ps: Vec<f64> = phi_s
        .iter()
        .map(|phi| {
            let exner = exner_ref - phi / (GAS.heat_capacity * theta0);
            1.0e5 * fms_math::pow(exner, 1.0 / GAS.kappa())
        })
        .collect();
    for levels in layer_sets() {
        let n = levels.count();
        let state = LayeredState::new(
            ps.clone(),
            vec![vec![0.0; m.n_edges()]; n],
            vec![vec![theta0; m.n_cells()]; n],
        );
        let mut model = model_with(m, g, phi_s.clone(), levels);
        let mut d = state.clone();
        model.tendency(&state, &mut d);
        let worst = d.u.iter().flatten().fold(0.0f64, |a, &x| a.max(x.abs()));
        assert!(
            worst < 1e-12,
            "{n} layers: spurious acceleration {worst:e} m s-2"
        );
        assert!(d.ps.iter().all(|&x| x == 0.0));
    }
}

#[test]
fn isothermal_atmosphere_at_rest_stays_nearly_at_rest() {
    // A stratified atmosphere over terrain is the classic weak spot of sigma
    // coordinates: the pressure gradient splits into two large terms that
    // must cancel. Run 5 days over the 2 km cone and measure the spurious
    // wind. Measured: 0.105 m/s (two and three layers), 0.103 (five).
    let (m, g) = mesh(16);
    let phi_s = mountain(m);
    let t0 = 260.0;
    let ps: Vec<f64> = phi_s
        .iter()
        .map(|phi| 1.0e5 * fms_math::exp(-phi / (GAS.gas_constant * t0)))
        .collect();
    for levels in layer_sets() {
        let n = levels.count();
        let mut model = model_with(m, g, phi_s.clone(), levels);
        let mut theta = vec![vec![0.0; m.n_cells()]; n];
        for (i, &p) in ps.iter().enumerate() {
            let exner = model.layer_exner(p);
            for k in 0..n {
                theta[k][i] = t0 / exner[k];
            }
        }
        let mut state = LayeredState::new(ps.clone(), vec![vec![0.0; m.n_edges()]; n], theta);
        let dt = model.stable_dt(&state, 0.5);
        let steps = (5.0 * 86_400.0 / dt).ceil() as u32;
        let dt = 5.0 * 86_400.0 / f64::from(steps);
        let mut worst = 0.0f64;
        for _ in 0..steps {
            model.step(&mut state, dt);
            worst = worst.max(
                state
                    .u
                    .iter()
                    .flatten()
                    .fold(0.0f64, |a, &x| a.max(x.abs())),
            );
        }
        eprintln!(
            "{n} layers, isothermal rest over a 2 km cone: max spurious wind {worst:.3} m/s over 5 days"
        );
        assert!(worst < 1.0, "{n} layers: spurious wind {worst} m/s");
    }
}

// ---------------------------------------------------------------------------
// Baroclinic flows
// ---------------------------------------------------------------------------

/// Earth-like sheared flow: upper layer `30·shear_scale` m/s faster than the
/// lower at the equator, 40 K static stability.
fn sheared(shear_scale: f64) -> ShearedFlow {
    ShearedFlow {
        speed: [10.0 + 30.0 * shear_scale, 10.0],
        stability: 40.0,
        theta_equator: 315.0,
        ps_equator: 1.0e5,
        omega: OMEGA,
    }
}

fn baroclinic_state(model: &LayeredModel, shear_scale: f64) -> LayeredState {
    sheared(shear_scale).state(model)
}

/// θ̂ = (θ₀ + θ₁)/2 at every cell.
fn theta_hat(s: &LayeredState) -> Vec<f64> {
    (0..s.ps.len())
        .map(|i| 0.5 * (s.theta_ps[0][i] + s.theta_ps[1][i]) / s.ps[i])
        .collect()
}

/// Amplitude (K) of zonal wavenumber `wave` in `a − b` over 30–60° N.
fn wave_amplitude(m: &SphereMesh, g: &MeshGeometry, a: &[f64], b: &[f64], wave: u32) -> f64 {
    let (mut re, mut im, mut area) = (0.0, 0.0, 0.0);
    for (i, &p) in m.cell_centers().iter().enumerate() {
        let (lat, lon) = to_lat_lon(p);
        if (30f64.to_radians()..60f64.to_radians()).contains(&lat) {
            let (s, c) = fms_math::sin_cos(f64::from(wave) * lon);
            let d = g.cell_areas[i] * (a[i] - b[i]);
            re += d * c;
            im += d * s;
            area += g.cell_areas[i];
        }
    }
    2.0 * (re * re + im * im).sqrt() / area
}

/// Growth rate (per day) of a 0.01 K wavenumber-`wave` θ seed at 45° N on
/// the sheared flow, measured against an unperturbed control run (which
/// removes the zonal mean's aliasing onto the irregular grid and the slow
/// drift of the background). Fitted from day 2 to day 8.
fn measured_growth(flow: ShearedFlow, wave: u32) -> f64 {
    let (m, g) = mesh(16);
    let mut control_model = model(m, g, vec![0.0; m.n_cells()]);
    let mut perturbed_model = model(m, g, vec![0.0; m.n_cells()]);
    let mut control = flow.state(&control_model);
    let mut perturbed = control.clone();
    for k in 0..2 {
        for (i, &p) in m.cell_centers().iter().enumerate() {
            let (lat, lon) = to_lat_lon(p);
            let bump = fms_math::exp(-((lat - 45f64.to_radians()) / 15f64.to_radians()).powi(2));
            perturbed.theta_ps[k][i] +=
                perturbed.ps[i] * 0.01 * bump * fms_math::cos(f64::from(wave) * lon);
        }
    }
    let dt = control_model.stable_dt(&control, 0.5);
    let per_day = (86_400.0 / dt).ceil() as u32;
    let dt = 86_400.0 / f64::from(per_day);
    let mut samples = Vec::new();
    for day in 0..=8 {
        if day > 0 {
            for _ in 0..per_day {
                control_model.step(&mut control, dt);
                perturbed_model.step(&mut perturbed, dt);
            }
        }
        let a = wave_amplitude(m, g, &theta_hat(&perturbed), &theta_hat(&control), wave);
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
    num / den
}

#[test]
fn baroclinic_instability_matches_linear_theory() {
    // Two-level quasi-geostrophic theory (Phillips 1954; Holton §8.2) on the
    // local β-plane at 45° N, with the model's own internal gravity-wave
    // speed setting the deformation radius. Below the critical shear nothing
    // grows; above it waves grow at close to the predicted rate. Local
    // theory ignores the waves' north–south structure, which on the sphere
    // (an unstable zone ~30° wide) slows growth by ~25%, so it is an upper
    // bound. Measured here (k = 16, days 2–8): 0.033 /day at 10 m/s;
    // 0.83 /day against 1.10 /day theory at 45 m/s (ratio 0.75). The
    // `baroclinic` example at k = 32 gives 0.79 (ratio 0.71), with zonal
    // wavenumber 6 the fastest-growing in both model and theory.
    let (m, g) = mesh(16);
    let reference = model(m, g, vec![0.0; m.n_cells()]);
    let lat = 45f64.to_radians();

    let stable = sheared(10.0 / 30.0);
    let (theory, critical) = stable.phillips_growth_rate(&reference, 6, lat);
    assert_eq!(
        theory, 0.0,
        "10 m/s should be below the critical {critical:.1} m/s"
    );
    let growth = measured_growth(stable, 6);
    eprintln!("shear 10 m/s (critical {critical:.1}): growth {growth:.3} /day");
    assert!(growth < 0.05, "subcritical wave grew at {growth} /day");

    let unstable = sheared(45.0 / 30.0);
    let (theory, _) = unstable.phillips_growth_rate(&reference, 6, lat);
    let growth = measured_growth(unstable, 6);
    let ratio = growth / (theory * 86_400.0);
    eprintln!(
        "shear 45 m/s: growth {growth:.3} /day, theory {:.3} /day, ratio {ratio:.2}",
        theory * 86_400.0
    );
    assert!((0.6..1.0).contains(&ratio), "growth ratio {ratio}");
}

#[test]
fn the_linearization_matches_the_tendencies() {
    // B against the model's own tendencies: rest, uniform p̄ and stratified
    // θ̄ₖ, a small divergent wind in one layer at a time. (With more than one
    // interface the gravity-wave matrix A = J B isn't self-adjoint in the
    // mass-weighted inner product, unlike with two layers: the arithmetic
    // interface θ makes it slightly non-normal. Checking B directly shows
    // that is the scheme's, not a slip in B.)
    let (m, g) = mesh(16);
    let ops = Operators::new(m, g, RADIUS);
    for levels in [Levels::new(&[0.4, 0.85]), Levels::equal_mass(5)] {
        let mut model = model_with(m, g, vec![0.0; m.n_cells()], levels);
        let n = model.layers();
        let p_bar = 1.0e5;
        let exner = model.layer_exner(p_bar);
        let theta: Vec<f64> = exner.iter().map(|e| 260.0 / e).collect();
        let b = model.linearization(p_bar, &theta).b;
        // A smooth potential χ; u = ε ∇χ at edges.
        let chi: Vec<f64> = m
            .cell_centers()
            .iter()
            .map(|p| p.x + 0.5 * p.y * p.z)
            .collect();
        let mut grad = vec![0.0; m.n_edges()];
        ops.gradient(&chi, &mut grad);
        let eps = 1e-6;
        let mut worst = 0.0f64;
        for l in 0..n {
            let mut u = vec![vec![0.0; m.n_edges()]; n];
            u[l] = grad.iter().map(|x| eps * x).collect();
            let state = LayeredState::new(
                vec![p_bar; m.n_cells()],
                u.clone(),
                theta.iter().map(|&t| vec![t; m.n_cells()]).collect(),
            );
            let mut d = state.clone();
            model.tendency(&state, &mut d);
            let mut delta = vec![0.0; m.n_cells()];
            ops.divergence(&u[l], &mut delta);
            for i in (0..m.n_cells()).step_by(97) {
                if delta[i].abs() < 1e-3 * delta.iter().fold(0.0f64, |a, x| a.max(x.abs())) {
                    continue;
                }
                let dps = d.ps[i];
                let rel = (dps / delta[i] - b[0][l]) / b[0][l].abs().max(1.0);
                worst = worst.max(rel.abs());
                for k in 0..n {
                    let dy = (d.theta_ps[k][i] - theta[k] * dps) / p_bar;
                    let got = dy / delta[i];
                    let scale = b
                        .iter()
                        .skip(1)
                        .flatten()
                        .fold(0.0f64, |a, x| a.max(x.abs()));
                    worst = worst.max(((got - b[k + 1][l]) / scale).abs());
                }
            }
        }
        assert!(worst < 1e-10, "{n} layers: B mismatch {worst:e}");
    }
}
