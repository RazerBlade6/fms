//! Semi-implicit leapfrog time stepping for the two-layer model, with the
//! RAW filter: the two-layer version of [`crate::semi_implicit`].
//!
//! # Splitting
//!
//! About a resting, flat reference state (`p̄`, `θ̄₀`, `θ̄₁`; by default the
//! initial global means), the gravity-wave terms are linear in
//! `y = (p_s, θ₀ᴸ, θ₁ᴸ)` with `θₖᴸ = (Θₖ − θ̄ₖ p_s)/p̄`, and the layer divergences
//! `δₖ = ∇·uₖ`:
//!
//! ```text
//! ∂y/∂t = B δ,     ∂uₖ/∂t = −∇Mₖ,   M = J y
//! ```
//!
//! ([`Linearization`]). These `L` terms are averaged over the old and new
//! time levels; the rest, `N = T − L` (advection, potential-vorticity flux,
//! topography, the nonlinear part of the pressure gradient), is explicit at
//! the middle level:
//!
//! ```text
//! xⁿ⁺¹ = xⁿ⁻¹ + 2Δt N(xⁿ) + Δt L(xⁿ⁺¹ + xⁿ⁻¹)
//! ```
//!
//! # Solving
//!
//! With `Y = yⁿ⁻¹ + 2Δt N_y + Δt B δⁿ⁻¹` and
//! `Uₖ = uₖⁿ⁻¹ + 2Δt N_u − Δt ∇Mₖ(yⁿ⁻¹)`, eliminating the new winds gives
//!
//! ```text
//! (I + Δt² A ∇²) Mⁿ⁺¹ = J Y + Δt A ∇·U,     A = J B
//! ```
//!
//! `A` has eigenvalues `−c²` for the external and internal gravity waves, so
//! in its eigenvectors the system splits into two scalar Helmholtz problems
//! `(I − Δt² c² ∇²) M̃ = r̃`, each solved by conjugate gradient. Then
//! `uⁿ⁺¹ = U − Δt ∇Mⁿ⁺¹` and `yⁿ⁺¹ = Y + Δt B ∇·uⁿ⁺¹`.
//!
//! Because the new surface pressure and `Θ` come from the divergence of the
//! new winds, which has no net source, total mass and total `Θ` are conserved
//! to round-off whatever the solver tolerance.
//!
//! # Diffusion
//!
//! The model's hyperdiffusion ([`crate::diffusion`]) and forcing
//! ([`crate::forcing`]) are evaluated at the old level and added to `N`:
//! `xⁿ⁺¹ = xⁿ⁻¹ + 2Δt D(xⁿ⁻¹) + …`. Centred in time, damping would be
//! unstable under leapfrog; lagged, it is stable while `Δt` times the fastest
//! damping rate is below 1. The energy each adds or removes is accumulated in
//! [`TwoLayerSemiImplicit::energy_budget`].
//!
//! # Humidity
//!
//! A moist state's humidity is carried by flux-corrected transport
//! ([`crate::transport`]) with the mass fluxes the step used, so it stays
//! positive and consistent with the surface pressure; the physics sources
//! (lagged, like the rest of the forcing) go in with it. The RAW filter acts
//! on `Q` too, and can leave a trace of negative humidity, which is filled
//! from the column's other layer where it can be and otherwise set to zero.
//! [`TwoLayerSemiImplicit::water_budget`] accounts for every kilogram:
//! because leapfrog advances two interleaved levels, the water it closes on
//! is the mean of the two current levels ([`TwoLayerSemiImplicit::water`]).

use fms_geo::Operators;
use fms_geo::operators::MIN_CHUNK;
use rayon::prelude::*;

use crate::forcing::EnergyTerms;
use crate::helmholtz::{Helmholtz, SolveReport, dot};
use crate::safety::{self, BlowUp, Health, Limits};
use crate::semi_implicit::SolverStats;
use crate::transport::{Transport, TransportStep};
use crate::two_layer::{Linearization, TwoLayer, TwoLayerState};

/// Settings for [`TwoLayerSemiImplicit`].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TwoLayerSemiImplicitOptions {
    /// Reference `(p̄, [θ̄₀, θ̄₁])` for the linearization; `None` uses the
    /// initial state's global means (area-weighted `p_s`, mass-weighted θ).
    pub reference: Option<(f64, [f64; 2])>,
    /// RAW filter strength `ν`.
    pub filter_strength: f64,
    /// RAW filter weight `α`.
    pub filter_alpha: f64,
    /// Helmholtz solver tolerance (relative residual).
    pub tolerance: f64,
}

impl Default for TwoLayerSemiImplicitOptions {
    fn default() -> Self {
        Self {
            reference: None,
            filter_strength: 0.05,
            filter_alpha: 0.53,
            tolerance: 1e-8,
        }
    }
}

/// Energy added since the start (J) by each non-adiabatic process. With it,
/// `energy(now) − energy(start) − diffusion − forcing` is the time-stepping
/// error.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct EnergyBudget {
    /// Hyperdiffusion (negative).
    pub diffusion: EnergyTerms,
    /// External forcing: for Held–Suarez, friction (kinetic) and Newtonian
    /// heating and cooling (thermal).
    pub forcing: EnergyTerms,
}

/// Water added since the start (kg), by cause. With it,
/// `water(now) − water(start) − sources − filter − fixer` is round-off
/// ([`TwoLayerSemiImplicit::water`]).
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct WaterBudget {
    /// The forcing's sources: evaporation minus precipitation.
    pub sources: f64,
    /// The RAW time filter (small; it is not exactly conservative when the
    /// levels' totals differ).
    pub filter: f64,
    /// Negative humidity filled with zero (a trace).
    pub fixer: f64,
    /// Water moved between layers to fill negative humidity (not added; for
    /// information).
    pub borrowed: f64,
    /// Most upwind substeps any transport step needed.
    pub max_substeps: u32,
}

/// The vertical modes: `A = V diag(−c²) V⁻¹`.
#[derive(Clone, Copy, Debug)]
struct Modes {
    /// Columns are the eigenvectors (external first).
    v: [[f64; 2]; 2],
    v_inv: [[f64; 2]; 2],
    /// `c²` per mode.
    speed2: [f64; 2],
}

impl Modes {
    fn of(a: [[f64; 2]; 2]) -> Self {
        let (tr, det) = (a[0][0] + a[1][1], a[0][0] * a[1][1] - a[0][1] * a[1][0]);
        let disc = (0.25 * tr * tr - det).max(0.0).sqrt();
        // Most negative eigenvalue (fastest wave, external) first.
        let lambda = [0.5 * tr - disc, 0.5 * tr + disc];
        let vector = |l: f64, fallback: [f64; 2]| -> [f64; 2] {
            let (x, y) = if a[0][1].abs() >= a[1][0].abs() && a[0][1] != 0.0 {
                (a[0][1], l - a[0][0])
            } else if a[1][0] != 0.0 {
                (l - a[1][1], a[1][0])
            } else {
                (fallback[0], fallback[1])
            };
            let n = (x * x + y * y).sqrt();
            [x / n, y / n]
        };
        let (e0, e1) = (vector(lambda[0], [1.0, 0.0]), vector(lambda[1], [0.0, 1.0]));
        let v = [[e0[0], e1[0]], [e0[1], e1[1]]];
        let d = v[0][0] * v[1][1] - v[0][1] * v[1][0];
        let v_inv = [[v[1][1] / d, -v[0][1] / d], [-v[1][0] / d, v[0][0] / d]];
        Self {
            v,
            v_inv,
            speed2: lambda.map(|l| (-l).max(0.0)),
        }
    }
}

/// Two-layer model stepped semi-implicitly.
pub struct TwoLayerSemiImplicit<'a> {
    model: TwoLayer<'a>,
    ops: Operators<'a>,
    helmholtz: [Helmholtz<'a>; 2],
    lin: Linearization,
    /// `A = J B`.
    a: [[f64; 2]; 2],
    modes: Modes,
    reference: (f64, [f64; 2]),
    options: TwoLayerSemiImplicitOptions,
    dt: f64,
    previous: TwoLayerState,
    current: TwoLayerState,
    steps: u64,
    stats: SolverStats,
    budget: EnergyBudget,
    water: WaterBudget,
    transport: Option<Transport>,
    // Scratch.
    tendency: TwoLayerState,
    /// Lagged tendencies: diffusion plus forcing.
    physics: TwoLayerState,
    forcing: TwoLayerState,
    next: TwoLayerState,
    big_u: [Vec<f64>; 2],
    m: [Vec<f64>; 2],
    modal: [Vec<f64>; 2],
    rhs: [Vec<f64>; 2],
    div: [Vec<f64>; 2],
    div_prev: [Vec<f64>; 2],
    grad: Vec<f64>,
    cells: Vec<f64>,
    mass_flux: [Vec<f64>; 2],
    zeros: Vec<f64>,
}

impl<'a> TwoLayerSemiImplicit<'a> {
    /// Start from `state` with step `dt` (s). Takes the first step at once,
    /// with RK4 substeps at Courant number `start_courant` (e.g. 0.5).
    ///
    /// # Panics
    /// If the model has hyperdiffusion too strong for `dt`:
    /// `dt · max_rate ≥ 1` ([`crate::Hyperdiffusion::max_rate`]).
    pub fn new(
        mut model: TwoLayer<'a>,
        state: TwoLayerState,
        dt: f64,
        options: TwoLayerSemiImplicitOptions,
        start_courant: f64,
    ) -> Self {
        if let Some(d) = model.hyperdiffusion() {
            assert!(
                dt * d.max_rate() < 1.0,
                "lagged hyperdiffusion is unstable at this step: Δt × fastest damping rate = {:.2} (must be < 1); lengthen the e-folding time",
                dt * d.max_rate()
            );
        }
        let ops = *model.operators();
        let reference = options
            .reference
            .unwrap_or_else(|| global_means(&ops, &state, model.dsigma()));
        let lin = model.linearization(reference.0, reference.1);
        let a = model.gravity_wave_matrix(reference.0, reference.1);
        let modes = Modes::of(a);
        let helmholtz = modes.speed2.map(|c2| {
            let mut h = Helmholtz::new(ops, dt * dt * c2);
            h.tolerance = options.tolerance;
            h
        });

        let substeps = (dt / model.stable_dt(&state, start_courant))
            .ceil()
            .max(1.0);
        let mut current = state.clone();
        let (mut physics, mut forcing) = (state.clone(), state.clone());
        let mut budget = EnergyBudget::default();
        model.diffusion_tendency(&state, &mut physics);
        budget
            .diffusion
            .add_scaled(dt, model.energy_tendency(&state, &physics));
        model.forcing_tendency(&state, &mut forcing);
        budget
            .forcing
            .add_scaled(dt, model.energy_tendency(&state, &forcing));
        for _ in 0..substeps as u64 {
            model.step(&mut current, dt / substeps);
        }
        let (nc, ne) = (state.ps.len(), state.u[0].len());
        let cells = || vec![0.0; nc];
        let edges = || vec![0.0; ne];
        let moist = state.is_moist();
        if moist {
            // The RK4 start's centred transport isn't monotone: keep its
            // humidity within the initial range nearby (the cell, its
            // neighbours and the other layer). What this changes is part of
            // the starting water, not the budget.
            let mesh = ops.mesh();
            for k in 0..2 {
                let q = |l: usize, i: usize| state.q_ps[l][i] / state.ps[i];
                for i in 0..nc {
                    let (mut lo, mut hi) = (q(k, i).min(q(1 - k, i)), q(k, i).max(q(1 - k, i)));
                    for &n in mesh.cell_neighbors(i as u32) {
                        lo = lo.min(q(k, n as usize));
                        hi = hi.max(q(k, n as usize));
                    }
                    let p = current.ps[i];
                    current.q_ps[k][i] = current.q_ps[k][i].clamp(p * lo.max(0.0), p * hi);
                }
            }
        }
        let mut solver = Self {
            model,
            ops,
            helmholtz,
            lin,
            a,
            modes,
            reference,
            options,
            dt,
            previous: state.clone(),
            current,
            steps: 1,
            stats: SolverStats::default(),
            budget,
            water: WaterBudget::default(),
            transport: moist.then(|| Transport::new(nc, ne)),
            tendency: state.dry(),
            physics,
            forcing,
            next: state,
            big_u: [edges(), edges()],
            m: [cells(), cells()],
            modal: [cells(), cells()],
            rhs: [cells(), cells()],
            div: [cells(), cells()],
            div_prev: [cells(), cells()],
            grad: edges(),
            cells: cells(),
            mass_flux: [edges(), edges()],
            zeros: cells(),
        };
        // First guess for the modal potentials: those of the current state.
        let y = solver.y_of(&solver.current.clone());
        let m = solver.potentials(&y);
        solver.modal = solver.to_modes(&m);
        solver
    }

    /// The current state.
    pub fn state(&self) -> &TwoLayerState {
        &self.current
    }

    /// The current state, mutably: for experiments such as injecting a
    /// fault. Only the current time level changes, not the previous one.
    pub fn state_mut(&mut self) -> &mut TwoLayerState {
        &mut self.current
    }

    /// Energy added by hyperdiffusion since the start (J; negative).
    pub fn dissipated_energy(&self) -> f64 {
        self.budget.diffusion.total()
    }

    /// Energy added by each non-adiabatic process since the start.
    pub fn energy_budget(&self) -> EnergyBudget {
        self.budget
    }

    /// Water added since the start, by cause.
    pub fn water_budget(&self) -> WaterBudget {
        self.water
    }

    /// Total water vapour (kg): the mean of the two current leapfrog levels,
    /// the quantity [`water_budget`](Self::water_budget) closes on exactly.
    pub fn water(&self) -> f64 {
        0.5 * (self.model.water(&self.previous) + self.model.water(&self.current))
    }

    /// Fastest wind (m s⁻¹) the step can advect stably: Courant number 1 on
    /// the shortest dual edge.
    pub fn courant_limit(&self) -> f64 {
        let min_d = self
            .ops
            .geometry()
            .dual_edge_lengths
            .iter()
            .fold(f64::MAX, |a, &b| a.min(b))
            * self.ops.radius();
        min_d / self.dt
    }

    /// Safety limits for this run: the Courant limit and bounds about the
    /// reference state ([`Limits::around`]).
    pub fn limits(&self) -> Limits {
        let (p_bar, th_bar) = self.reference;
        Limits::around(self.courant_limit(), p_bar, th_bar)
    }

    /// Check the current state against [`limits`](Self::limits).
    pub fn check(&self) -> Result<Health, BlowUp> {
        safety::check(&self.current, &self.limits())
    }

    /// The underlying model.
    pub fn model(&self) -> &TwoLayer<'a> {
        &self.model
    }

    /// Steps taken, including the RK4 start.
    pub fn steps(&self) -> u64 {
        self.steps
    }

    /// The time step (s).
    pub fn dt(&self) -> f64 {
        self.dt
    }

    /// Reference `(p̄, [θ̄₀, θ̄₁])` of the linearization.
    pub fn reference(&self) -> (f64, [f64; 2]) {
        self.reference
    }

    /// Gravity-wave speeds `[external, internal]` of the reference state
    /// (m s⁻¹).
    pub fn wave_speeds(&self) -> [f64; 2] {
        self.modes.speed2.map(f64::sqrt)
    }

    /// Helmholtz solver totals (both modes).
    pub fn solver_stats(&self) -> SolverStats {
        self.stats
    }

    /// The longest stable step (s) for winds up to `max_wind` (m s⁻¹) at
    /// Courant number `courant`: advection and the inertial limit, as for the
    /// shallow-water scheme ([`crate::SemiImplicit::stable_dt`]).
    pub fn stable_dt(model: &TwoLayer, max_wind: f64, courant: f64) -> f64 {
        let ops = model.operators();
        let min_d = ops
            .geometry()
            .dual_edge_lengths
            .iter()
            .fold(f64::MAX, |a, &b| a.min(b))
            * ops.radius();
        let max_f = model.coriolis().iter().fold(0.0f64, |a, &f| a.max(f.abs()));
        let advective = if max_wind > 0.0 {
            min_d / max_wind
        } else {
            f64::MAX
        };
        let inertial = if max_f > 0.0 { 1.0 / max_f } else { f64::MAX };
        courant * advective.min(inertial)
    }

    /// `y = (p_s, θ₀ᴸ, θ₁ᴸ)` of a state.
    fn y_of(&self, s: &TwoLayerState) -> [Vec<f64>; 3] {
        let (p_bar, th_bar) = self.reference;
        let theta_l = |k: usize| -> Vec<f64> {
            s.theta_ps[k]
                .par_iter()
                .with_min_len(MIN_CHUNK)
                .zip(&s.ps)
                .map(|(t, p)| (t - th_bar[k] * p) / p_bar)
                .collect()
        };
        [s.ps.clone(), theta_l(0), theta_l(1)]
    }

    /// `M = J y` per layer.
    fn potentials(&self, y: &[Vec<f64>; 3]) -> [Vec<f64>; 2] {
        let j = self.lin.j;
        std::array::from_fn(|k| {
            (0..y[0].len())
                .into_par_iter()
                .with_min_len(MIN_CHUNK)
                .map(|i| j[k][0] * y[0][i] + j[k][1] * y[1][i] + j[k][2] * y[2][i])
                .collect()
        })
    }

    fn to_modes(&self, m: &[Vec<f64>; 2]) -> [Vec<f64>; 2] {
        let w = self.modes.v_inv;
        std::array::from_fn(|r| {
            m[0].par_iter()
                .with_min_len(MIN_CHUNK)
                .zip(&m[1])
                .map(|(a, b)| w[r][0] * a + w[r][1] * b)
                .collect()
        })
    }

    /// Advance one step. Returns the solve reports of the external and
    /// internal modes.
    // Loops over the two layers index several per-layer arrays together.
    #[allow(clippy::needless_range_loop)]
    pub fn step(&mut self) -> [SolveReport; 2] {
        let ops = self.ops;
        let dt = self.dt;
        let (p_bar, th_bar) = self.reference;
        let (j, b, a) = (self.lin.j, self.lin.b, self.a);

        // Full tendency at the middle level.
        let empty = TwoLayerState::new(
            Vec::new(),
            [Vec::new(), Vec::new()],
            [Vec::new(), Vec::new()],
        );
        let mut current = std::mem::replace(&mut self.current, empty);
        self.model.tendency(&current, &mut self.tendency);
        // Hyperdiffusion and forcing, lagged to the old level.
        let diffuse = self.model.hyperdiffusion().is_some();
        let forced = self.model.has_forcing();
        if diffuse {
            self.model
                .diffusion_tendency(&self.previous, &mut self.physics);
            let rate = self.model.energy_tendency(&self.previous, &self.physics);
            self.budget.diffusion.add_scaled(dt, rate);
        }
        if forced {
            self.model
                .forcing_tendency(&self.previous, &mut self.forcing);
            let rate = self.model.energy_tendency(&self.previous, &self.forcing);
            self.budget.forcing.add_scaled(dt, rate);
            self.water.sources += dt * self.model.water_tendency(&self.forcing);
            if diffuse {
                for (p, f) in self
                    .physics
                    .parts_mut()
                    .into_iter()
                    .zip(self.forcing.parts())
                {
                    p.par_iter_mut()
                        .with_min_len(MIN_CHUNK)
                        .zip(f.par_iter())
                        .for_each(|(p, f)| *p += f);
                }
            } else {
                std::mem::swap(&mut self.physics, &mut self.forcing);
            }
        }
        let diffuse = diffuse || forced;
        let t = &self.tendency;
        let d = &self.physics;

        // Divergences of uⁿ and uⁿ⁻¹.
        for k in 0..2 {
            ops.divergence(&current.u[k], &mut self.div[k]);
            ops.divergence(&self.previous.u[k], &mut self.div_prev[k]);
        }
        let (div, div_prev) = (&self.div, &self.div_prev);

        // Y = yⁿ⁻¹ + 2Δt (T_y − B δⁿ) + Δt B δⁿ⁻¹, in y-space.
        let y_prev = self.y_of(&self.previous);
        let nc = current.ps.len();
        let mut big_y = [vec![0.0; nc], vec![0.0; nc], vec![0.0; nc]];
        for (r, out) in big_y.iter_mut().enumerate() {
            out.par_iter_mut()
                .with_min_len(MIN_CHUNK)
                .enumerate()
                .for_each(|(i, o)| {
                    let t_y = if r == 0 {
                        t.ps[i]
                    } else {
                        let diffusion = if diffuse { d.theta_ps[r - 1][i] } else { 0.0 };
                        (t.theta_ps[r - 1][i] + diffusion - th_bar[r - 1] * t.ps[i]) / p_bar
                    };
                    let b_now = b[r][0] * div[0][i] + b[r][1] * div[1][i];
                    let b_old = b[r][0] * div_prev[0][i] + b[r][1] * div_prev[1][i];
                    *o = y_prev[r][i] + 2.0 * dt * (t_y - b_now) + dt * b_old;
                });
        }

        // Uₖ = uₖⁿ⁻¹ + 2Δt (T_u + ∇Mₖ(yⁿ)) − Δt ∇Mₖ(yⁿ⁻¹).
        let y_now = self.y_of(&current);
        let (m_now, m_prev) = (self.potentials(&y_now), self.potentials(&y_prev));
        for k in 0..2 {
            // ∇(2 Mₖ(yⁿ) − Mₖ(yⁿ⁻¹)) combines both linear terms in one gradient.
            self.cells
                .par_iter_mut()
                .with_min_len(MIN_CHUNK)
                .enumerate()
                .for_each(|(i, c)| *c = 2.0 * m_now[k][i] - m_prev[k][i]);
            ops.gradient(&self.cells, &mut self.grad);
            let (prev_u, t_u, d_u, grad) = (&self.previous.u[k], &t.u[k], &d.u[k], &self.grad);
            self.big_u[k]
                .par_iter_mut()
                .with_min_len(MIN_CHUNK)
                .enumerate()
                .for_each(|(e, u)| {
                    let diffusion = if diffuse { d_u[e] } else { 0.0 };
                    *u = prev_u[e] + 2.0 * dt * (t_u[e] + diffusion) + dt * grad[e];
                });
        }

        // rhs = J Y + Δt A ∇·U.
        let mut div_u = [vec![0.0; nc], vec![0.0; nc]];
        for (u, d) in self.big_u.iter().zip(&mut div_u) {
            ops.divergence(u, d);
        }
        for (k, out) in self.rhs.iter_mut().enumerate() {
            out.par_iter_mut()
                .with_min_len(MIN_CHUNK)
                .enumerate()
                .for_each(|(i, r)| {
                    *r = j[k][0] * big_y[0][i]
                        + j[k][1] * big_y[1][i]
                        + j[k][2] * big_y[2][i]
                        + dt * (a[k][0] * div_u[0][i] + a[k][1] * div_u[1][i]);
                });
        }

        // Two scalar Helmholtz problems in the vertical modes.
        let rhs_modal = self.to_modes(&self.rhs);
        let mut reports = [SolveReport {
            iterations: 0,
            relative_residual: 0.0,
        }; 2];
        for (mode, report) in reports.iter_mut().enumerate() {
            *report = self.helmholtz[mode].solve(&rhs_modal[mode], &mut self.modal[mode]);
        }
        let v = self.modes.v;
        for (k, out) in self.m.iter_mut().enumerate() {
            let modal = &self.modal;
            out.par_iter_mut()
                .with_min_len(MIN_CHUNK)
                .enumerate()
                .for_each(|(i, m)| *m = v[k][0] * modal[0][i] + v[k][1] * modal[1][i]);
        }

        // uⁿ⁺¹ = U − Δt ∇Mⁿ⁺¹; yⁿ⁺¹ = Y + Δt B ∇·uⁿ⁺¹.
        let next = &mut self.next;
        for k in 0..2 {
            ops.gradient(&self.m[k], &mut self.grad);
            let (big_u, grad) = (&self.big_u[k], &self.grad);
            next.u[k]
                .par_iter_mut()
                .with_min_len(MIN_CHUNK)
                .enumerate()
                .for_each(|(e, u)| *u = big_u[e] - dt * grad[e]);
            ops.divergence(&next.u[k], &mut self.div[k]);
        }
        let div = &self.div;
        let y_new: [Vec<f64>; 3] = std::array::from_fn(|r| {
            (0..nc)
                .into_par_iter()
                .with_min_len(MIN_CHUNK)
                .map(|i| big_y[r][i] + dt * (b[r][0] * div[0][i] + b[r][1] * div[1][i]))
                .collect()
        });
        next.ps.copy_from_slice(&y_new[0]);
        for k in 0..2 {
            let (ps, yk) = (&y_new[0], &y_new[k + 1]);
            next.theta_ps[k]
                .par_iter_mut()
                .with_min_len(MIN_CHUNK)
                .enumerate()
                .for_each(|(i, t)| *t = p_bar * yk[i] + th_bar[k] * ps[i]);
        }

        // Humidity: flux-corrected transport with the step's mass fluxes,
        // F̄ₖ = Fₖⁿ + p̄ δσₖ (uₖⁿ⁻¹ + uₖⁿ⁺¹ − 2uₖⁿ)/2.
        if let Some(transport) = &mut self.transport {
            let ds = self.model.dsigma();
            ops.cell_to_edge(&current.ps, &mut self.grad);
            for k in 0..2 {
                let (ps_e, now, old, new) =
                    (&self.grad, &current.u[k], &self.previous.u[k], &next.u[k]);
                self.mass_flux[k]
                    .par_iter_mut()
                    .with_min_len(MIN_CHUNK)
                    .enumerate()
                    .for_each(|(e, f)| {
                        *f = ds[k]
                            * (ps_e[e] * now[e] + 0.5 * p_bar * (old[e] + new[e] - 2.0 * now[e]));
                    });
            }
            let mut d0 = vec![0.0; nc];
            ops.divergence(&self.mass_flux[0], &mut d0);
            ops.divergence(&self.mass_flux[1], &mut self.cells);
            let interface: Vec<f64> = d0
                .par_iter()
                .with_min_len(MIN_CHUNK)
                .zip(self.cells.par_iter())
                .map(|(d0, d1)| ds[0] * d1 - ds[1] * d0)
                .collect();
            let q_mid: [Vec<f64>; 2] = std::array::from_fn(|k| current.humidity(k));
            let ratio = self.model.interface_ratios(&current);
            let source = if diffuse {
                [d.q_ps[0].as_slice(), d.q_ps[1].as_slice()]
            } else {
                [self.zeros.as_slice(), self.zeros.as_slice()]
            };
            let input = TransportStep {
                dsigma: ds,
                interval: 2.0 * dt,
                ps_old: &self.previous.ps,
                ps_new: &next.ps,
                flux: [&self.mass_flux[0], &self.mass_flux[1]],
                interface: &interface,
                q_old: [&self.previous.q_ps[0], &self.previous.q_ps[1]],
                q_mid: [&q_mid[0], &q_mid[1]],
                source,
                interface_ratio: [&ratio[0], &ratio[1]],
            };
            let [q0, q1] = &mut next.q_ps;
            let report = transport.step(&ops, &input, [q0, q1]);
            self.water.max_substeps = self.water.max_substeps.max(report.substeps);
        }

        // RAW filter on xⁿ and xⁿ⁺¹, then rotate the three levels (the old
        // xⁿ⁻¹ buffer becomes scratch for the next step).
        let (nu, alpha) = (self.options.filter_strength, self.options.filter_alpha);
        // Humidity's range over the three levels, per cell and layer, and the
        // water of the two levels the filter changes.
        let moist = self.transport.is_some();
        let (bounds, water_before) = if moist {
            let bounds: [Vec<(f64, f64)>; 2] = std::array::from_fn(|k| {
                (0..nc)
                    .into_par_iter()
                    .with_min_len(MIN_CHUNK)
                    .map(|i| {
                        let q = [
                            self.previous.q_ps[k][i] / self.previous.ps[i],
                            current.q_ps[k][i] / current.ps[i],
                            next.q_ps[k][i] / next.ps[i],
                        ];
                        (q[0].min(q[1]).min(q[2]), q[0].max(q[1]).max(q[2]))
                    })
                    .collect()
            });
            let w = self.model.water(&current) + self.model.water(next);
            (Some(bounds), w)
        } else {
            (None, 0.0)
        };
        raw_filter(&self.previous, &mut current, next, nu, alpha);
        if let Some(bounds) = bounds {
            // Keep the filtered humidity within the range of the three levels
            // (the filter is linear, so this never acts on a uniform q), and
            // account for the water the filter changed.
            for level in [&mut current, &mut *next] {
                for k in 0..2 {
                    let (ps, b) = (&level.ps, &bounds[k]);
                    level.q_ps[k]
                        .par_iter_mut()
                        .with_min_len(MIN_CHUNK)
                        .enumerate()
                        .for_each(|(i, q)| *q = q.clamp(ps[i] * b[i].0, ps[i] * b[i].1));
                }
            }
            let after = self.model.water(&current) + self.model.water(next);
            self.water.filter += 0.5 * (after - water_before);
            let (ds, g) = (self.model.dsigma(), self.model.gravity());
            for level in [&mut current, &mut *next] {
                let (added, borrowed) = fill_negative(&ops, ds, g, level);
                self.water.fixer += 0.5 * added;
                self.water.borrowed += 0.5 * borrowed;
            }
        }
        let old_previous = std::mem::replace(&mut self.previous, current);
        self.current = std::mem::replace(&mut self.next, old_previous);
        self.model.advance_forcing(&self.current, dt);

        self.steps += 1;
        for r in reports {
            self.stats.solves += 1;
            self.stats.iterations += u64::from(r.iterations);
            self.stats.worst_residual = self.stats.worst_residual.max(r.relative_residual);
        }
        reports
    }
}

/// RAW filter (Williams 2009) on every prognostic field.
fn raw_filter(
    old: &TwoLayerState,
    now: &mut TwoLayerState,
    next: &mut TwoLayerState,
    nu: f64,
    alpha: f64,
) {
    let apply = |old: &[f64], now: &mut [f64], next: &mut [f64]| {
        now.par_iter_mut()
            .with_min_len(MIN_CHUNK)
            .zip(next.par_iter_mut().with_min_len(MIN_CHUNK))
            .zip(old.par_iter())
            .for_each(|((now, next), &old)| {
                let d = 0.5 * nu * (old - 2.0 * *now + *next);
                *now += alpha * d;
                *next -= (1.0 - alpha) * d;
            });
    };
    apply(&old.ps, &mut now.ps, &mut next.ps);
    for k in 0..2 {
        apply(&old.u[k], &mut now.u[k], &mut next.u[k]);
        apply(
            &old.theta_ps[k],
            &mut now.theta_ps[k],
            &mut next.theta_ps[k],
        );
        apply(&old.q_ps[k], &mut now.q_ps[k], &mut next.q_ps[k]);
    }
}

/// Fill negative humidity in `state`: from the same column's other layer
/// where it has water to spare, else with zero. Returns the water added and
/// the water moved between layers (kg).
fn fill_negative(
    ops: &Operators,
    ds: [f64; 2],
    gravity: f64,
    state: &mut TwoLayerState,
) -> (f64, f64) {
    let areas = &ops.geometry().cell_areas;
    let scale = ops.radius() * ops.radius() / gravity;
    let [q0, q1] = &mut state.q_ps;
    let (mut added, mut borrowed) = (0.0, 0.0);
    for i in 0..q0.len() {
        // Column contents δσₖ Qₖ (Pa).
        let mut c = [ds[0] * q0[i], ds[1] * q1[i]];
        for k in 0..2 {
            if c[k] < 0.0 {
                let take = (-c[k]).min(c[1 - k].max(0.0));
                c[1 - k] -= take;
                c[k] += take;
                borrowed += take * areas[i] * scale;
                if c[k] < 0.0 {
                    added -= c[k] * areas[i] * scale;
                    c[k] = 0.0;
                }
            }
        }
        q0[i] = c[0] / ds[0];
        q1[i] = c[1] / ds[1];
    }
    (added, borrowed)
}

/// Area-weighted mean surface pressure and mass-weighted mean θ per layer.
fn global_means(ops: &Operators, state: &TwoLayerState, _dsigma: [f64; 2]) -> (f64, [f64; 2]) {
    let areas = &ops.geometry().cell_areas;
    let total_area: f64 = areas.iter().sum();
    let mass = dot(areas, &state.ps);
    let p_bar = mass / total_area;
    let th = [0, 1].map(|k| dot(areas, &state.theta_ps[k]) / mass);
    (p_bar, th)
}
