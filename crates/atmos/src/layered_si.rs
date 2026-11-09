//! Semi-implicit leapfrog time stepping for the layered model, with the
//! RAW filter: the layered version of [`crate::semi_implicit`].
//!
//! # Splitting
//!
//! About a resting, flat reference state (`p̄` and `θ̄ₖ`; see [Reference
//! state](#reference-state)), the gravity-wave terms are linear in
//! `y = (p_s, θ₀ᴸ, …, θ_{N−1}ᴸ)` with `θₖᴸ = (Θₖ − θ̄ₖ p_s)/p̄`, and the layer
//! divergences `δₖ = ∇·uₖ`:
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
//! `A` has eigenvalues `−c²` for the external and `N − 1` internal gravity
//! waves, so in its eigenvectors the system splits into `N` scalar Helmholtz
//! problems
//! `(I − Δt² c² ∇²) M̃ = r̃`, each solved by conjugate gradient. Then
//! `uⁿ⁺¹ = U − Δt ∇Mⁿ⁺¹` and `yⁿ⁺¹ = Y + Δt B ∇·uⁿ⁺¹`.
//!
//! Because the new surface pressure and `Θ` come from the divergence of the
//! new winds, which has no net source, total mass and total `Θ` are conserved
//! to round-off whatever the solver tolerance.
//!
//! # Reference state
//!
//! The split is stable only if the reference's gravity waves are at least as
//! fast as the atmosphere's: otherwise the excess is stepped explicitly and
//! grows at long steps (Simmons, Hoskins & Burridge 1978). So by default the
//! reference is isothermal (more stable than any troposphere) and warm: at
//! [`REFERENCE_TEMPERATURE`], or at the starting state's warmest layer if
//! that is warmer, with `p̄` the global-mean surface pressure. The starting
//! state's own mean profile would do for two layers, but with more, a
//! weakly stratified start (Held–Suarez begins with the layers' mean
//! equilibrium θ only a few kelvin apart near the ground) blows up within a
//! month as the forcing stratifies it.
//!
//! # Diffusion
//!
//! The model's hyperdiffusion ([`crate::diffusion`]) and forcing
//! ([`crate::forcing`]) are evaluated at the old level and added to `N`:
//! `xⁿ⁺¹ = xⁿ⁻¹ + 2Δt D(xⁿ⁻¹) + …`. Centred in time, damping would be
//! unstable under leapfrog; lagged, it is stable while `Δt` times the fastest
//! damping rate is below 1. The energy each adds or removes is accumulated in
//! [`LayeredSemiImplicit::energy_budget`].
//!
//! # Humidity
//!
//! A moist state's humidity is carried by flux-corrected transport
//! ([`crate::transport`]) with the mass fluxes the step used, so it stays
//! positive and consistent with the surface pressure; the physics sources
//! (lagged, like the rest of the forcing) go in with it. The RAW filter acts
//! on `Q` too, and can leave a trace of negative humidity, which is filled
//! from the column's other layers where they can spare it (the nearest
//! first) and otherwise set to zero.
//! [`LayeredSemiImplicit::water_budget`] accounts for every kilogram:
//! because leapfrog advances two interleaved levels, the water it closes on
//! is the mean of the two current levels ([`LayeredSemiImplicit::water`]).

use fms_grid::Operators;
use fms_grid::operators::MIN_CHUNK;
use rayon::prelude::*;

use crate::forcing::EnergyTerms;
use crate::helmholtz::{Helmholtz, SolveReport, dot};
use crate::layered::{LayeredModel, LayeredState, Linearization};
use crate::levels::Eigen;
use crate::safety::{self, BlowUp, Health, Limits};
use crate::semi_implicit::SolverStats;
use crate::transport::{Transport, TransportStep};

/// Lowest temperature (K) of the default isothermal reference state (see
/// the module docs).
pub const REFERENCE_TEMPERATURE: f64 = 300.0;

/// Settings for [`LayeredSemiImplicit`].
#[derive(Clone, Debug, PartialEq)]
pub struct LayeredSemiImplicitOptions {
    /// Reference `(p̄, θ̄ₖ)` for the linearization; `None` uses an isothermal
    /// reference (see the module docs).
    pub reference: Option<(f64, Vec<f64>)>,
    /// RAW filter strength `ν`.
    pub filter_strength: f64,
    /// RAW filter weight `α`.
    pub filter_alpha: f64,
    /// Helmholtz solver tolerance (relative residual).
    pub tolerance: f64,
}

impl Default for LayeredSemiImplicitOptions {
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
/// ([`LayeredSemiImplicit::water`]).
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

/// The vertical modes: `A = V diag(−c²) V⁻¹`, fastest (external) first.
#[derive(Clone, Debug)]
struct Modes {
    /// Columns are the eigenvectors.
    v: Vec<Vec<f64>>,
    v_inv: Vec<Vec<f64>>,
    /// `c²` per mode.
    speed2: Vec<f64>,
}

impl Modes {
    fn of(a: &[Vec<f64>]) -> Self {
        // Eigenvalues increase: the most negative (fastest wave) first.
        let e = Eigen::of(a);
        Self {
            v: e.vectors,
            v_inv: e.inverse,
            speed2: e.values.iter().map(|l| (-l).max(0.0)).collect(),
        }
    }
}

/// The layered model stepped semi-implicitly.
pub struct LayeredSemiImplicit<'a> {
    model: LayeredModel<'a>,
    ops: Operators<'a>,
    helmholtz: Vec<Helmholtz<'a>>,
    lin: Linearization,
    /// `A = J B`.
    a: Vec<Vec<f64>>,
    modes: Modes,
    reference: (f64, Vec<f64>),
    options: LayeredSemiImplicitOptions,
    dt: f64,
    previous: LayeredState,
    current: LayeredState,
    steps: u64,
    stats: SolverStats,
    budget: EnergyBudget,
    water: WaterBudget,
    transport: Option<Transport>,
    // Scratch.
    tendency: LayeredState,
    /// Lagged tendencies: diffusion plus forcing.
    physics: LayeredState,
    forcing: LayeredState,
    next: LayeredState,
    big_u: Vec<Vec<f64>>,
    m: Vec<Vec<f64>>,
    modal: Vec<Vec<f64>>,
    rhs: Vec<Vec<f64>>,
    div: Vec<Vec<f64>>,
    div_prev: Vec<Vec<f64>>,
    grad: Vec<Vec<f64>>,
    cells: Vec<Vec<f64>>,
    mass_flux: Vec<Vec<f64>>,
    zeros: Vec<f64>,
}

impl<'a> LayeredSemiImplicit<'a> {
    /// Start from `state` with step `dt` (s). Takes the first step at once,
    /// with RK4 substeps at Courant number `start_courant` (e.g. 0.5).
    ///
    /// # Panics
    /// If the model has hyperdiffusion too strong for `dt`:
    /// `dt · max_rate ≥ 1` ([`crate::Hyperdiffusion::max_rate`]).
    pub fn new(
        mut model: LayeredModel<'a>,
        state: LayeredState,
        dt: f64,
        options: LayeredSemiImplicitOptions,
        start_courant: f64,
    ) -> Self {
        if let Some(d) = model.hyperdiffusion() {
            assert!(
                dt * d.max_rate() < 1.0,
                "lagged hyperdiffusion is unstable at this step: Δt × fastest damping rate = {:.2} (must be < 1); lengthen the e-folding time",
                dt * d.max_rate()
            );
        }
        let n = model.layers();
        assert_eq!(state.layers(), n, "the state must have the model's layers");
        let ops = *model.operators();
        let reference = options
            .reference
            .clone()
            .unwrap_or_else(|| isothermal_reference(&model, &state));
        let lin = model.linearization(reference.0, &reference.1);
        let a = model.gravity_wave_matrix(reference.0, &reference.1);
        let modes = Modes::of(&a);
        let helmholtz = modes
            .speed2
            .iter()
            .map(|&c2| {
                let mut h = Helmholtz::new(ops, dt * dt * c2);
                h.tolerance = options.tolerance;
                h
            })
            .collect();

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
        let cells = || vec![vec![0.0; nc]; n];
        let edges = || vec![vec![0.0; ne]; n];
        let moist = state.is_moist();
        if moist {
            // The RK4 start's centred transport isn't monotone: keep its
            // humidity within the initial range nearby (the cell, its
            // neighbours and the adjacent layers). What this changes is part
            // of the starting water, not the budget.
            let mesh = ops.mesh();
            for k in 0..n {
                let q = |l: usize, i: usize| state.q_ps[l][i] / state.ps[i];
                for i in 0..nc {
                    let (mut lo, mut hi) = (q(k, i), q(k, i));
                    for l in adjacent(k, n) {
                        lo = lo.min(q(l, i));
                        hi = hi.max(q(l, i));
                    }
                    for &c in mesh.cell_neighbors(i as u32) {
                        lo = lo.min(q(k, c as usize));
                        hi = hi.max(q(k, c as usize));
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
            transport: moist.then(|| Transport::new(n, nc, ne)),
            tendency: state.dry(),
            physics,
            forcing,
            next: state,
            big_u: edges(),
            m: cells(),
            modal: cells(),
            rhs: cells(),
            div: cells(),
            div_prev: cells(),
            grad: edges(),
            cells: cells(),
            mass_flux: edges(),
            zeros: vec![0.0; nc],
        };
        // First guess for the modal potentials: those of the current state.
        let y = solver.y_of(&solver.current.clone());
        let m = solver.potentials(&y);
        solver.modal = solver.to_modes(&m);
        solver
    }

    /// The current state.
    pub fn state(&self) -> &LayeredState {
        &self.current
    }

    /// The current state, mutably: for experiments such as injecting a
    /// fault. Only the current time level changes, not the previous one.
    pub fn state_mut(&mut self) -> &mut LayeredState {
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
        let (p_bar, th_bar) = &self.reference;
        Limits::around(self.courant_limit(), *p_bar, th_bar)
    }

    /// Check the current state against [`limits`](Self::limits).
    pub fn check(&self) -> Result<Health, BlowUp> {
        safety::check(&self.current, &self.limits())
    }

    /// The underlying model.
    pub fn model(&self) -> &LayeredModel<'a> {
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

    /// Reference `(p̄, θ̄ₖ)` of the linearization.
    pub fn reference(&self) -> (f64, &[f64]) {
        (self.reference.0, &self.reference.1)
    }

    /// Gravity-wave speeds of the reference state (m s⁻¹), external first.
    pub fn wave_speeds(&self) -> Vec<f64> {
        self.modes.speed2.iter().map(|c| c.sqrt()).collect()
    }

    /// Helmholtz solver totals (all modes).
    pub fn solver_stats(&self) -> SolverStats {
        self.stats
    }

    /// The longest stable step (s) for winds up to `max_wind` (m s⁻¹) at
    /// Courant number `courant`: advection and the inertial limit, as for the
    /// shallow-water scheme ([`crate::SemiImplicit::stable_dt`]).
    pub fn stable_dt(model: &LayeredModel, max_wind: f64, courant: f64) -> f64 {
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

    /// `y = (p_s, θ₀ᴸ, …)` of a state.
    fn y_of(&self, s: &LayeredState) -> Vec<Vec<f64>> {
        let (p_bar, th_bar) = (self.reference.0, &self.reference.1);
        let theta_l = |k: usize| -> Vec<f64> {
            s.theta_ps[k]
                .par_iter()
                .with_min_len(MIN_CHUNK)
                .zip(&s.ps)
                .map(|(t, p)| (t - th_bar[k] * p) / p_bar)
                .collect()
        };
        let mut y = vec![s.ps.clone()];
        y.extend((0..s.layers()).map(theta_l));
        y
    }

    /// `M = J y` per layer.
    fn potentials(&self, y: &[Vec<f64>]) -> Vec<Vec<f64>> {
        self.lin.j.par_iter().map(|row| combine(row, y)).collect()
    }

    fn to_modes(&self, m: &[Vec<f64>]) -> Vec<Vec<f64>> {
        self.modes
            .v_inv
            .par_iter()
            .map(|row| combine(row, m))
            .collect()
    }

    /// Advance one step. Returns the solve report of each vertical mode
    /// (external first).
    // Loops over the layers index several per-layer arrays together.
    #[allow(clippy::needless_range_loop)]
    pub fn step(&mut self) -> Vec<SolveReport> {
        let ops = self.ops;
        let dt = self.dt;
        let n = self.model.layers();
        let (p_bar, th_bar) = (self.reference.0, self.reference.1.clone());
        let (j, b, a) = (&self.lin.j, &self.lin.b, &self.a);

        // Full tendency at the middle level.
        let empty = LayeredState::new(Vec::new(), Vec::new(), Vec::new());
        let mut current = std::mem::replace(&mut self.current, empty);
        // Hyperdiffusion and forcing are lagged to the old level.
        let diffuse = self.model.hyperdiffusion().is_some();
        let forced = self.model.has_forcing();
        let (diffusion_rate, forcing_rate) = self.model.step_tendencies(
            &current,
            &self.previous,
            &mut self.tendency,
            &mut self.physics,
            &mut self.forcing,
        );
        if let Some(rate) = diffusion_rate {
            self.budget.diffusion.add_scaled(dt, rate);
        }
        if forced {
            if let Some(rate) = forcing_rate {
                self.budget.forcing.add_scaled(dt, rate);
            }
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
        {
            let (now, old) = (&current.u, &self.previous.u);
            self.div
                .par_iter_mut()
                .zip(self.div_prev.par_iter_mut())
                .enumerate()
                .for_each(|(k, (d, d_prev))| {
                    ops.divergence(&now[k], d);
                    ops.divergence(&old[k], d_prev);
                });
        }
        let (div, div_prev) = (&self.div, &self.div_prev);

        // Y = yⁿ⁻¹ + 2Δt (T_y − B δⁿ) + Δt B δⁿ⁻¹, in y-space.
        let y_prev = self.y_of(&self.previous);
        let nc = current.ps.len();
        let mut big_y = vec![vec![0.0; nc]; n + 1];
        big_y.par_iter_mut().enumerate().for_each(|(r, out)| {
            let th_bar = &th_bar;
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
                    let (mut b_now, mut b_old) = (b[r][0] * div[0][i], b[r][0] * div_prev[0][i]);
                    for l in 1..n {
                        b_now += b[r][l] * div[l][i];
                        b_old += b[r][l] * div_prev[l][i];
                    }
                    *o = y_prev[r][i] + 2.0 * dt * (t_y - b_now) + dt * b_old;
                });
        });

        // Uₖ = uₖⁿ⁻¹ + 2Δt (T_u + ∇Mₖ(yⁿ)) − Δt ∇Mₖ(yⁿ⁻¹).
        let y_now = self.y_of(&current);
        let (m_now, m_prev) = (self.potentials(&y_now), self.potentials(&y_prev));
        let mut div_u = vec![vec![0.0; nc]; n];
        {
            let (m_now, m_prev, previous) = (&m_now, &m_prev, &self.previous);
            self.big_u
                .par_iter_mut()
                .zip(self.cells.par_iter_mut())
                .zip(self.grad.par_iter_mut())
                .zip(div_u.par_iter_mut())
                .enumerate()
                .for_each(|(k, (((big_u, cells), grad), div_u))| {
                    // ∇(2 Mₖ(yⁿ) − Mₖ(yⁿ⁻¹)) combines both linear terms in one gradient.
                    cells
                        .par_iter_mut()
                        .with_min_len(MIN_CHUNK)
                        .enumerate()
                        .for_each(|(i, c)| *c = 2.0 * m_now[k][i] - m_prev[k][i]);
                    ops.gradient(cells, grad);
                    let (prev_u, t_u, d_u, grad) = (&previous.u[k], &t.u[k], &d.u[k], &*grad);
                    big_u
                        .par_iter_mut()
                        .with_min_len(MIN_CHUNK)
                        .enumerate()
                        .for_each(|(e, u)| {
                            let diffusion = if diffuse { d_u[e] } else { 0.0 };
                            *u = prev_u[e] + 2.0 * dt * (t_u[e] + diffusion) + dt * grad[e];
                        });
                    ops.divergence(big_u, div_u);
                });
        }

        // rhs = J Y + Δt A ∇·U.
        self.rhs.par_iter_mut().enumerate().for_each(|(k, out)| {
            let (big_y, div_u) = (&big_y, &div_u);
            out.par_iter_mut()
                .with_min_len(MIN_CHUNK)
                .enumerate()
                .for_each(|(i, r)| {
                    let mut jy = j[k][0] * big_y[0][i];
                    for c in 1..=n {
                        jy += j[k][c] * big_y[c][i];
                    }
                    let mut ad = a[k][0] * div_u[0][i];
                    for c in 1..n {
                        ad += a[k][c] * div_u[c][i];
                    }
                    *r = jy + dt * ad;
                });
        });

        // N scalar Helmholtz problems in the vertical modes.
        let rhs_modal = self.to_modes(&self.rhs);
        let reports: Vec<SolveReport> = self
            .helmholtz
            .par_iter_mut()
            .zip(self.modal.par_iter_mut())
            .zip(rhs_modal.par_iter())
            .map(|((solver, x), b)| solver.solve(b, x))
            .collect();
        {
            let (v, modal) = (&self.modes.v, &self.modal);
            self.m
                .par_iter_mut()
                .enumerate()
                .for_each(|(k, out)| *out = combine(&v[k], modal));
        }

        // uⁿ⁺¹ = U − Δt ∇Mⁿ⁺¹; yⁿ⁺¹ = Y + Δt B ∇·uⁿ⁺¹.
        let next = &mut self.next;
        {
            let (m, big_u) = (&self.m, &self.big_u);
            next.u
                .par_iter_mut()
                .zip(self.grad.par_iter_mut())
                .zip(self.div.par_iter_mut())
                .enumerate()
                .for_each(|(k, ((u, grad), div))| {
                    ops.gradient(&m[k], grad);
                    let (big_u, grad) = (&big_u[k], &*grad);
                    u.par_iter_mut()
                        .with_min_len(MIN_CHUNK)
                        .enumerate()
                        .for_each(|(e, u)| *u = big_u[e] - dt * grad[e]);
                    ops.divergence(u, div);
                });
        }
        let div = &self.div;
        let y_new: Vec<Vec<f64>> = (0..=n)
            .into_par_iter()
            .map(|r| {
                let big_y = &big_y;
                (0..nc)
                    .into_par_iter()
                    .with_min_len(MIN_CHUNK)
                    .map(|i| {
                        let mut bd = b[r][0] * div[0][i];
                        for l in 1..n {
                            bd += b[r][l] * div[l][i];
                        }
                        big_y[r][i] + dt * bd
                    })
                    .collect()
            })
            .collect();
        next.ps.copy_from_slice(&y_new[0]);
        for k in 0..n {
            let (ps, yk, th) = (&y_new[0], &y_new[k + 1], th_bar[k]);
            next.theta_ps[k]
                .par_iter_mut()
                .with_min_len(MIN_CHUNK)
                .enumerate()
                .for_each(|(i, t)| *t = p_bar * yk[i] + th * ps[i]);
        }

        // Humidity: flux-corrected transport with the step's mass fluxes,
        // F̄ₖ = Fₖⁿ + p̄ δσₖ (uₖⁿ⁻¹ + uₖⁿ⁺¹ − 2uₖⁿ)/2.
        if let Some(transport) = &mut self.transport {
            let ds = self.model.dsigma();
            ops.cell_to_edge(&current.ps, &mut self.grad[0]);
            let mut flux_div = vec![vec![0.0; nc]; n];
            {
                let (ps_e, previous, next) = (&self.grad[0], &self.previous, &*next);
                self.mass_flux
                    .par_iter_mut()
                    .zip(flux_div.par_iter_mut())
                    .enumerate()
                    .for_each(|(k, (flux, div))| {
                        let (now, old, new) = (&current.u[k], &previous.u[k], &next.u[k]);
                        flux.par_iter_mut()
                            .with_min_len(MIN_CHUNK)
                            .enumerate()
                            .for_each(|(e, f)| {
                                *f = ds[k]
                                    * (ps_e[e] * now[e]
                                        + 0.5 * p_bar * (old[e] + new[e] - 2.0 * now[e]));
                            });
                        ops.divergence(flux, div);
                    });
            }
            // The interface fluxes they imply (as in the dynamics).
            let interface: Vec<Vec<f64>> = self
                .model
                .levels()
                .interior()
                .iter()
                .enumerate()
                .map(|(j, &s)| {
                    let flux_div = &flux_div;
                    (0..nc)
                        .into_par_iter()
                        .with_min_len(MIN_CHUNK)
                        .map(|i| {
                            let mut above = flux_div[0][i];
                            for d in &flux_div[1..=j] {
                                above += d[i];
                            }
                            let mut below = flux_div[j + 1][i];
                            for d in &flux_div[j + 2..] {
                                below += d[i];
                            }
                            s * below - (1.0 - s) * above
                        })
                        .collect()
                })
                .collect();
            let q_mid: Vec<Vec<f64>> = (0..n).map(|k| current.humidity(k)).collect();
            let ratio = self.model.interface_ratios(&current);
            let source: Vec<&[f64]> = if diffuse {
                d.q_ps.iter().map(Vec::as_slice).collect()
            } else {
                vec![self.zeros.as_slice(); n]
            };
            let input = TransportStep {
                dsigma: ds,
                interval: 2.0 * dt,
                ps_old: &self.previous.ps,
                ps_new: &next.ps,
                flux: self.mass_flux.iter().map(Vec::as_slice).collect(),
                interface: interface.iter().map(Vec::as_slice).collect(),
                q_old: self.previous.q_ps.iter().map(Vec::as_slice).collect(),
                q_mid: q_mid.iter().map(Vec::as_slice).collect(),
                source,
                interface_ratio: ratio
                    .iter()
                    .map(|[up, down]| [up.as_slice(), down.as_slice()])
                    .collect(),
            };
            let outs: Vec<&mut [f64]> = next.q_ps.iter_mut().map(Vec::as_mut_slice).collect();
            let report = transport.step(&ops, &input, outs);
            self.water.max_substeps = self.water.max_substeps.max(report.substeps);
        }

        // RAW filter on xⁿ and xⁿ⁺¹, then rotate the three levels (the old
        // xⁿ⁻¹ buffer becomes scratch for the next step).
        let (nu, alpha) = (self.options.filter_strength, self.options.filter_alpha);
        // Humidity's range over the three levels, per cell and layer, and the
        // water of the two levels the filter changes.
        let moist = self.transport.is_some();
        let (bounds, water_before) = if moist {
            let bounds: Vec<Vec<(f64, f64)>> = (0..n)
                .map(|k| {
                    let (previous, current, next) = (&self.previous, &current, &*next);
                    (0..nc)
                        .into_par_iter()
                        .with_min_len(MIN_CHUNK)
                        .map(|i| {
                            let q = [
                                previous.q_ps[k][i] / previous.ps[i],
                                current.q_ps[k][i] / current.ps[i],
                                next.q_ps[k][i] / next.ps[i],
                            ];
                            (q[0].min(q[1]).min(q[2]), q[0].max(q[1]).max(q[2]))
                        })
                        .collect()
                })
                .collect();
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
                for k in 0..n {
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
            let g = self.model.gravity();
            let ds = self.model.dsigma();
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
        for r in &reports {
            self.stats.solves += 1;
            self.stats.iterations += u64::from(r.iterations);
            self.stats.worst_residual = self.stats.worst_residual.max(r.relative_residual);
        }
        reports
    }
}

/// `Σᵢ cᵢ xᵢ` per cell, left to right.
fn combine(c: &[f64], x: &[Vec<f64>]) -> Vec<f64> {
    (0..x[0].len())
        .into_par_iter()
        .with_min_len(MIN_CHUNK)
        .map(|i| {
            let mut s = c[0] * x[0][i];
            for (c, x) in c.iter().zip(x).skip(1) {
                s += c * x[i];
            }
            s
        })
        .collect()
}

/// The layers next to `k` of `n`.
fn adjacent(k: usize, n: usize) -> impl Iterator<Item = usize> {
    [k.checked_sub(1), (k + 1 < n).then_some(k + 1)]
        .into_iter()
        .flatten()
}

/// RAW filter (Williams 2009) on every prognostic field.
fn raw_filter(
    old: &LayeredState,
    now: &mut LayeredState,
    next: &mut LayeredState,
    nu: f64,
    alpha: f64,
) {
    for ((old, now), next) in old
        .parts()
        .into_iter()
        .zip(now.parts_mut())
        .zip(next.parts_mut())
    {
        now.par_iter_mut()
            .with_min_len(MIN_CHUNK)
            .zip(next.par_iter_mut().with_min_len(MIN_CHUNK))
            .zip(old.par_iter())
            .for_each(|((now, next), &old)| {
                let d = 0.5 * nu * (old - 2.0 * *now + *next);
                *now += alpha * d;
                *next -= (1.0 - alpha) * d;
            });
    }
}

/// Fill negative humidity in `state`: from the same column's other layers
/// where they have water to spare (the nearest first, the one below before
/// the one above), else with zero. Returns the water added and the water
/// moved between layers (kg).
// Cells and layers index several arrays together.
#[allow(clippy::needless_range_loop)]
fn fill_negative(
    ops: &Operators,
    ds: &[f64],
    gravity: f64,
    state: &mut LayeredState,
) -> (f64, f64) {
    let areas = &ops.geometry().cell_areas;
    let scale = ops.radius() * ops.radius() / gravity;
    let n = ds.len();
    let (mut added, mut borrowed) = (0.0, 0.0);
    let mut c = vec![0.0; n];
    for i in 0..state.ps.len() {
        // Column contents δσₖ Qₖ (Pa).
        for k in 0..n {
            c[k] = ds[k] * state.q_ps[k][i];
        }
        for k in 0..n {
            if c[k] >= 0.0 {
                continue;
            }
            for distance in 1..n {
                let donors = [Some(k + distance), k.checked_sub(distance)];
                for l in donors.into_iter().flatten().filter(|&l| l < n) {
                    if c[k] < 0.0 {
                        let take = (-c[k]).min(c[l].max(0.0));
                        c[l] -= take;
                        c[k] += take;
                        borrowed += take * areas[i] * scale;
                    }
                }
            }
            if c[k] < 0.0 {
                added -= c[k] * areas[i] * scale;
                c[k] = 0.0;
            }
        }
        for k in 0..n {
            state.q_ps[k][i] = c[k] / ds[k];
        }
    }
    (added, borrowed)
}

/// The default reference: global-mean surface pressure, and an isothermal
/// atmosphere at [`REFERENCE_TEMPERATURE`] or the warmest layer of `state`,
/// whichever is warmer.
fn isothermal_reference(model: &LayeredModel, state: &LayeredState) -> (f64, Vec<f64>) {
    let ops = model.operators();
    let areas = &ops.geometry().cell_areas;
    let p_bar = dot(areas, &state.ps) / areas.iter().sum::<f64>();
    let mut warmest = REFERENCE_TEMPERATURE;
    for (k, theta_ps) in state.theta_ps.iter().enumerate() {
        for (tp, &ps) in theta_ps.iter().zip(&state.ps) {
            warmest = warmest.max(tp / ps * model.layer_exner(ps)[k]);
        }
    }
    let exner = model.layer_exner(p_bar);
    (p_bar, exner.iter().map(|e| warmest / e).collect())
}

/// The starting state's global means: area-weighted surface pressure and
/// mass-weighted θ per layer (the default reference before M2 step 9.5; a
/// stable choice for two layers only).
pub fn global_means(ops: &Operators, state: &LayeredState) -> (f64, Vec<f64>) {
    let areas = &ops.geometry().cell_areas;
    let total_area: f64 = areas.iter().sum();
    let mass = dot(areas, &state.ps);
    let p_bar = mass / total_area;
    let th = state
        .theta_ps
        .iter()
        .map(|t| dot(areas, t) / mass)
        .collect();
    (p_bar, th)
}
