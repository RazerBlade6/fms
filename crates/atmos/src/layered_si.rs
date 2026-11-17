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
//!
//! # Restarts
//!
//! [`LayeredSemiImplicit::save_restart`] saves what the stepper carries
//! between steps, and its forcing's state after it
//! ([`Forcing::save_restart`](crate::Forcing::save_restart)), and
//! [`LayeredSemiImplicit::load_restart`] loads both over a freshly built
//! stepper's; with the same settings the run continues bit for bit.

use std::time::{Duration, Instant};

use fms_grid::Operators;
use fms_grid::operators::{MIN_CHUNK, for_each_mut, for_each_mut2};
use fms_restart::{Restart, RestartError, count};
use rayon::prelude::*;

use crate::forcing::EnergyTerms;
use crate::helmholtz::{Helmholtz, SolveReport, dot};
use crate::layered::{LayeredModel, LayeredState, Linearization};
use crate::levels::{Eigen, MAX_LAYERS};
use crate::safety::{self, BlowUp, Health, Limits};
use crate::semi_implicit::SolverStats;
use crate::timing::{Phase, Stopwatch, Timings};
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
    /// Time spent in the step's phases, and in each vertical mode's solves
    /// (not saved in restarts).
    timings: Timings,
    solve_time: Vec<Duration>,
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
            timings: Timings::default(),
            solve_time: vec![Duration::ZERO; n],
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

    /// Time spent in each phase of the steps since this stepper was built
    /// (not the steps before a restart), the forcing's own phases included
    /// ([`Forcing::timings`](crate::Forcing::timings)), which count from
    /// the forcing's construction and include the RK4 start: take
    /// differences ([`Timings::since`]).
    pub fn timings(&self) -> Timings {
        let mut t = self.timings;
        if let Some(f) = self.model.forcing() {
            t.merge(&f.timings());
        }
        t
    }

    /// Time spent in each vertical mode's Helmholtz solves since this
    /// stepper was built, external first. The modes are solved side by
    /// side, so these overlap.
    pub fn solve_times(&self) -> &[Duration] {
        &self.solve_time
    }

    /// A restart holding everything the stepper carries between steps,
    /// under `dynamics.`: both leapfrog levels, the modal potentials (the
    /// solver's first guess, which the next step's bits depend on), the
    /// energy and water budgets, the solver totals and the step count, with
    /// the σ interfaces to check against on loading and the time step to
    /// re-base the step count by; then the forcing's own state
    /// ([`Forcing::save_restart`](crate::Forcing::save_restart)).
    ///
    /// Settings are not loaded: the time step, the reference state and the
    /// operators built from it, the filter and solver options all come from
    /// the stepper that loads it.
    pub fn save_restart(&self) -> Restart {
        let mut r = Restart::new(self.ops.mesh().fingerprint(), self.steps);
        r.put(SIGMA, self.model.levels().interfaces());
        for (prefix, level) in [(CURRENT, &self.current), (PREVIOUS, &self.previous)] {
            for (name, values) in level_arrays(prefix, level) {
                r.put(&name, values);
            }
        }
        for (k, m) in self.modal.iter().enumerate() {
            r.put(&format!("{MODAL}{k}"), m);
        }
        let (e, w, s) = (&self.budget, &self.water, &self.stats);
        r.put(
            ENERGY_BUDGET,
            &[
                e.diffusion.kinetic,
                e.diffusion.thermal,
                e.forcing.kinetic,
                e.forcing.thermal,
            ],
        );
        r.put(
            WATER_BUDGET,
            &[
                w.sources,
                w.filter,
                w.fixer,
                w.borrowed,
                f64::from(w.max_substeps),
            ],
        );
        r.put(
            SOLVER_STATS,
            &[s.solves as f64, s.iterations as f64, s.worst_residual],
        );
        r.put(DT, &[self.dt]);
        if let Some(forcing) = self.model.forcing() {
            forcing.save_restart(&mut r);
        }
        r
    }

    /// Load a state saved by [`save_restart`](Self::save_restart) over this
    /// stepper's, including the step count. All or nothing: on an error the
    /// stepper is unchanged.
    ///
    /// Refused if the file is from another mesh, other σ levels, or a moist
    /// state into a dry stepper (or the reverse), or if the forcing refuses
    /// its part ([`Forcing::load_restart`](crate::Forcing::load_restart)),
    /// which it loads only once the stepper's part has been checked. The
    /// reference state stays the one this stepper built from its own
    /// starting state, not the file's, so the continuation is bit for bit
    /// only when this stepper was built as the saved one was (same model,
    /// starting state, time step and options); otherwise it is a warm start.
    ///
    /// Returns a sentence for each way the state was adapted to this
    /// stepper's settings, its own then the forcing's; none means nothing
    /// was. With another time step the two leapfrog levels are loaded as
    /// they are and the step count is re-based so that steps × step stays
    /// the model time.
    pub fn load_restart(&mut self, from: &Restart) -> Result<Vec<String>, RestartError> {
        let run = self.ops.mesh().fingerprint();
        if from.grid() != run {
            return Err(RestartError::DifferentGrid {
                file: from.grid(),
                run,
            });
        }
        let sigma = self.model.levels().interfaces();
        let same_levels = match from.get(SIGMA, sigma.len()) {
            Ok(saved) => saved
                .iter()
                .zip(sigma)
                .all(|(a, b)| a.to_bits() == b.to_bits()),
            Err(RestartError::WrongLength { .. }) => false,
            Err(e) => return Err(e),
        };
        if !same_levels {
            return Err(RestartError::Mismatch(format!(
                "its σ levels differ from this run's {sigma:?} (a state on other levels would need vertical regridding)"
            )));
        }
        let moist = self.transport.is_some();
        if from.contains(&format!("{CURRENT}q_ps.0")) != moist {
            let (file, run) = if moist {
                ("dry", "moist")
            } else {
                ("moist", "dry")
            };
            return Err(RestartError::Mismatch(format!(
                "its atmosphere is {file} and this run's is {run}"
            )));
        }

        // Check everything before changing anything.
        let mut levels = Vec::new();
        for (prefix, level) in [(CURRENT, &self.current), (PREVIOUS, &self.previous)] {
            let saved = level_arrays(prefix, level)
                .into_iter()
                .map(|(name, values)| from.get(&name, values.len()))
                .collect::<Result<Vec<_>, _>>()?;
            levels.push(saved);
        }
        let modal = self
            .modal
            .iter()
            .enumerate()
            .map(|(k, m)| from.get(&format!("{MODAL}{k}"), m.len()))
            .collect::<Result<Vec<_>, _>>()?;
        let e = from.get(ENERGY_BUDGET, 4)?;
        let w = from.get(WATER_BUDGET, 5)?;
        let s = from.get(SOLVER_STATS, 3)?;
        let max_substeps = count(WATER_BUDGET, w[4])?;
        let (solves, iterations) = (count(SOLVER_STATS, s[0])?, count(SOLVER_STATS, s[1])?);
        let [saved_dt] = from.get_array(DT)?;
        let mut notes = Vec::new();
        let steps = if saved_dt.to_bits() == self.dt.to_bits() {
            from.step()
        } else {
            if !(saved_dt > 0.0 && saved_dt.is_finite()) {
                return Err(RestartError::Mismatch(format!(
                    "{DT} holds {saved_dt}, not a time step"
                )));
            }
            let steps = (from.step() as f64 * saved_dt / self.dt).round() as u64;
            notes.push(format!(
                "the time step changed from {saved_dt} s to {} s: both leapfrog levels are loaded as they are, \
                 and the {} steps taken are counted as {steps} of the new step, the nearest to the same model time",
                self.dt,
                from.step()
            ));
            steps
        };
        // The forcing checks its own part and loads it, all or nothing; the
        // copies below cannot fail.
        if let Some(forcing) = self.model.forcing_mut() {
            notes.extend(forcing.load_restart(from)?);
        }

        for (level, saved) in [&mut self.current, &mut self.previous]
            .into_iter()
            .zip(levels)
        {
            for (to, from) in level_arrays_mut(level).into_iter().zip(saved) {
                to.copy_from_slice(from);
            }
        }
        for (to, from) in self.modal.iter_mut().zip(modal) {
            to.copy_from_slice(from);
        }
        self.budget = EnergyBudget {
            diffusion: EnergyTerms {
                kinetic: e[0],
                thermal: e[1],
            },
            forcing: EnergyTerms {
                kinetic: e[2],
                thermal: e[3],
            },
        };
        self.water = WaterBudget {
            sources: w[0],
            filter: w[1],
            fixer: w[2],
            borrowed: w[3],
            max_substeps,
        };
        self.stats = SolverStats {
            solves,
            iterations,
            worst_residual: s[2],
        };
        self.steps = steps;
        Ok(notes)
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
        // On a pool thread, so the step's many small parallel loops start
        // without a hand-off from the caller's thread each time.
        let start = Instant::now();
        let reports = rayon::scope(|_| self.step_in_pool());
        self.timings.add(Phase::Step, start.elapsed());
        reports
    }

    // Loops over the layers index several per-layer arrays together.
    #[allow(clippy::needless_range_loop)]
    fn step_in_pool(&mut self) -> Vec<SolveReport> {
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
        // Laps partition the phases below; one that doesn't run leaves its
        // moment to the next.
        let mut clock = Stopwatch::start();
        let (diffusion_rate, forcing_rate, times) = self.model.step_tendencies(
            &current,
            &self.previous,
            &mut self.tendency,
            &mut self.physics,
            &mut self.forcing,
        );
        self.timings.add(Phase::Tendencies, clock.lap());
        if let Some(time) = times.physics {
            self.timings.add(Phase::Physics, time);
        }
        self.timings.add(Phase::Dynamics, times.dynamics);
        if let Some(time) = times.hyperdiffusion {
            self.timings.add(Phase::Hyperdiffusion, time);
        }
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
                    for_each_mut(p, |i, p| *p += f[i]);
                }
            } else {
                std::mem::swap(&mut self.physics, &mut self.forcing);
            }
        }
        let diffuse = diffuse || forced;
        let t = &self.tendency;
        let d = &self.physics;
        self.timings.add(Phase::BudgetsAndLagged, clock.lap());

        // Per layer: the divergences of uⁿ and uⁿ⁻¹, and the θᴸ rows of yⁿ⁻¹
        // and yⁿ (whose first row is p_s itself).
        let nc = current.ps.len();
        let th_bar = &th_bar;
        let mut y_prev = vec![vec![0.0; nc]; n];
        let mut y_now = vec![vec![0.0; nc]; n];
        {
            let (now, old) = (&current, &self.previous);
            let theta_l = |s: &LayeredState, k: usize, out: &mut [f64]| {
                for_each_mut(out, |i, y| {
                    *y = (s.theta_ps[k][i] - th_bar[k] * s.ps[i]) / p_bar
                });
            };
            self.div
                .par_iter_mut()
                .zip(self.div_prev.par_iter_mut())
                .zip(y_prev.par_iter_mut())
                .zip(y_now.par_iter_mut())
                .enumerate()
                .for_each(|(k, (((div, div_prev), y_prev), y_now))| {
                    ops.divergence(&now.u[k], div);
                    ops.divergence(&old.u[k], div_prev);
                    theta_l(old, k, y_prev);
                    theta_l(now, k, y_now);
                });
        }
        let (div, div_prev) = (&self.div, &self.div_prev);
        let y_prev = rows(&self.previous.ps, &y_prev);
        let y_now = rows(&current.ps, &y_now);
        let (y_prev, y_now) = (&y_prev[..=n], &y_now[..=n]);

        // Per layer, with every layer's divergences and y at hand:
        // Y = yⁿ⁻¹ + 2Δt (T_y − B δⁿ) + Δt B δⁿ⁻¹, in y-space (row k + 1, and
        // row 0 with layer 0), then
        // Uₖ = uₖⁿ⁻¹ + 2Δt (T_u + ∇Mₖ(yⁿ)) − Δt ∇Mₖ(yⁿ⁻¹) and its divergence.
        let mut big_y = vec![vec![0.0; nc]; n + 1];
        let mut div_u = vec![vec![0.0; nc]; n];
        {
            let y_row = |r: usize, out: &mut [f64]| {
                for_each_mut(out, |i, o| {
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
            };
            let (big_y_0, big_y_rows) = big_y.split_first_mut().expect("y has rows");
            let mut first: [Option<&mut Vec<f64>>; MAX_LAYERS] = std::array::from_fn(|_| None);
            first[0] = Some(big_y_0);
            let previous = &self.previous;
            big_y_rows
                .par_iter_mut()
                .zip(first[..n].par_iter_mut())
                .zip(self.big_u.par_iter_mut())
                .zip(self.cells.par_iter_mut())
                .zip(self.grad.par_iter_mut())
                .zip(div_u.par_iter_mut())
                .enumerate()
                .for_each(|(k, (((((big_y, big_y_0), big_u), cells), grad), div_u))| {
                    if let Some(big_y_0) = big_y_0.as_deref_mut() {
                        y_row(0, big_y_0);
                    }
                    y_row(k + 1, big_y);
                    // ∇(2 Mₖ(yⁿ) − Mₖ(yⁿ⁻¹)) combines both linear terms in one gradient.
                    for_each_mut(cells, |i, c| {
                        *c = 2.0 * combine_at(&j[k], y_now, i) - combine_at(&j[k], y_prev, i);
                    });
                    ops.gradient(cells, grad);
                    let (prev_u, t_u, d_u, grad) = (&previous.u[k], &t.u[k], &d.u[k], &*grad);
                    for_each_mut(big_u, |e, u| {
                        let diffusion = if diffuse { d_u[e] } else { 0.0 };
                        *u = prev_u[e] + 2.0 * dt * (t_u[e] + diffusion) + dt * grad[e];
                    });
                    ops.divergence(big_u, div_u);
                });
        }

        // rhs = J Y + Δt A ∇·U, per layer.
        {
            let (big_y, div_u) = (&big_y, &div_u);
            self.rhs.iter_mut().enumerate().for_each(|(k, out)| {
                for_each_mut(out, |i, r| {
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
        }
        self.timings.add(Phase::ImplicitRhs, clock.lap());

        // N scalar Helmholtz problems in the vertical modes, each mode's
        // right-hand side taken from rhs in its own task. Beside them, the
        // humidity's interface ratios for the transport below (they need
        // only the middle level, whose Exner functions the dynamics' tendency
        // left).
        let (helmholtz, modal, model) = (&mut self.helmholtz, &mut self.modal, &self.model);
        let (v_inv, rhs) = (&self.modes.v_inv, &self.rhs);
        let solve_time = &mut self.solve_time;
        let moist = self.transport.is_some();
        let mut rhs_modal = vec![vec![0.0; nc]; n];
        let (reports, ratio): (Vec<SolveReport>, _) = rayon::join(
            || {
                helmholtz
                    .par_iter_mut()
                    .zip(modal.par_iter_mut())
                    .zip(rhs_modal.par_iter_mut())
                    .zip(solve_time.par_iter_mut())
                    .enumerate()
                    .map(|(mode, (((solver, x), b), time))| {
                        combine_into(&v_inv[mode], rhs, b);
                        let start = Instant::now();
                        let report = solver.solve(b, x);
                        *time += start.elapsed();
                        report
                    })
                    .collect()
            },
            || {
                moist.then(|| {
                    let start = Instant::now();
                    (
                        model.interface_ratios_of_tendency_state(&current),
                        start.elapsed(),
                    )
                })
            },
        );
        self.timings.add(Phase::HelmholtzSolves, clock.lap());
        let ratio = ratio.map(|(ratio, time)| {
            self.timings.add(Phase::InterfaceRatios, time);
            ratio
        });

        // Per layer: Mₖⁿ⁺¹ from the modes, uₖⁿ⁺¹ = Uₖ − Δt ∇Mₖⁿ⁺¹ and its
        // divergence.
        let next = &mut self.next;
        {
            let (v, modal, big_u) = (&self.modes.v, &self.modal, &self.big_u);
            next.u
                .par_iter_mut()
                .zip(self.m.par_iter_mut())
                .zip(self.grad.par_iter_mut())
                .zip(self.div.par_iter_mut())
                .enumerate()
                .for_each(|(k, (((u, m), grad), div))| {
                    combine_into(&v[k], modal, m);
                    ops.gradient(m, grad);
                    let (big_u, grad) = (&big_u[k], &*grad);
                    for_each_mut(u, |e, u| *u = big_u[e] - dt * grad[e]);
                    ops.divergence(u, div);
                });
        }
        // yⁿ⁺¹ = Y + Δt B ∇·uⁿ⁺¹: cheap, and across the layers, so on this
        // thread.
        {
            let (div, big_y) = (&self.div, &big_y);
            let y_new = |r: usize, i: usize| {
                let mut bd = b[r][0] * div[0][i];
                for l in 1..n {
                    bd += b[r][l] * div[l][i];
                }
                big_y[r][i] + dt * bd
            };
            for_each_mut(&mut next.ps, |i, p| *p = y_new(0, i));
            let ps = &next.ps;
            for (k, theta_ps) in next.theta_ps.iter_mut().enumerate() {
                let th = th_bar[k];
                for_each_mut(theta_ps, |i, t| *t = p_bar * y_new(k + 1, i) + th * ps[i]);
            }
        }
        self.timings.add(Phase::BackSubstitution, clock.lap());

        // Humidity: flux-corrected transport with the step's mass fluxes,
        // F̄ₖ = Fₖⁿ + p̄ δσₖ (uₖⁿ⁻¹ + uₖⁿ⁺¹ − 2uₖⁿ)/2.
        if let Some(transport) = &mut self.transport {
            let ds = self.model.dsigma();
            ops.cell_to_edge(&current.ps, &mut self.grad[0]);
            let mut flux_div = vec![vec![0.0; nc]; n];
            let mut q_mid = vec![vec![0.0; nc]; n];
            {
                let (ps_e, previous, next, current) =
                    (&self.grad[0], &self.previous, &*next, &current);
                self.mass_flux
                    .par_iter_mut()
                    .zip(flux_div.par_iter_mut())
                    .zip(q_mid.par_iter_mut())
                    .enumerate()
                    .for_each(|(k, ((flux, div), q_mid))| {
                        let (now, old, new) = (&current.u[k], &previous.u[k], &next.u[k]);
                        for_each_mut(flux, |e, f| {
                            *f = ds[k]
                                * (ps_e[e] * now[e]
                                    + 0.5 * p_bar * (old[e] + new[e] - 2.0 * now[e]));
                        });
                        ops.divergence(flux, div);
                        let (q, ps) = (&current.q_ps[k], &current.ps);
                        for_each_mut(q_mid, |i, x| *x = q[i] / ps[i]);
                    });
            }
            // The interface fluxes they imply (as in the dynamics): cheap,
            // and across the layers, so on this thread.
            let mut interface = vec![vec![0.0; nc]; n - 1];
            for (j, (&s, out)) in self
                .model
                .levels()
                .interior()
                .iter()
                .zip(&mut interface)
                .enumerate()
            {
                let flux_div = &flux_div;
                for_each_mut(out, |i, x| {
                    let mut above = flux_div[0][i];
                    for d in &flux_div[1..=j] {
                        above += d[i];
                    }
                    let mut below = flux_div[j + 1][i];
                    for d in &flux_div[j + 2..] {
                        below += d[i];
                    }
                    *x = s * below - (1.0 - s) * above;
                });
            }
            let ratio = ratio.expect("a moist state has interface ratios");
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
            let fluxes = clock.lap();
            let report = transport.step(&ops, &input, outs);
            self.water.max_substeps = self.water.max_substeps.max(report.substeps);
            let fct = clock.lap();
            self.timings.add(Phase::MassFluxes, fluxes);
            self.timings.add(Phase::FluxCorrectedTransport, fct);
            self.timings.add(Phase::HumidityTransport, fluxes + fct);
        }

        // RAW filter on xⁿ and xⁿ⁺¹, then rotate the three levels (the old
        // xⁿ⁻¹ buffer becomes scratch for the next step).
        let (nu, alpha) = (self.options.filter_strength, self.options.filter_alpha);
        let moist = self.transport.is_some();
        let mut filtered = None;
        if moist {
            // Per layer: humidity's range over the three levels, the water of
            // the two levels the filter changes, the filter, and the filtered
            // humidity kept within that range (the filter is linear, so this
            // never acts on a uniform q) against the filtered surface
            // pressure, and the water after. Surface pressure is filtered
            // after the layers, which take its filtered value as the filter
            // gives it.
            let raw = |old: f64, now: f64, next: f64| {
                let d = 0.5 * nu * (old - 2.0 * now + next);
                (now + alpha * d, next - (1.0 - alpha) * d)
            };
            let mut bounds = vec![vec![(0.0, 0.0); nc]; n];
            let mut water = [[0.0; 4]; MAX_LAYERS];
            {
                let previous = &self.previous;
                let model = &self.model;
                let LayeredState {
                    ps: now_ps,
                    u: now_u,
                    theta_ps: now_theta,
                    q_ps: now_q,
                } = &mut current;
                let LayeredState {
                    ps: next_ps,
                    u: next_u,
                    theta_ps: next_theta,
                    q_ps: next_q,
                } = &mut *next;
                let (now_ps, next_ps) = (&*now_ps, &*next_ps);
                now_u
                    .par_iter_mut()
                    .zip(now_theta.par_iter_mut())
                    .zip(now_q.par_iter_mut())
                    .zip(next_u.par_iter_mut())
                    .zip(next_theta.par_iter_mut())
                    .zip(next_q.par_iter_mut())
                    .zip(bounds.par_iter_mut())
                    .zip(water[..n].par_iter_mut())
                    .enumerate()
                    .for_each(
                        |(
                            k,
                            (
                                ((((((now_u, now_theta), now_q), next_u), next_theta), next_q), b),
                                water,
                            ),
                        )| {
                            let old_q = &previous.q_ps[k];
                            {
                                let (now_q, next_q) = (&*now_q, &*next_q);
                                for_each_mut(b, |i, b| {
                                    let q = [
                                        old_q[i] / previous.ps[i],
                                        now_q[i] / now_ps[i],
                                        next_q[i] / next_ps[i],
                                    ];
                                    *b = (q[0].min(q[1]).min(q[2]), q[0].max(q[1]).max(q[2]));
                                });
                            }
                            water[0] = model.layer_water(k, now_q);
                            water[1] = model.layer_water(k, next_q);
                            raw_filter_part(&previous.u[k], now_u, next_u, nu, alpha);
                            raw_filter_part(
                                &previous.theta_ps[k],
                                now_theta,
                                next_theta,
                                nu,
                                alpha,
                            );
                            raw_filter_part(old_q, now_q, next_q, nu, alpha);
                            let (b, old_ps) = (&*b, &previous.ps);
                            for_each_mut2(now_q, next_q, |i, now, next| {
                                let (p_now, p_next) = raw(old_ps[i], now_ps[i], next_ps[i]);
                                *now = now.clamp(p_now * b[i].0, p_now * b[i].1);
                                *next = next.clamp(p_next * b[i].0, p_next * b[i].1);
                            });
                            water[2] = model.layer_water(k, now_q);
                            water[3] = model.layer_water(k, next_q);
                        },
                    );
            }
            raw_filter_part(&self.previous.ps, &mut current.ps, &mut next.ps, nu, alpha);
            let total = |level: usize| water[..n].iter().map(|w| w[level]).sum::<f64>();
            let before = total(0) + total(1);
            let after = total(2) + total(3);
            self.water.filter += 0.5 * (after - before);
            filtered = Some(clock.lap());
            let g = self.model.gravity();
            let ds = self.model.dsigma();
            for level in [&mut current, &mut *next] {
                let (added, borrowed) = fill_negative(&ops, ds, g, level);
                self.water.fixer += 0.5 * added;
                self.water.borrowed += 0.5 * borrowed;
            }
        } else {
            raw_filter(&self.previous, &mut current, next, nu, alpha);
        }
        let old_previous = std::mem::replace(&mut self.previous, current);
        self.current = std::mem::replace(&mut self.next, old_previous);
        let last = clock.lap();
        let fixers = match filtered {
            Some(filter) => {
                self.timings.add(Phase::Filter, filter);
                self.timings.add(Phase::NegativeHumidityFill, last);
                filter + last
            }
            None => {
                self.timings.add(Phase::Filter, last);
                last
            }
        };
        self.timings.add(Phase::TimeFilterAndFixers, fixers);
        self.model.advance_forcing(&self.current, dt);
        if forced {
            self.timings.add(Phase::PhysicsAdvance, clock.lap());
        }

        self.steps += 1;
        for r in &reports {
            self.stats.solves += 1;
            self.stats.iterations += u64::from(r.iterations);
            self.stats.worst_residual = self.stats.worst_residual.max(r.relative_residual);
        }
        reports
    }
}

// Restart array names.
const SIGMA: &str = "dynamics.sigma";
const CURRENT: &str = "dynamics.";
const PREVIOUS: &str = "dynamics.previous.";
const MODAL: &str = "dynamics.modal.";
const ENERGY_BUDGET: &str = "dynamics.energy_budget";
const WATER_BUDGET: &str = "dynamics.water_budget";
const SOLVER_STATS: &str = "dynamics.solver_stats";
const DT: &str = "dynamics.dt";

/// A level's arrays and their restart names: `ps`, then `u.{k}`,
/// `theta_ps.{k}` and (moist) `q_ps.{k}` per layer, top first.
fn level_arrays<'s>(prefix: &str, s: &'s LayeredState) -> Vec<(String, &'s [f64])> {
    let mut arrays = vec![(format!("{prefix}ps"), s.ps.as_slice())];
    for (field, layers) in [("u", &s.u), ("theta_ps", &s.theta_ps), ("q_ps", &s.q_ps)] {
        for (k, values) in layers.iter().enumerate() {
            arrays.push((format!("{prefix}{field}.{k}"), values.as_slice()));
        }
    }
    arrays
}

/// [`level_arrays`]' arrays, mutably and in the same order.
fn level_arrays_mut(s: &mut LayeredState) -> Vec<&mut [f64]> {
    let mut arrays = vec![s.ps.as_mut_slice()];
    for layers in [&mut s.u, &mut s.theta_ps, &mut s.q_ps] {
        arrays.extend(layers.iter_mut().map(Vec::as_mut_slice));
    }
    arrays
}

/// `Σᵢ cᵢ xᵢ` per cell, left to right.
fn combine(c: &[f64], x: &[Vec<f64>]) -> Vec<f64> {
    let mut out = vec![0.0; x[0].len()];
    combine_into(c, x, &mut out);
    out
}

/// [`combine`] into `out`.
fn combine_into(c: &[f64], x: &[Vec<f64>], out: &mut [f64]) {
    for_each_mut(out, |i, o| *o = combine_at(c, x, i));
}

/// `Σᵢ cᵢ xᵢ` at cell `i`, left to right.
#[inline]
fn combine_at<X: AsRef<[f64]>>(c: &[f64], x: &[X], i: usize) -> f64 {
    let mut s = c[0] * x[0].as_ref()[i];
    for (c, x) in c.iter().zip(x).skip(1) {
        s += c * x.as_ref()[i];
    }
    s
}

/// The rows of `y = (p_s, θ₀ᴸ, …)`: `ps`, then `theta_l`'s (the first
/// `theta_l.len() + 1` are set).
fn rows<'r>(ps: &'r [f64], theta_l: &'r [Vec<f64>]) -> [&'r [f64]; MAX_LAYERS + 1] {
    let mut rows: [&[f64]; MAX_LAYERS + 1] = [&[]; MAX_LAYERS + 1];
    rows[0] = ps;
    for (row, t) in rows[1..].iter_mut().zip(theta_l) {
        *row = t;
    }
    rows
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
        raw_filter_part(old, now, next, nu, alpha);
    }
}

/// The RAW filter on one field.
fn raw_filter_part(old: &[f64], now: &mut [f64], next: &mut [f64], nu: f64, alpha: f64) {
    for_each_mut2(now, next, |i, now, next| {
        let d = 0.5 * nu * (old[i] - 2.0 * *now + *next);
        *now += alpha * d;
        *next -= (1.0 - alpha) * d;
    });
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
