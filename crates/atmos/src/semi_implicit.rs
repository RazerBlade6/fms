//! Semi-implicit leapfrog time stepping for the shallow-water model, with the
//! Robert–Asselin–Williams (RAW) filter.
//!
//! The tendency is split into the linear gravity-wave terms about a reference
//! depth `H̄`, `L = (−H̄ ∇·u, −g ∇h)`, and the rest, `N = T − L` (advection,
//! the potential-vorticity flux, and the gravity-wave terms' departure from
//! `H̄`). `N` is evaluated at the middle time level; `L` is averaged over the
//! old and new levels (Robert 1969):
//!
//! ```text
//! xⁿ⁺¹ = xⁿ⁻¹ + 2Δt N(xⁿ) + Δt L(xⁿ⁺¹ + xⁿ⁻¹)
//! ```
//!
//! Eliminating `uⁿ⁺¹` leaves a Helmholtz equation for the thickness,
//! `(I − Δt² g H̄ ∇²) hⁿ⁺¹ = rhs`, solved by preconditioned conjugate gradient
//! ([`Helmholtz`]). Gravity waves no longer limit the step, only advection, so
//! the step can be several times longer than the explicit one. `H̄` should be
//! at least the deepest fluid, or the remaining explicit gravity waves can go
//! unstable. What still limits the step is advection (Courant number below
//! about 0.8) and the explicit Coriolis term (`f Δt < 1`); see
//! [`SemiImplicit::stable_dt`].
//!
//! The Laplacian has no net source, so the exact Helmholtz solution has the
//! same area integral as its right-hand side. After each solve the thickness
//! is shifted by the uniform constant that restores it. That keeps mass
//! conserved to round-off whatever the solver tolerance, and leaves the
//! velocity untouched (a uniform shift has no gradient).
//!
//! After each step the RAW filter (Williams 2009) damps leapfrog's
//! computational mode: with `d = ν/2 (xⁿ⁻¹ − 2xⁿ + xⁿ⁺¹)`,
//! `xⁿ += α d` and `xⁿ⁺¹ −= (1 − α) d`. With `α = ½` it changes no mean and
//! keeps third-order accuracy for the physical mode; `α = 1` is the classic
//! Robert–Asselin filter.
//!
//! Leapfrog needs two starting levels: the first step is taken with RK4
//! substeps short enough to be explicitly stable.

use fms_grid::operators::MIN_CHUNK;
use rayon::prelude::*;

use crate::helmholtz::{Helmholtz, SolveReport, dot};
use crate::shallow_water::{ShallowWater, State};

/// Settings for [`SemiImplicit`].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SemiImplicitOptions {
    /// Reference depth `H̄` for the implicit gravity-wave terms (m). Use at
    /// least the deepest fluid.
    pub reference_depth: f64,
    /// RAW filter strength `ν` (0 disables it). Typical: 0.01 to 0.2.
    pub filter_strength: f64,
    /// RAW filter weight `α` (½ to 1; 0.53 is Williams' recommendation).
    pub filter_alpha: f64,
    /// Helmholtz solver tolerance (relative residual).
    pub tolerance: f64,
}

impl SemiImplicitOptions {
    /// Defaults for a state: `H̄` the deepest fluid, `ν = 0.05`, `α = 0.53`,
    /// tolerance `1e-8`. Tighter tolerances cost iterations without changing
    /// the Williamson errors measurably (mass is exact either way).
    pub fn for_state(state: &State) -> Self {
        Self {
            reference_depth: state.h.iter().fold(0.0f64, |a, &h| a.max(h)),
            filter_strength: 0.05,
            filter_alpha: 0.53,
            tolerance: 1e-8,
        }
    }
}

/// Totals for the solves so far.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct SolverStats {
    /// Number of Helmholtz solves.
    pub solves: u64,
    /// Total conjugate-gradient iterations.
    pub iterations: u64,
    /// Worst final relative residual.
    pub worst_residual: f64,
}

/// A shallow-water model stepped semi-implicitly.
pub struct SemiImplicit<'a> {
    model: ShallowWater<'a>,
    helmholtz: Helmholtz<'a>,
    options: SemiImplicitOptions,
    dt: f64,
    /// `xⁿ⁻¹` and `xⁿ` (filtered).
    previous: State,
    current: State,
    steps: u64,
    stats: SolverStats,
    /// Σ Aᵢ on the unit sphere.
    total_area: f64,
    // Scratch.
    tendency: State,
    a: State,
    rhs: Vec<f64>,
    cells: Vec<f64>,
    edges: Vec<f64>,
}

impl<'a> SemiImplicit<'a> {
    /// Start from `state` with step `dt` (s). Takes the first step at once,
    /// with RK4 substeps at Courant number `start_courant` (e.g. 0.5).
    pub fn new(
        mut model: ShallowWater<'a>,
        state: State,
        dt: f64,
        options: SemiImplicitOptions,
        start_courant: f64,
    ) -> Self {
        let ops = *model.operators();
        let c = dt * dt * model.gravity() * options.reference_depth;
        let mut helmholtz = Helmholtz::new(ops, c);
        helmholtz.tolerance = options.tolerance;

        let substeps = (dt / model.stable_dt(&state, start_courant))
            .ceil()
            .max(1.0);
        let mut current = state.clone();
        for _ in 0..substeps as u64 {
            model.step(&mut current, dt / substeps);
        }
        let (nc, ne) = (state.h.len(), state.u.len());
        let zeros = || State {
            h: vec![0.0; nc],
            u: vec![0.0; ne],
        };
        Self {
            model,
            helmholtz,
            options,
            dt,
            previous: state,
            current,
            steps: 1,
            stats: SolverStats::default(),
            total_area: ops.geometry().cell_areas.iter().sum(),
            tendency: zeros(),
            a: zeros(),
            rhs: vec![0.0; nc],
            cells: vec![0.0; nc],
            edges: vec![0.0; ne],
        }
    }

    /// The longest stable step (s) for winds up to `max_wind` (m s⁻¹) at
    /// Courant number `courant` (0.6 leaves a margin; the Williamson cases go
    /// unstable near 0.8): the smaller of the advective limit
    /// `courant · min dₑ / max_wind` and the inertial limit `courant / max |f|`.
    /// Gravity waves don't limit it.
    ///
    /// Pass the fastest wind expected over the run, not the current one: a
    /// flow that accelerates (e.g. Williamson case 5, from 20 to 40 m s⁻¹)
    /// would otherwise outrun the step.
    pub fn stable_dt(model: &ShallowWater, max_wind: f64, courant: f64) -> f64 {
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

    /// The current state.
    pub fn state(&self) -> &State {
        &self.current
    }

    /// The underlying model (for diagnostics).
    pub fn model(&self) -> &ShallowWater<'a> {
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

    /// Helmholtz solver totals.
    pub fn solver_stats(&self) -> SolverStats {
        self.stats
    }

    /// Advance one step.
    pub fn step(&mut self) -> SolveReport {
        let ops = *self.model.operators();
        let (dt, gravity, depth) = (self.dt, self.model.gravity(), self.options.reference_depth);
        let Self {
            model,
            helmholtz,
            previous,
            current,
            tendency,
            a,
            rhs,
            cells,
            edges,
            ..
        } = self;

        // N(xⁿ) = T(xⁿ) − L(xⁿ).
        model.tendency(current, tendency);
        ops.divergence(&current.u, cells);
        ops.gradient(&current.h, edges);
        add_scaled(&mut tendency.h, depth, cells);
        add_scaled(&mut tendency.u, gravity, edges);

        // Explicit parts: Aᵤ = uⁿ⁻¹ + 2Δt Nᵤ − Δt g ∇hⁿ⁻¹,
        //                 Aₕ = hⁿ⁻¹ + 2Δt Nₕ − Δt H̄ ∇·uⁿ⁻¹.
        ops.gradient(&previous.h, edges);
        ops.divergence(&previous.u, cells);
        a.u.par_iter_mut()
            .with_min_len(MIN_CHUNK)
            .enumerate()
            .for_each(|(e, a)| {
                *a = previous.u[e] + 2.0 * dt * tendency.u[e] - dt * gravity * edges[e];
            });
        a.h.par_iter_mut()
            .with_min_len(MIN_CHUNK)
            .enumerate()
            .for_each(|(i, a)| {
                *a = previous.h[i] + 2.0 * dt * tendency.h[i] - dt * depth * cells[i];
            });

        // (I − Δt² g H̄ ∇²) hⁿ⁺¹ = Aₕ − Δt H̄ ∇·Aᵤ.
        ops.divergence(&a.u, cells);
        rhs.par_iter_mut()
            .with_min_len(MIN_CHUNK)
            .enumerate()
            .for_each(|(i, r)| *r = a.h[i] - dt * depth * cells[i]);
        // Guess hⁿ⁺¹ ≈ 2hⁿ − hⁿ⁻¹; store the new state in `tendency`.
        tendency
            .h
            .par_iter_mut()
            .with_min_len(MIN_CHUNK)
            .enumerate()
            .for_each(|(i, h)| *h = 2.0 * current.h[i] - previous.h[i]);
        let report = helmholtz.solve(rhs, &mut tendency.h);
        // Restore the exact area integral (mass) lost to the solver tolerance.
        let areas = &ops.geometry().cell_areas;
        let shift = (dot(areas, rhs) - dot(areas, &tendency.h)) / self.total_area;
        tendency
            .h
            .par_iter_mut()
            .with_min_len(MIN_CHUNK)
            .for_each(|h| *h += shift);

        // uⁿ⁺¹ = Aᵤ − Δt g ∇hⁿ⁺¹.
        ops.gradient(&tendency.h, edges);
        tendency
            .u
            .par_iter_mut()
            .with_min_len(MIN_CHUNK)
            .enumerate()
            .for_each(|(e, u)| *u = a.u[e] - dt * gravity * edges[e]);

        // RAW filter on xⁿ and xⁿ⁺¹, then shift the levels.
        let (nu, alpha) = (self.options.filter_strength, self.options.filter_alpha);
        let filter = |old: &[f64], now: &mut [f64], next: &mut [f64]| {
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
        filter(&previous.h, &mut current.h, &mut tendency.h);
        filter(&previous.u, &mut current.u, &mut tendency.u);
        std::mem::swap(previous, current);
        std::mem::swap(current, tendency);

        self.steps += 1;
        self.stats.solves += 1;
        self.stats.iterations += u64::from(report.iterations);
        self.stats.worst_residual = self.stats.worst_residual.max(report.relative_residual);
        report
    }
}

fn add_scaled(x: &mut [f64], a: f64, y: &[f64]) {
    x.par_iter_mut()
        .with_min_len(MIN_CHUNK)
        .zip(y)
        .for_each(|(x, &y)| *x += a * y);
}
