//! Wall-clock time spent in each phase of a run's time loop, for the timing
//! report. The steppers and forcings add to a [`Timings`] as they go: an
//! array index per phase, cheap beside a step. Phases nest ([`Phase::level`]);
//! consecutive phases are timed with one [`Stopwatch`]'s laps, so they
//! partition their parent's time, except those that run side by side
//! ([`Phase::side_by_side`]), whose times overlap.
//!
//! Timings are not saved in restart files: they cover the process's own work.

use std::time::{Duration, Instant};

/// A timed phase of the time loop, in report order: each after its parent.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
#[repr(u8)]
pub enum Phase {
    /// Copying the state before a step, for the blow-up dump.
    KeepLastGood,
    /// One semi-implicit step.
    Step,
    /// The step's tendencies, computed side by side.
    Tendencies,
    /// The forcing's tendency and its energy rate.
    Physics,
    /// Column physics: convection, condensation and clouds.
    Moisture,
    /// Column physics: the heating, moistening and drag the moisture,
    /// radiation and surface fluxes imply.
    ColumnTendencies,
    /// The adiabatic dynamics' tendency.
    Dynamics,
    /// The hyperdiffusion's tendency and its energy rate.
    Hyperdiffusion,
    /// Energy and water budgets, and the sum of the lagged tendencies.
    BudgetsAndLagged,
    /// The Helmholtz problems' right-hand side.
    ImplicitRhs,
    /// The Helmholtz solves of the vertical modes.
    HelmholtzSolves,
    /// Humidity's interface ratios, beside the solves.
    InterfaceRatios,
    /// The new winds, surface pressure and `Θ` from the solution.
    BackSubstitution,
    /// Humidity's flux-corrected transport and the fluxes it needs.
    HumidityTransport,
    /// The step's mass fluxes and their interface fluxes.
    MassFluxes,
    /// The flux-corrected transport itself.
    FluxCorrectedTransport,
    /// The RAW filter and the humidity fixers.
    TimeFilterAndFixers,
    /// The filter, humidity's bounds and the water it changed.
    Filter,
    /// Filling negative humidity.
    NegativeHumidityFill,
    /// The forcing's work after a step ([`crate::Forcing::advance`]).
    PhysicsAdvance,
    /// The ocean's surface heat, sea ice and wind stress (or a prescribed
    /// sea's uptake).
    OceanSurface,
    /// The ocean's transport (every few steps).
    OceanTransport,
    /// Land, snow and soil, and the planet's budgets.
    LandAndBudgets,
    /// A spin-up year's end and its ocean jump.
    SpinUpJump,
    /// Sea ice rearranged before radiation.
    SeaIceRearrangement,
    /// The surface fluxes for the coming step.
    SurfaceFluxes,
    /// A radiation call (every few steps).
    Radiation,
    /// The blow-up check.
    HealthCheck,
    /// A diagnostics sample.
    Sampling,
    /// The global diagnostics.
    Diagnostics,
    /// The run's own sampling (monthly means, climatology …).
    RunsOwn,
    /// GIF frames.
    GifFrames,
}

/// Number of phases.
pub const PHASES: usize = Phase::ALL.len();

/// Name, nesting level and side-by-side flag of each phase.
const INFO: [(&str, u8, bool); PHASES] = [
    ("keep the last good state", 0, false),
    ("step", 0, false),
    ("tendencies", 1, false),
    ("physics", 2, true),
    ("moisture", 3, false),
    ("column tendencies", 3, false),
    ("dynamics", 2, true),
    ("hyperdiffusion", 2, true),
    ("budgets and lagged tendencies", 1, false),
    ("implicit right-hand side", 1, false),
    ("Helmholtz solves", 1, false),
    ("interface ratios", 2, true),
    ("back-substitution", 1, false),
    ("humidity transport", 1, false),
    ("mass fluxes", 2, false),
    ("flux-corrected transport", 2, false),
    ("time filter and fixers", 1, false),
    ("filter", 2, false),
    ("negative-humidity fill", 2, false),
    ("physics advance", 1, false),
    ("ocean surface heat and sea ice", 2, false),
    ("ocean transport", 2, false),
    ("land and budgets", 2, false),
    ("spin-up jump", 2, false),
    ("sea-ice rearrangement", 2, false),
    ("surface fluxes", 2, false),
    ("radiation", 2, false),
    ("health check", 0, false),
    ("sampling", 0, false),
    ("diagnostics", 1, false),
    ("the run's own", 1, false),
    ("GIF frames", 0, false),
];

impl Phase {
    /// Every phase, in report order.
    pub const ALL: [Phase; 32] = [
        Phase::KeepLastGood,
        Phase::Step,
        Phase::Tendencies,
        Phase::Physics,
        Phase::Moisture,
        Phase::ColumnTendencies,
        Phase::Dynamics,
        Phase::Hyperdiffusion,
        Phase::BudgetsAndLagged,
        Phase::ImplicitRhs,
        Phase::HelmholtzSolves,
        Phase::InterfaceRatios,
        Phase::BackSubstitution,
        Phase::HumidityTransport,
        Phase::MassFluxes,
        Phase::FluxCorrectedTransport,
        Phase::TimeFilterAndFixers,
        Phase::Filter,
        Phase::NegativeHumidityFill,
        Phase::PhysicsAdvance,
        Phase::OceanSurface,
        Phase::OceanTransport,
        Phase::LandAndBudgets,
        Phase::SpinUpJump,
        Phase::SeaIceRearrangement,
        Phase::SurfaceFluxes,
        Phase::Radiation,
        Phase::HealthCheck,
        Phase::Sampling,
        Phase::Diagnostics,
        Phase::RunsOwn,
        Phase::GifFrames,
    ];

    /// Name for reports: short, plain, and without commas (it goes in CSV).
    pub fn name(self) -> &'static str {
        INFO[self as usize].0
    }

    /// Nesting depth: 0 for the time loop's own phases.
    pub fn level(self) -> u8 {
        INFO[self as usize].1
    }

    /// Whether it runs side by side with its siblings (inside a
    /// `rayon::join`): its time overlaps theirs, so siblings don't add up
    /// to the parent.
    pub fn side_by_side(self) -> bool {
        INFO[self as usize].2
    }
}

/// Total time and calls of each phase.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Timings {
    time: [Duration; PHASES],
    calls: [u64; PHASES],
}

impl Default for Timings {
    fn default() -> Self {
        Self {
            time: [Duration::ZERO; PHASES],
            calls: [0; PHASES],
        }
    }
}

impl Timings {
    /// Count one call of `phase` that took `duration`.
    pub fn add(&mut self, phase: Phase, duration: Duration) {
        self.time[phase as usize] += duration;
        self.calls[phase as usize] += 1;
    }

    /// Add `other`'s times and calls to these.
    pub fn merge(&mut self, other: &Self) {
        for p in 0..PHASES {
            self.time[p] += other.time[p];
            self.calls[p] += other.calls[p];
        }
    }

    /// What was added since `earlier`, a copy of these taken before.
    pub fn since(&self, earlier: &Self) -> Self {
        let mut out = *self;
        for p in 0..PHASES {
            out.time[p] = out.time[p].saturating_sub(earlier.time[p]);
            out.calls[p] = out.calls[p].saturating_sub(earlier.calls[p]);
        }
        out
    }

    /// Total time spent in `phase`.
    pub fn time(&self, phase: Phase) -> Duration {
        self.time[phase as usize]
    }

    /// Times `phase` ran.
    pub fn calls(&self, phase: Phase) -> u64 {
        self.calls[phase as usize]
    }
}

/// Times consecutive phases: each [`lap`](Self::lap) is the time since the
/// last, so the laps partition the time since the start.
#[derive(Clone, Copy, Debug)]
pub struct Stopwatch {
    last: Instant,
}

impl Stopwatch {
    /// Start timing now.
    pub fn start() -> Self {
        Self {
            last: Instant::now(),
        }
    }

    /// The time since the last lap (or the start), restarting from now.
    pub fn lap(&mut self) -> Duration {
        let now = Instant::now();
        let lap = now - self.last;
        self.last = now;
        lap
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const MS: Duration = Duration::from_millis(1);

    #[test]
    fn the_table_matches_the_phases() {
        for (i, p) in Phase::ALL.iter().enumerate() {
            assert_eq!(*p as usize, i, "{p:?} is out of order");
            assert!(!p.name().contains(','), "{p:?}'s name has a comma");
        }
        // Every phase below the top sits under a parent one level up.
        let mut depth = 0;
        for p in Phase::ALL {
            assert!(p.level() <= depth + 1, "{p:?} skips a level");
            depth = p.level();
        }
    }

    #[test]
    fn add_counts_calls_and_time() {
        let mut t = Timings::default();
        t.add(Phase::Radiation, 3 * MS);
        t.add(Phase::Radiation, 2 * MS);
        t.add(Phase::Step, MS);
        assert_eq!(t.time(Phase::Radiation), 5 * MS);
        assert_eq!(t.calls(Phase::Radiation), 2);
        assert_eq!(t.time(Phase::Step), MS);
        assert_eq!(t.calls(Phase::Step), 1);
        assert_eq!(t.calls(Phase::Dynamics), 0);
        assert_eq!(t.time(Phase::Dynamics), Duration::ZERO);
    }

    #[test]
    fn merge_sums_and_since_subtracts() {
        let mut a = Timings::default();
        a.add(Phase::Step, 4 * MS);
        a.add(Phase::Moisture, MS);
        let mut b = Timings::default();
        b.add(Phase::Step, 2 * MS);
        b.add(Phase::Radiation, 7 * MS);
        let earlier = a;
        a.merge(&b);
        assert_eq!(a.time(Phase::Step), 6 * MS);
        assert_eq!(a.calls(Phase::Step), 2);
        assert_eq!(a.time(Phase::Radiation), 7 * MS);
        assert_eq!(a.calls(Phase::Moisture), 1);
        assert_eq!(a.since(&earlier), b);
        assert_eq!(a.since(&a), Timings::default());
    }

    #[test]
    fn laps_partition_the_elapsed_time() {
        let outer = Instant::now();
        let mut s = Stopwatch::start();
        let mut total = Duration::ZERO;
        for _ in 0..3 {
            std::thread::sleep(MS);
            let lap = s.lap();
            assert!(lap >= MS);
            total += lap;
        }
        assert!(total >= 3 * MS);
        assert!(total <= outer.elapsed());
    }
}
