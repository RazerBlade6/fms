//! Deterministic, parallel execution of systems.

use std::fmt;
use std::time::{Duration, Instant};

use fms_fields::{FieldRegistry, validate_name};
use rayon::prelude::*;

use crate::clock::{Calendar, Date};
use crate::system::{Rate, StepContext, System, SystemDesc, Tier};

/// Handle to a system added to a [`Scheduler`].
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct SystemId(u32);

impl SystemId {
    /// Registration index.
    pub fn index(self) -> usize {
        self.0 as usize
    }
}

/// Errors from adding a system.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum SchedError {
    /// A system with this name already exists.
    DuplicateName(&'static str),
    /// The name is not `namespace.system` snake_case.
    InvalidName(&'static str),
    /// `Rate::Every(0)`.
    ZeroRate(&'static str),
}

impl fmt::Display for SchedError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            SchedError::DuplicateName(n) => write!(f, "system `{n}` is already added"),
            SchedError::InvalidName(n) => {
                write!(
                    f,
                    "system name `{n}` must be `namespace.system`, both parts snake_case"
                )
            }
            SchedError::ZeroRate(n) => write!(f, "system `{n}` has rate Every(0)"),
        }
    }
}

impl std::error::Error for SchedError {}

/// What happened in one step.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct StepReport {
    /// The date the step ran at (before advancing).
    pub date: Date,
    /// Systems that ran, grouped into stages. Stages ran in order; systems within
    /// a stage ran in parallel.
    pub stages: Vec<Vec<SystemId>>,
    /// Wall-clock time each system took this step, in the order they ran.
    /// Diagnostic only: never feeds back into the simulation.
    pub timings: Vec<(SystemId, Duration)>,
    /// Wall-clock time of the whole step, including scheduling overhead.
    pub elapsed: Duration,
}

/// Accumulated wall-clock cost of one system.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SystemProfile {
    /// The system.
    pub id: SystemId,
    /// Its name.
    pub name: &'static str,
    /// Number of times it ran.
    pub runs: u64,
    /// Total wall-clock time across those runs.
    pub total: Duration,
}

struct Slot {
    system: Box<dyn System>,
    desc: SystemDesc,
    last_run: Option<u64>,
    requested: bool,
    runs: u64,
    total: Duration,
    last_elapsed: Duration,
}

/// Runs systems at their rates, in a fixed order, in parallel where safe.
///
/// **Ordering contract:** the result of a step is always identical to running
/// the due systems one after another in the order they were added. To get
/// parallelism without breaking that, each due system goes in the earliest
/// stage after every earlier-added due system it conflicts with (see
/// [`Access::conflicts_with`](crate::Access::conflicts_with)). Systems in the same
/// stage touch disjoint data, so running them concurrently commutes. The stage
/// plan does not depend on thread count, so results are identical on any number
/// of threads (and, with `fms-math`, on any platform).
pub struct Scheduler {
    calendar: Calendar,
    step: u64,
    tier: Tier,
    slots: Vec<Slot>,
}

impl Scheduler {
    /// A scheduler at step 0 in the [`Tier::Dynamics`] tier.
    pub fn new(calendar: Calendar) -> Self {
        Self {
            calendar,
            step: 0,
            tier: Tier::Dynamics,
            slots: Vec::new(),
        }
    }

    /// Add a system. Order matters: see the ordering contract on [`Scheduler`].
    pub fn add(&mut self, system: impl System + 'static) -> Result<SystemId, SchedError> {
        let desc = system.describe();
        if !validate_name(desc.name) {
            return Err(SchedError::InvalidName(desc.name));
        }
        if self.slots.iter().any(|s| s.desc.name == desc.name) {
            return Err(SchedError::DuplicateName(desc.name));
        }
        if desc.rate == Rate::Every(0) {
            return Err(SchedError::ZeroRate(desc.name));
        }
        let id = SystemId(self.slots.len() as u32);
        self.slots.push(Slot {
            system: Box::new(system),
            desc,
            last_run: None,
            requested: false,
            runs: 0,
            total: Duration::ZERO,
            last_elapsed: Duration::ZERO,
        });
        Ok(id)
    }

    /// The calendar.
    pub fn calendar(&self) -> &Calendar {
        &self.calendar
    }

    /// Steps taken so far (the step the next [`Scheduler::step`] will run).
    pub fn current_step(&self) -> u64 {
        self.step
    }

    /// Date of the next step.
    pub fn date(&self) -> Date {
        self.calendar.date(self.step)
    }

    /// Active tier.
    pub fn tier(&self) -> Tier {
        self.tier
    }

    /// Switch tier; systems not active in it are skipped until switched back.
    pub fn set_tier(&mut self, tier: Tier) {
        self.tier = tier;
    }

    /// Descriptor of a system.
    pub fn describe(&self, id: SystemId) -> &SystemDesc {
        &self.slots[id.index()].desc
    }

    /// Ask an [`Rate::OnDemand`] (or any) system to run on the next step, once.
    pub fn request(&mut self, id: SystemId) {
        self.slots[id.index()].requested = true;
    }

    /// Stage plan for the next step, without running it.
    pub fn plan(&self) -> Vec<Vec<SystemId>> {
        let date = self.date();
        let due: Vec<usize> = (0..self.slots.len())
            .filter(|&i| self.is_due(i, &date))
            .collect();
        let mut stage_of = vec![0usize; self.slots.len()];
        let mut stages: Vec<Vec<SystemId>> = Vec::new();
        for (n, &i) in due.iter().enumerate() {
            let stage = due[..n]
                .iter()
                .filter(|&&j| {
                    self.slots[j]
                        .desc
                        .access
                        .conflicts_with(&self.slots[i].desc.access)
                })
                .map(|&j| stage_of[j] + 1)
                .max()
                .unwrap_or(0);
            stage_of[i] = stage;
            if stages.len() <= stage {
                stages.resize_with(stage + 1, Vec::new);
            }
            stages[stage].push(SystemId(i as u32));
        }
        stages
    }

    fn is_due(&self, i: usize, date: &Date) -> bool {
        let slot = &self.slots[i];
        if !slot.desc.tiers.contains(self.tier) {
            return false;
        }
        slot.requested
            || match slot.desc.rate {
                Rate::Every(n) => date.step.is_multiple_of(u64::from(n)),
                Rate::Daily => date.is_day_start(),
                Rate::Monthly => date.is_month_start(&self.calendar),
                Rate::Yearly => date.is_year_start(),
                Rate::OnDemand => false,
            }
    }

    fn nominal_interval(&self, rate: Rate) -> f64 {
        let day = self.calendar.day_seconds();
        match rate {
            Rate::Every(n) => f64::from(n) * self.calendar.dt(),
            Rate::Daily => day,
            Rate::Monthly => day * self.calendar.days_per_year() / 12.0,
            Rate::Yearly => day * self.calendar.days_per_year(),
            Rate::OnDemand => 0.0,
        }
    }

    /// Run every due system for the current step, then advance the clock by one
    /// base step.
    pub fn step(&mut self, fields: &FieldRegistry) -> StepReport {
        let step_start = Instant::now();
        let date = self.date();
        let stages = self.plan();
        let dt = self.calendar.dt();

        // Per-system interval: time since its last run, or its nominal interval.
        let intervals: Vec<f64> = self
            .slots
            .iter()
            .map(|s| match s.last_run {
                Some(last) => (date.step - last) as f64 * dt,
                None => self.nominal_interval(s.desc.rate),
            })
            .collect();

        let calendar = self.calendar;
        for stage in &stages {
            let mut members: Vec<(&mut Slot, f64)> = self
                .slots
                .iter_mut()
                .zip(&intervals)
                .enumerate()
                .filter(|(i, _)| stage.contains(&SystemId(*i as u32)))
                .map(|(_, (slot, &interval))| (slot, interval))
                .collect();
            let run = |(slot, interval): &mut (&mut Slot, f64)| {
                let Slot {
                    system,
                    desc,
                    last_elapsed,
                    ..
                } = &mut **slot;
                let ctx = StepContext {
                    fields,
                    access: &desc.access,
                    name: desc.name,
                    calendar: &calendar,
                    date,
                    interval_seconds: *interval,
                };
                let start = Instant::now();
                system.run(&ctx);
                *last_elapsed = start.elapsed();
            };
            if members.len() == 1 {
                members.iter_mut().for_each(run);
            } else {
                members.par_iter_mut().for_each(run);
            }
        }

        let mut timings = Vec::new();
        for &id in stages.iter().flatten() {
            let slot = &mut self.slots[id.index()];
            slot.last_run = Some(date.step);
            slot.requested = false;
            slot.runs += 1;
            slot.total += slot.last_elapsed;
            timings.push((id, slot.last_elapsed));
        }
        self.step += 1;
        StepReport {
            date,
            stages,
            timings,
            elapsed: step_start.elapsed(),
        }
    }

    /// Accumulated cost of every system since it was added (or since
    /// [`Scheduler::reset_profile`]), in registration order.
    pub fn profile(&self) -> Vec<SystemProfile> {
        self.slots
            .iter()
            .enumerate()
            .map(|(i, s)| SystemProfile {
                id: SystemId(i as u32),
                name: s.desc.name,
                runs: s.runs,
                total: s.total,
            })
            .collect()
    }

    /// Zero the accumulated profile.
    pub fn reset_profile(&mut self) {
        for s in &mut self.slots {
            s.runs = 0;
            s.total = Duration::ZERO;
        }
    }

    /// Run `n` steps.
    pub fn run_steps(&mut self, fields: &FieldRegistry, n: u64) {
        for _ in 0..n {
            self.step(fields);
        }
    }
}
