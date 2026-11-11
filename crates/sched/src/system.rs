//! The `System` trait and what a system declares about itself.

use fms_fields::{FieldId, FieldKey, FieldMut, FieldRef, FieldRegistry, FieldType};

use crate::clock::{Calendar, Date};

/// How often a system runs.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Rate {
    /// Every `n` base steps (`n ≥ 1`), on steps that are multiples of `n`.
    Every(u32),
    /// On the first step of each solar day.
    Daily,
    /// On the first step of each month.
    Monthly,
    /// On the first step of each year.
    Yearly,
    /// Only when requested with [`Scheduler::request`](crate::Scheduler::request),
    /// on the next step (e.g. after a terrain edit).
    OnDemand,
}

/// Which speed tier(s) a system runs in (see `docs/climate/PLAN.md` §4.6).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Tier {
    /// Speeds 1–2: the full GCM.
    Dynamics,
    /// Speeds 3–5: climatology-driven stochastic weather.
    ClimateMode,
}

/// Set of tiers a system is active in.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Tiers {
    /// Active in [`Tier::Dynamics`].
    pub dynamics: bool,
    /// Active in [`Tier::ClimateMode`].
    pub climate_mode: bool,
}

impl Tiers {
    /// Active in every tier.
    pub const ALL: Self = Self {
        dynamics: true,
        climate_mode: true,
    };
    /// Active only with the full GCM.
    pub const DYNAMICS: Self = Self {
        dynamics: true,
        climate_mode: false,
    };
    /// Active only in climate mode.
    pub const CLIMATE_MODE: Self = Self {
        dynamics: false,
        climate_mode: true,
    };

    /// Whether `tier` is in the set.
    pub fn contains(self, tier: Tier) -> bool {
        match tier {
            Tier::Dynamics => self.dynamics,
            Tier::ClimateMode => self.climate_mode,
        }
    }
}

/// Fields a system reads and writes. Writing implies reading.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Access {
    reads: Vec<FieldKey>,
    writes: Vec<FieldKey>,
}

impl Access {
    /// No access.
    pub fn new() -> Self {
        Self::default()
    }

    /// Declare a read.
    pub fn read<T>(mut self, id: FieldId<T>) -> Self {
        insert(&mut self.reads, id.key());
        self
    }

    /// Declare a write (read-write).
    pub fn write<T>(mut self, id: FieldId<T>) -> Self {
        insert(&mut self.writes, id.key());
        self
    }

    /// Fields read but not written.
    pub fn reads(&self) -> impl Iterator<Item = FieldKey> + '_ {
        self.reads
            .iter()
            .copied()
            .filter(|k| !self.writes.contains(k))
    }

    /// Fields written.
    pub fn writes(&self) -> &[FieldKey] {
        &self.writes
    }

    /// Whether `key` may be read.
    pub fn can_read(&self, key: FieldKey) -> bool {
        self.reads.contains(&key) || self.writes.contains(&key)
    }

    /// Whether `key` may be written.
    pub fn can_write(&self, key: FieldKey) -> bool {
        self.writes.contains(&key)
    }

    /// Whether two systems with these accesses must not run concurrently: one
    /// writes a field the other reads or writes.
    pub fn conflicts_with(&self, other: &Access) -> bool {
        self.writes.iter().any(|&k| other.can_read(k))
            || other.writes.iter().any(|&k| self.can_read(k))
    }
}

fn insert(keys: &mut Vec<FieldKey>, key: FieldKey) {
    if let Err(pos) = keys.binary_search(&key) {
        keys.insert(pos, key);
    }
}

/// What a system declares about itself, once, when it is added.
#[derive(Clone, Debug)]
pub struct SystemDesc {
    /// Unique name, e.g. `"atmos.dynamics"`.
    pub name: &'static str,
    /// Fields read and written.
    pub access: Access,
    /// How often it runs.
    pub rate: Rate,
    /// Tiers it is active in.
    pub tiers: Tiers,
}

/// A unit of simulation work: radiation, dynamics, hydrology, a civ step, ….
///
/// Systems must touch fields only through the [`StepContext`], which checks every
/// access against the declared [`Access`]. Given that, the scheduler can run
/// non-conflicting systems in parallel and still produce results identical to
/// running them one after another in registration order.
pub trait System: Send {
    /// Declare name, access, rate and tiers. Called once, by
    /// [`Scheduler::add`](crate::Scheduler::add).
    fn describe(&self) -> SystemDesc;

    /// Advance this system's state by one of its intervals.
    fn run(&mut self, ctx: &StepContext<'_>);
}

/// What a system sees while running: checked field access and the time.
pub struct StepContext<'a> {
    pub(crate) fields: &'a FieldRegistry,
    pub(crate) access: &'a Access,
    pub(crate) name: &'static str,
    pub(crate) calendar: &'a Calendar,
    pub(crate) date: Date,
    pub(crate) interval_seconds: f64,
}

impl StepContext<'_> {
    /// Read a declared field.
    ///
    /// # Panics
    /// If the field is not in the system's declared reads or writes.
    pub fn read<T: FieldType>(&self, id: FieldId<T>) -> FieldRef<'_, T> {
        self.check(self.access.can_read(id.key()), id.key(), "read");
        self.fields.read(id)
    }

    /// Write a declared field in place.
    ///
    /// # Panics
    /// If the field is not in the system's declared writes.
    pub fn write<T: FieldType>(&self, id: FieldId<T>) -> FieldMut<'_, T> {
        self.check(self.access.can_write(id.key()), id.key(), "write");
        self.fields.write(id)
    }

    /// Read the back buffer of a declared double-buffered field.
    ///
    /// # Panics
    /// If the field is not declared, or is single buffered.
    pub fn read_back<T: FieldType>(&self, id: FieldId<T>) -> FieldRef<'_, T> {
        self.check(self.access.can_read(id.key()), id.key(), "read");
        self.fields.read_back(id)
    }

    /// Write the back buffer of a declared double-buffered field.
    ///
    /// # Panics
    /// If the field is not in the declared writes, or is single buffered.
    pub fn write_back<T: FieldType>(&self, id: FieldId<T>) -> FieldMut<'_, T> {
        self.check(self.access.can_write(id.key()), id.key(), "write");
        self.fields.write_back(id)
    }

    /// Swap the buffers of a declared double-buffered field.
    ///
    /// # Panics
    /// If the field is not in the declared writes, or is single buffered.
    pub fn swap<T: FieldType>(&self, id: FieldId<T>) {
        self.check(self.access.can_write(id.key()), id.key(), "write");
        self.fields.swap(id.key());
    }

    /// The current date (before this step's advance).
    pub fn date(&self) -> Date {
        self.date
    }

    /// The world's calendar.
    pub fn calendar(&self) -> &Calendar {
        self.calendar
    }

    /// Base step `dt` in seconds.
    pub fn dt(&self) -> f64 {
        self.calendar.dt()
    }

    /// Simulated seconds this run covers: the time since the system last ran, or
    /// its nominal interval on the first run. Use this, not `dt`, as the timestep
    /// for anything integrated in time.
    pub fn interval_seconds(&self) -> f64 {
        self.interval_seconds
    }

    fn check(&self, allowed: bool, key: fms_fields::FieldKey, what: &str) {
        if !allowed {
            panic!(
                "system `{}` tried to {what} field `{}` without declaring it",
                self.name,
                self.fields.meta(key).spec.name
            );
        }
    }
}
