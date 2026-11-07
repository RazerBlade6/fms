//! System scheduler: read/write sets, tick rates, deterministic execution, and
//! the simulation clock.
//!
//! A [`System`] declares the fields it reads and writes, how often it runs
//! ([`Rate`]) and in which speed tiers ([`Tiers`]). The [`Scheduler`] advances
//! time in fixed base steps of `dt`, runs whichever systems are due, and
//! parallelizes those that don't conflict, so the result is always identical to
//! running them one after another in the order they were added.
//!
//! See `docs/PLAN.md` §3.1 and §4.6 for the design.

mod clock;
mod scheduler;
mod system;

pub use clock::{Calendar, CalendarError, Date};
pub use scheduler::{SchedError, Scheduler, StepReport, SystemId, SystemProfile};
pub use system::{Access, Rate, StepContext, System, SystemDesc, Tier, Tiers};
