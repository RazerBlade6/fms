//! Two-layer atmospheric general circulation model.
//!
//! So far: the single-layer shallow-water model ([`shallow_water`]) the
//! dynamical core is built from, and the Williamson test cases that check it
//! ([`williamson`]).
//!
//! See `docs/PLAN.md` for the design.

pub mod shallow_water;
pub mod williamson;

pub use shallow_water::{Diagnostics, ShallowWater, State};
