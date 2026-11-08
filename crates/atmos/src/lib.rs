//! Two-layer atmospheric general circulation model.
//!
//! So far: the single-layer shallow-water model ([`shallow_water`]) the
//! dynamical core is built from, stepped explicitly (RK4) or semi-implicitly
//! ([`semi_implicit`], with the [`helmholtz`] solver), and the Williamson test
//! cases that check it ([`williamson`]); and the two-layer dry primitive-
//! equation core ([`two_layer`]).
//!
//! See `docs/PLAN.md` for the design.

pub mod baroclinic;
pub mod helmholtz;
pub mod semi_implicit;
pub mod shallow_water;
pub mod two_layer;
pub mod williamson;

pub use baroclinic::ShearedFlow;
pub use helmholtz::{Helmholtz, SolveReport};
pub use semi_implicit::{SemiImplicit, SemiImplicitOptions, SolverStats};
pub use shallow_water::{Diagnostics, ShallowWater, State};
pub use two_layer::{Gas, TwoLayer, TwoLayerDiagnostics, TwoLayerState};
