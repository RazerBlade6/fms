//! Two-layer atmospheric general circulation model.
//!
//! So far: the single-layer shallow-water model ([`shallow_water`]) the
//! dynamical core is built from, stepped explicitly (RK4) or semi-implicitly
//! ([`semi_implicit`], with the [`helmholtz`] solver), and the Williamson test
//! cases that check it ([`williamson`]); and the two-layer dry primitive-
//! equation core ([`two_layer`]), stepped explicitly or semi-implicitly
//! ([`two_layer_si`]), with scale-selective [`diffusion`], blow-up
//! detection ([`safety`]) and external [`forcing`] such as [`held_suarez`].
//!
//! See `docs/PLAN.md` for the design.

pub mod baroclinic;
pub mod climatology;
pub mod diffusion;
pub mod forcing;
pub mod held_suarez;
pub mod helmholtz;
pub mod safety;
pub mod semi_implicit;
pub mod shallow_water;
pub mod two_layer;
pub mod two_layer_si;
pub mod williamson;

pub use baroclinic::ShearedFlow;
pub use climatology::{ZonalClimatology, ZonalMeans};
pub use diffusion::{Hyperdiffusion, HyperdiffusionOptions};
pub use forcing::{EnergyTerms, Forcing};
pub use held_suarez::{HeldSuarez, HeldSuarezParams};
pub use helmholtz::{Helmholtz, SolveReport};
pub use safety::{BlowUp, Fault, Health, Limits, Variable};
pub use semi_implicit::{SemiImplicit, SemiImplicitOptions, SolverStats};
pub use shallow_water::{Diagnostics, ShallowWater, State};
pub use two_layer::{Gas, Linearization, TwoLayer, TwoLayerDiagnostics, TwoLayerState};
pub use two_layer_si::{EnergyBudget, TwoLayerSemiImplicit, TwoLayerSemiImplicitOptions};
