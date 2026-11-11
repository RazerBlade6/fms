//! Layered atmospheric general circulation model.
//!
//! The single-layer shallow-water model ([`shallow_water`]) the dynamical
//! core is built from, stepped explicitly (RK4) or semi-implicitly
//! ([`semi_implicit`], with the [`helmholtz`] solver), and the Williamson test
//! cases that check it ([`williamson`]); the layered primitive-equation core
//! ([`layered`], any number of σ layers: [`levels`]), stepped explicitly or
//! semi-implicitly ([`layered_si`]), with scale-selective [`diffusion`],
//! blow-up detection ([`safety`]) and external [`forcing`]: [`held_suarez`],
//! or the radiation, moisture and surface of [`column_physics`].
//!
//! See `docs/DESIGN.md` for the design.

pub mod baroclinic;
pub mod climatology;
pub mod column_cases;
pub mod column_physics;
pub mod diffusion;
pub mod dissipation;
pub mod forcing;
pub mod held_suarez;
pub mod helmholtz;
pub mod layered;
pub mod layered_si;
pub mod levels;
pub mod moisture;
pub mod radiation;
pub mod safety;
pub mod semi_implicit;
pub mod shallow_water;
pub mod surface_layer;
pub mod transport;
pub mod williamson;

pub use baroclinic::ShearedFlow;
pub use climatology::{ZonalClimatology, ZonalMeans};
pub use column_physics::{
    ColumnPhysics, Equilibrium, OrographyParams, PhysicsParams, SingleColumn,
};
pub use diffusion::{Hyperdiffusion, HyperdiffusionOptions};
pub use fms_ocean::{EARTH_ROTATION, EddyParams, EkmanParams, OceanParams};
pub use forcing::{CellField, EnergyTerms, Forcing, ForcingBudget};
pub use held_suarez::{HeldSuarez, HeldSuarezParams};
pub use helmholtz::{Helmholtz, SolveReport};
pub use layered::{Gas, LayeredDiagnostics, LayeredModel, LayeredState, Linearization};
pub use layered_si::{
    EnergyBudget, LayeredSemiImplicit, LayeredSemiImplicitOptions, REFERENCE_TEMPERATURE,
    WaterBudget,
};
pub use levels::Levels;
pub use moisture::{CloudParams, MoistTendency, MoistureParams};
pub use radiation::{Cloud, Column, Fluxes, Humidity, RadiationParams};
pub use safety::{BlowUp, Fault, Health, Limits, Variable};
pub use semi_implicit::{SemiImplicit, SemiImplicitOptions, SolverStats};
pub use shallow_water::{Diagnostics, ShallowWater, State};
