//! External forcing of the dynamical core: heating, cooling and friction from
//! outside the adiabatic dynamics. Held–Suarez ([`crate::held_suarez`]) is
//! the first; radiation and surface fluxes replace it later.

use crate::two_layer::TwoLayerState;

/// A source of tendencies for the winds and `Θ`.
///
/// Both steppers treat forcing explicitly: RK4 at every stage, the
/// semi-implicit scheme lagged to the old time level. So rates must stay
/// well below `1/Δt`.
pub trait Forcing: Send + Sync {
    /// Add the forcing tendencies of `state` to `out` (winds and `Θ`; leave
    /// `out.ps` alone, so mass is conserved).
    fn add_tendency(&self, state: &TwoLayerState, out: &mut TwoLayerState);

    /// Called once after every completed time step of `dt` seconds with the
    /// new state: for forcing with state of its own (a surface) or work done
    /// every few steps (radiation). Default: nothing.
    fn advance(&mut self, _state: &TwoLayerState, _dt: f64) {}

    /// Energy exchanged with space and stored outside the atmosphere, for
    /// the planet's energy budget. `None` if the forcing has no such notion.
    fn budget(&self) -> Option<ForcingBudget> {
        None
    }

    /// Cell fields worth showing (surface temperature, outgoing longwave …).
    fn cell_fields(&self) -> Vec<CellField> {
        Vec::new()
    }
}

/// Energy a forcing exchanges outside the atmosphere (J, summed over the
/// planet since the start, except `stored`).
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ForcingBudget {
    /// Sunlight absorbed by the planet.
    pub absorbed: f64,
    /// Longwave emitted to space.
    pub outgoing: f64,
    /// Heat held outside the atmosphere now (e.g. a surface slab), relative
    /// to a fixed zero.
    pub stored: f64,
}

/// A named cell field.
#[derive(Clone, Debug, PartialEq)]
pub struct CellField {
    /// `namespace.field` name.
    pub name: &'static str,
    /// Description for legends.
    pub description: &'static str,
    /// Units.
    pub units: &'static str,
    /// One value per cell.
    pub values: Vec<f64>,
}

/// An energy rate (W) or amount (J) split by where it acts.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct EnergyTerms {
    /// Through the winds.
    pub kinetic: f64,
    /// Through `Θ` (enthalpy).
    pub thermal: f64,
}

impl EnergyTerms {
    /// Kinetic plus thermal.
    pub fn total(&self) -> f64 {
        self.kinetic + self.thermal
    }

    /// `self + a · other`.
    pub fn add_scaled(&mut self, a: f64, other: Self) {
        self.kinetic += a * other.kinetic;
        self.thermal += a * other.thermal;
    }
}
