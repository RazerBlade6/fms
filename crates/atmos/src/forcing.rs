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
