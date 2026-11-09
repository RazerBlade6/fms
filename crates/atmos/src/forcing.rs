//! External forcing of the dynamical core: heating, cooling and friction from
//! outside the adiabatic dynamics. Held–Suarez ([`crate::held_suarez`]) is
//! the first; radiation and surface fluxes replace it later.

use crate::layered::LayeredState;

/// A source of tendencies for the winds, `Θ` and humidity.
///
/// Both steppers treat forcing explicitly: RK4 at every stage, the
/// semi-implicit scheme lagged to the old time level. So rates must stay
/// well below `1/Δt`.
pub trait Forcing: Send + Sync {
    /// Add the forcing tendencies of `state` to `out` (winds, `Θ` and, for a
    /// moist state, `Q`; leave `out.ps` alone, so mass is conserved). The
    /// semi-implicit stepper calls this once per step, so a forcing may
    /// record what it applied (e.g. precipitation) for
    /// [`advance`](Self::advance) to account.
    fn add_tendency(&mut self, state: &LayeredState, out: &mut LayeredState);

    /// Called once after every completed time step of `dt` seconds with the
    /// new state: for forcing with state of its own (a surface) or work done
    /// every few steps (radiation). Default: nothing.
    fn advance(&mut self, _state: &LayeredState, _dt: f64) {}

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
    /// Water evaporated into the atmosphere (kg).
    pub evaporated: f64,
    /// Water precipitated out of the atmosphere (kg).
    pub precipitated: f64,
    /// Kinetic energy the forcing removed from the winds and returned as
    /// heat (J; [`crate::dissipation`]).
    pub dissipated: f64,
    /// Water held in the soil now (kg).
    pub soil_water: f64,
    /// Water that has run off the land (kg).
    pub runoff: f64,
    /// Evaporation minus precipitation over the ocean (kg).
    pub ocean_exchange: f64,
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
