//! Blow-up detection: catch a numerical instability while the state is still
//! readable, before it turns into NaNs.
//!
//! [`check`] scans a state once and either summarizes it ([`Health`]) or
//! reports the worst problem it found ([`BlowUp`]): a value that isn't finite,
//! a wind fast enough to break the time step's Courant limit, or a surface
//! pressure or potential temperature outside generous physical bounds.

use std::fmt;

use crate::two_layer::TwoLayerState;

/// Bounds a healthy state stays within.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Limits {
    /// Largest stable wind speed (m s⁻¹): the advective Courant limit of the
    /// time step.
    pub max_wind: f64,
    /// Surface pressure range (Pa).
    pub surface_pressure: (f64, f64),
    /// Potential temperature range (K).
    pub theta: (f64, f64),
}

impl Limits {
    /// Limits for a step stable up to `max_wind`, about a reference surface
    /// pressure `p̄` and layer potential temperatures `θ̄` (upper, lower):
    /// `p_s` within `[p̄/4, 2p̄]` and θ within `[θ̄₁/3, 3θ̄₀]`. These are
    /// far outside anything a planet's weather produces; crossing them means
    /// the numerics have failed.
    pub fn around(max_wind: f64, p_bar: f64, theta_bar: [f64; 2]) -> Self {
        Self {
            max_wind,
            surface_pressure: (0.25 * p_bar, 2.0 * p_bar),
            theta: (theta_bar[1] / 3.0, 3.0 * theta_bar[0]),
        }
    }
}

/// A prognostic variable.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Variable {
    /// Surface pressure (cells).
    SurfacePressure,
    /// Normal wind of a layer (edges); 0 is the upper layer.
    Wind(usize),
    /// Potential temperature of a layer (cells).
    Theta(usize),
}

impl Variable {
    /// Whether the variable lives at edges (else cells).
    pub fn at_edges(self) -> bool {
        matches!(self, Self::Wind(_))
    }
}

impl fmt::Display for Variable {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let layer = |k: &usize| if *k == 0 { "upper" } else { "lower" };
        match self {
            Self::SurfacePressure => write!(f, "surface pressure"),
            Self::Wind(k) => write!(f, "{} wind", layer(k)),
            Self::Theta(k) => write!(f, "{} θ", layer(k)),
        }
    }
}

/// What went wrong.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Fault {
    /// NaN or infinite.
    NotFinite,
    /// Wind beyond the Courant limit.
    TooFast,
    /// Outside the physical bounds.
    OutOfRange,
}

/// The worst problem in a state.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BlowUp {
    /// The variable at fault.
    pub variable: Variable,
    /// What is wrong with it.
    pub fault: Fault,
    /// Cell or edge index (see [`Variable::at_edges`]).
    pub index: usize,
    /// The offending value (m s⁻¹, Pa or K).
    pub value: f64,
    /// The limit it crossed (NaN for [`Fault::NotFinite`]).
    pub limit: f64,
}

impl fmt::Display for BlowUp {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let place = if self.variable.at_edges() {
            "edge"
        } else {
            "cell"
        };
        let (v, i) = (self.variable, self.index);
        match self.fault {
            Fault::NotFinite => write!(f, "{v} is {} at {place} {i}", self.value),
            Fault::TooFast => write!(
                f,
                "{v} reached {:.1} m/s at {place} {i}, beyond the time step's Courant limit of {:.1} m/s",
                self.value, self.limit
            ),
            Fault::OutOfRange => write!(
                f,
                "{v} is {:.6e} at {place} {i}, outside the physical bound {:.6e}",
                self.value, self.limit
            ),
        }
    }
}

impl std::error::Error for BlowUp {}

/// Summary of a healthy state.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Health {
    /// Fastest normal wind (m s⁻¹).
    pub max_wind: f64,
    /// Edge and layer of the fastest wind.
    pub max_wind_at: (usize, usize),
    /// Surface pressure range (Pa).
    pub surface_pressure: (f64, f64),
    /// Potential temperature range over both layers (K).
    pub theta: (f64, f64),
}

/// Check `state` against `limits`. Non-finite values are reported first,
/// then the fastest wind if it's too fast, then the first value outside the
/// bounds. Deterministic: a fixed scan order.
pub fn check(state: &TwoLayerState, limits: &Limits) -> Result<Health, BlowUp> {
    let not_finite = |variable, index: usize, value: f64| BlowUp {
        variable,
        fault: Fault::NotFinite,
        index,
        value,
        limit: f64::NAN,
    };
    let mut health = Health {
        max_wind: 0.0,
        max_wind_at: (0, 0),
        surface_pressure: (f64::INFINITY, f64::NEG_INFINITY),
        theta: (f64::INFINITY, f64::NEG_INFINITY),
    };
    for (i, &p) in state.ps.iter().enumerate() {
        if !p.is_finite() {
            return Err(not_finite(Variable::SurfacePressure, i, p));
        }
        health.surface_pressure = (
            health.surface_pressure.0.min(p),
            health.surface_pressure.1.max(p),
        );
    }
    for k in 0..2 {
        for (e, &u) in state.u[k].iter().enumerate() {
            if !u.is_finite() {
                return Err(not_finite(Variable::Wind(k), e, u));
            }
            if u.abs() > health.max_wind {
                health.max_wind = u.abs();
                health.max_wind_at = (e, k);
            }
        }
        for (i, (&t, &p)) in state.theta_ps[k].iter().zip(&state.ps).enumerate() {
            let theta = t / p;
            if !theta.is_finite() {
                return Err(not_finite(Variable::Theta(k), i, theta));
            }
            health.theta = (health.theta.0.min(theta), health.theta.1.max(theta));
        }
    }
    if health.max_wind > limits.max_wind {
        let (e, k) = health.max_wind_at;
        return Err(BlowUp {
            variable: Variable::Wind(k),
            fault: Fault::TooFast,
            index: e,
            value: state.u[k][e],
            limit: limits.max_wind,
        });
    }
    let out_of_range = |variable, index, value: f64, (lo, hi): (f64, f64)| {
        (value < lo || value > hi).then_some(BlowUp {
            variable,
            fault: Fault::OutOfRange,
            index,
            value,
            limit: if value < lo { lo } else { hi },
        })
    };
    for (i, &p) in state.ps.iter().enumerate() {
        if let Some(b) = out_of_range(Variable::SurfacePressure, i, p, limits.surface_pressure) {
            return Err(b);
        }
    }
    for k in 0..2 {
        for (i, (&t, &p)) in state.theta_ps[k].iter().zip(&state.ps).enumerate() {
            if let Some(b) = out_of_range(Variable::Theta(k), i, t / p, limits.theta) {
                return Err(b);
            }
        }
    }
    Ok(health)
}
