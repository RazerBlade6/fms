//! Blow-up detection: catch a numerical instability while the state is still
//! readable, before it turns into NaNs.
//!
//! [`check`] scans a state once and either summarizes it ([`Health`]) or
//! reports the worst problem it found ([`BlowUp`]): a value that isn't finite,
//! a wind fast enough to break the time step's Courant limit, or a surface
//! pressure or potential temperature outside generous physical bounds.

use std::fmt;

use fms_grid::operators::{SPLIT_CHUNK, splits};
use rayon::prelude::*;

use crate::layered::LayeredState;
use crate::levels::MAX_LAYERS;

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
    /// pressure `p̄` and layer potential temperatures `θ̄` (top first):
    /// `p_s` within `[p̄/4, 2p̄]` and θ within a third of the lowest layer's
    /// and three times the top layer's. These are
    /// far outside anything a planet's weather produces; crossing them means
    /// the numerics have failed.
    pub fn around(max_wind: f64, p_bar: f64, theta_bar: &[f64]) -> Self {
        Self {
            max_wind,
            surface_pressure: (0.25 * p_bar, 2.0 * p_bar),
            theta: (theta_bar[theta_bar.len() - 1] / 3.0, 3.0 * theta_bar[0]),
        }
    }
}

/// A prognostic variable.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Variable {
    /// Surface pressure (cells).
    SurfacePressure,
    /// Normal wind of a layer (edges); 0 is the top layer.
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
        match self {
            Self::SurfacePressure => write!(f, "surface pressure"),
            Self::Wind(k) => write!(f, "layer {k} wind"),
            Self::Theta(k) => write!(f, "layer {k} θ"),
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
///
/// On a grid whose loops are split between threads, the fields are first
/// summarized side by side; only a state with something wrong is scanned
/// again in order, to find what to report.
pub fn check(state: &LayeredState, limits: &Limits) -> Result<Health, BlowUp> {
    if splits(state.ps.len())
        && let Some(health) = healthy(state, limits)
    {
        return Ok(health);
    }
    scan(state, limits)
}

/// One field's summary: whether every value is finite and within its
/// bounds (the wind's: within its limit), the range of the values (cells),
/// and the largest magnitude and the first index that has it (edges).
#[derive(Clone, Copy, Debug)]
struct Summary {
    sound: bool,
    range: (f64, f64),
    fastest: (f64, usize),
}

impl Summary {
    const EMPTY: Self = Self {
        sound: true,
        range: (f64::INFINITY, f64::NEG_INFINITY),
        fastest: (0.0, 0),
    };

    /// Both parts' summary. The range of sound values doesn't depend on
    /// the order (they hold no NaN, and bounds that exclude zero hold no
    /// zero of either sign), and the fastest is the first largest, so
    /// parts may be joined in any order.
    fn join(self, other: Self) -> Self {
        let (a, b) = (self.fastest, other.fastest);
        Self {
            sound: self.sound && other.sound,
            range: (
                self.range.0.min(other.range.0),
                self.range.1.max(other.range.1),
            ),
            fastest: if b.0 > a.0 || (b.0 == a.0 && b.1 < a.1) {
                b
            } else {
                a
            },
        }
    }
}

/// The state's [`Health`] if nothing in it is wrong: [`scan`]'s answer,
/// from summaries of the fields computed side by side. `None` if anything
/// is: not finite, too fast or out of bounds.
fn healthy(state: &LayeredState, limits: &Limits) -> Option<Health> {
    // Bounds that hold zero would let the ranges depend on the order.
    if !(limits.surface_pressure.0 > 0.0 && limits.theta.0 > 0.0) {
        return None;
    }
    let n = state.layers();
    // Surface pressure, then each layer's wind and θ.
    let mut summaries = [Summary::EMPTY; 1 + 2 * MAX_LAYERS];
    summaries[..1 + 2 * n]
        .par_iter_mut()
        .enumerate()
        .for_each(|(j, s)| {
            *s = match j {
                0 => cells(state.ps.len(), limits.surface_pressure, |i| state.ps[i]),
                j if j <= n => edges(&state.u[j - 1]),
                j => {
                    let (theta_ps, ps) = (&state.theta_ps[j - 1 - n], &state.ps);
                    cells(ps.len(), limits.theta, |i| theta_ps[i] / ps[i])
                }
            };
        });
    let (ps, layers) = summaries[..1 + 2 * n]
        .split_first()
        .expect("surface pressure");
    let (winds, thetas) = layers.split_at(n);
    // The fastest wind: the first layer's of the largest.
    let mut fastest = (0.0, (0, 0));
    for (k, w) in winds.iter().enumerate() {
        if w.fastest.0 > fastest.0 {
            fastest = (w.fastest.0, (w.fastest.1, k));
        }
    }
    let theta = thetas.iter().fold(Summary::EMPTY, |a, &b| a.join(b));
    let sound = ps.sound && theta.sound && winds.iter().all(|w| w.sound);
    (sound && fastest.0 <= limits.max_wind).then_some(Health {
        max_wind: fastest.0,
        max_wind_at: fastest.1,
        surface_pressure: ps.range,
        theta: theta.range,
    })
}

/// The [`Summary`] of the `len` values `value(i)` of a cell field, within
/// `bounds`.
fn cells(len: usize, (lo, hi): (f64, f64), value: impl Fn(usize) -> f64 + Sync) -> Summary {
    (0..len)
        .into_par_iter()
        .with_min_len(SPLIT_CHUNK)
        .fold(
            || Summary::EMPTY,
            |s, i| {
                let v = value(i);
                Summary {
                    sound: s.sound && v.is_finite() && !(v < lo || v > hi),
                    range: (s.range.0.min(v), s.range.1.max(v)),
                    ..s
                }
            },
        )
        .reduce(|| Summary::EMPTY, Summary::join)
}

/// The [`Summary`] of a wind field `u`.
fn edges(u: &[f64]) -> Summary {
    u.par_iter()
        .with_min_len(SPLIT_CHUNK)
        .enumerate()
        .fold(
            || Summary::EMPTY,
            |s, (e, &u)| Summary {
                sound: s.sound && u.is_finite(),
                fastest: if u.abs() > s.fastest.0 {
                    (u.abs(), e)
                } else {
                    s.fastest
                },
                ..s
            },
        )
        .reduce(|| Summary::EMPTY, Summary::join)
}

/// [`check`] by one scan of the state in order.
fn scan(state: &LayeredState, limits: &Limits) -> Result<Health, BlowUp> {
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
    for k in 0..state.layers() {
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
    for k in 0..state.layers() {
        for (i, (&t, &p)) in state.theta_ps[k].iter().zip(&state.ps).enumerate() {
            if let Some(b) = out_of_range(Variable::Theta(k), i, t / p, limits.theta) {
                return Err(b);
            }
        }
    }
    Ok(health)
}
