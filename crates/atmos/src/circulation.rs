//! Circulation diagnostics: the meridional overturning and the measures of
//! the general circulation read from it and from the zonal winds — the
//! overturning cells (Hadley, Ferrel, polar) and their edges, where the trade
//! winds end, and the jets.
//!
//! The analysis functions are pure functions of latitude profiles, so they
//! work on any run's time means (or on analytic profiles in tests).

use crate::climatology::ZonalMeans;

/// Overturning runs weaker than this fraction of the profile's largest |ψ|
/// are noise, not cells.
pub const CELL_FLOOR: f64 = 0.03;
/// The smallest topographic prominence of a jet (m s⁻¹).
pub const JET_PROMINENCE: f64 = 2.0;
/// The smallest prominence of a jet as a fraction of the profile's strongest
/// westerly.
pub const JET_FRACTION: f64 = 0.15;

/// Time-mean meridional overturning: the mass streamfunction on band boundaries.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Overturning {
    /// Latitudes of the band boundaries (degrees), south pole to north pole.
    pub latitude: Vec<f64>,
    /// σ at the bottom of each layer, top first; the last is 1 (the surface).
    pub sigma: Vec<f64>,
    /// `streamfunction[k][b]`: northward mass flux (kg/s) above σ `sigma[k]`
    /// across boundary `b`. Positive is northward aloft and southward below,
    /// so the northern Hadley cell is positive and the southern negative.
    /// The last row (the surface) is the mass-closure residual, near zero
    /// in a time mean.
    pub streamfunction: Vec<Vec<f64>>,
}

impl Overturning {
    /// The index of the row at the interior interface nearest `sigma` (the
    /// upper one on a tie); the surface row only if there is one layer.
    pub fn nearest_row(&self, sigma: f64) -> usize {
        let interior = self.sigma.len().saturating_sub(1).max(1);
        (0..interior)
            .min_by(|&a, &b| {
                (self.sigma[a] - sigma)
                    .abs()
                    .total_cmp(&(self.sigma[b] - sigma).abs())
            })
            .unwrap_or(0)
    }

    /// The streamfunction at the interior interface nearest `sigma` (never
    /// the surface row, unless there is only one layer).
    pub fn at(&self, sigma: f64) -> &[f64] {
        &self.streamfunction[self.nearest_row(sigma)]
    }

    /// The largest |ψ| at the surface (kg/s): how far the time mean is from
    /// closing the column's mass budget.
    pub fn residual(&self) -> f64 {
        self.streamfunction
            .last()
            .map_or(0.0, |row| row.iter().fold(0.0f64, |m, v| m.max(v.abs())))
    }

    /// The north–south mirror image: each row reversed and negated
    /// (latitudes keep their order), to match [`ZonalMeans::mirrored`].
    pub fn mirrored(&self) -> Self {
        Self {
            latitude: self.latitude.clone(),
            sigma: self.sigma.clone(),
            streamfunction: self
                .streamfunction
                .iter()
                .map(|row| row.iter().rev().map(|v| -v).collect())
                .collect(),
        }
    }
}

/// One overturning cell.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Cell {
    /// Southern edge (degrees): a zero of the streamfunction, or −90.
    pub south: f64,
    /// Northern edge (degrees): a zero of the streamfunction, or +90.
    pub north: f64,
    /// Latitude of the streamfunction's extremum (degrees).
    pub latitude: f64,
    /// The signed extremum of the streamfunction (kg/s): positive turns
    /// clockwise seen from the east (northward aloft).
    pub strength: f64,
}

/// A maximum of the zonal-mean eastward wind.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Jet {
    /// Latitude (degrees).
    pub latitude: f64,
    /// Eastward wind (m/s).
    pub speed: f64,
}

/// The overturning cells of a streamfunction profile `psi` on `latitude`
/// (degrees, increasing), south to north.
///
/// Each maximal run of one sign is a cell (exact zeros belong to none),
/// with its edges at the zero crossings, interpolated linearly, or at ±90°
/// where it reaches the end of the profile. Runs weaker than [`CELL_FLOOR`]
/// of the largest |ψ| are noise and dropped; where a dropped run sat
/// between two cells of the same sign, those two are one cell.
pub fn cells(latitude: &[f64], psi: &[f64]) -> Vec<Cell> {
    assert_eq!(latitude.len(), psi.len(), "one ψ per latitude");
    let largest = psi.iter().fold(0.0f64, |m, v| m.max(v.abs()));
    if largest == 0.0 {
        return Vec::new();
    }
    let zero = |i: usize, j: usize| {
        let (a, b) = (psi[i], psi[j]);
        latitude[i] + (latitude[j] - latitude[i]) * a / (a - b)
    };
    let last = psi.len() - 1;
    // Every run of one sign, with whether it clears the floor.
    let mut runs: Vec<(Cell, bool)> = Vec::new();
    let mut i = 0;
    while i < psi.len() {
        if psi[i] == 0.0 {
            i += 1;
            continue;
        }
        let positive = psi[i] > 0.0;
        let start = i;
        while i < psi.len() && psi[i] != 0.0 && (psi[i] > 0.0) == positive {
            i += 1;
        }
        let end = i - 1;
        let peak = (start..=end)
            .max_by(|&a, &b| psi[a].abs().total_cmp(&psi[b].abs()))
            .expect("a run is not empty");
        let cell = Cell {
            south: if start == 0 {
                -90.0
            } else {
                zero(start - 1, start)
            },
            north: if end == last {
                90.0
            } else {
                zero(end, end + 1)
            },
            latitude: latitude[peak],
            strength: psi[peak],
        };
        runs.push((cell, psi[peak].abs() >= CELL_FLOOR * largest));
    }
    let mut out: Vec<Cell> = Vec::new();
    let mut dropped_since_last = false;
    for (cell, kept) in runs {
        if !kept {
            dropped_since_last = true;
            continue;
        }
        match out.last_mut() {
            Some(prev)
                if dropped_since_last && prev.strength.signum() == cell.strength.signum() =>
            {
                prev.north = cell.north;
                if cell.strength.abs() > prev.strength.abs() {
                    prev.latitude = cell.latitude;
                    prev.strength = cell.strength;
                }
            }
            _ => out.push(cell),
        }
        dropped_since_last = false;
    }
    out
}

/// Where the trade winds end in one hemisphere, as a positive latitude
/// (degrees), from band-centre `latitude` and the lowest layer's eastward
/// `wind`.
///
/// Walks poleward from the equator, past any westerlies there, to the first
/// easterlies (u < 0), and returns where the wind turns westerly again
/// (u ≥ 0), interpolated linearly. `None` if there are no easterlies or they
/// reach the pole.
pub fn trades_edge(latitude: &[f64], wind: &[f64], north: bool) -> Option<f64> {
    assert_eq!(latitude.len(), wind.len(), "one wind per latitude");
    let poleward: Vec<usize> = if north {
        (0..latitude.len())
            .filter(|&i| latitude[i] >= 0.0)
            .collect()
    } else {
        (0..latitude.len())
            .rev()
            .filter(|&i| latitude[i] <= 0.0)
            .collect()
    };
    let first = poleward.iter().position(|&i| wind[i] < 0.0)?;
    poleward[first..].windows(2).find_map(|w| {
        let (i, j) = (w[0], w[1]);
        (wind[j] >= 0.0).then(|| {
            let (a, b) = (latitude[i].abs(), latitude[j].abs());
            a + (b - a) * wind[i] / (wind[i] - wind[j])
        })
    })
}

/// The jets of an eastward-wind profile on `latitude` (degrees, increasing),
/// south to north.
///
/// A jet is a local maximum (a plateau counts once) of westerly wind whose
/// topographic prominence — its height above the higher of the lowest points
/// separating it from higher ground on either side, or above the profile's
/// minimum for the highest peak — is at least [`JET_PROMINENCE`] and
/// [`JET_FRACTION`] of the strongest westerly. Latitude and speed are refined
/// by a parabola through the maximum and its neighbours.
pub fn jets(latitude: &[f64], wind: &[f64]) -> Vec<Jet> {
    assert_eq!(latitude.len(), wind.len(), "one wind per latitude");
    let n = wind.len();
    let strongest = wind.iter().fold(0.0f64, |m, &v| m.max(v));
    if strongest <= 0.0 {
        return Vec::new();
    }
    let floor = JET_PROMINENCE.max(JET_FRACTION * strongest);
    let lowest = wind.iter().fold(f64::INFINITY, |m, &v| m.min(v));
    let mut out = Vec::new();
    let mut i = 0;
    while i < n {
        // The plateau [i, j] of equal values.
        let mut j = i;
        while j + 1 < n && wind[j + 1] == wind[i] {
            j += 1;
        }
        let h = wind[i];
        let left_lower = i == 0 || wind[i - 1] < h;
        let right_lower = j == n - 1 || wind[j + 1] < h;
        if h > 0.0 && left_lower && right_lower {
            let base = match (col(wind, h, (0..i).rev()), col(wind, h, j + 1..n)) {
                (Some(a), Some(b)) => a.max(b),
                (Some(a), None) | (None, Some(a)) => a,
                (None, None) => lowest,
            };
            if h - base >= floor {
                out.push(refine(latitude, wind, i, j));
            }
        }
        i = j + 1;
    }
    out
}

/// The lowest point between a peak of height `h` and the first higher
/// point along `path`, or `None` if there is no higher ground that way.
fn col(wind: &[f64], h: f64, mut path: impl Iterator<Item = usize>) -> Option<f64> {
    let mut low = h;
    path.find(|&k| {
        low = low.min(wind[k]);
        wind[k] > h
    })
    .map(|_| low)
}

/// The peak of the maximum spanning `[i, j]`: a parabola through a single
/// interior maximum and its neighbours; the point itself at an end; the
/// plateau's middle for a plateau.
fn refine(latitude: &[f64], wind: &[f64], i: usize, j: usize) -> Jet {
    if i != j {
        return Jet {
            latitude: 0.5 * (latitude[i] + latitude[j]),
            speed: wind[i],
        };
    }
    if i == 0 || i == wind.len() - 1 {
        return Jet {
            latitude: latitude[i],
            speed: wind[i],
        };
    }
    let (x0, x1, x2) = (latitude[i - 1], latitude[i], latitude[i + 1]);
    let (y0, y1, y2) = (wind[i - 1], wind[i], wind[i + 1]);
    // y = y1 + b (x − x1) + c (x − x1)², from divided differences.
    let (d0, d2) = (x0 - x1, x2 - x1);
    let (s0, s2) = ((y0 - y1) / d0, (y2 - y1) / d2);
    let c = (s2 - s0) / (d2 - d0);
    let b = s0 - c * d0;
    if c >= 0.0 {
        return Jet {
            latitude: x1,
            speed: y1,
        };
    }
    let dx = (-b / (2.0 * c)).clamp(d0, d2);
    Jet {
        latitude: x1 + dx,
        speed: y1 + b * dx + c * dx * dx,
    }
}

/// The equator-to-pole contrast of a band profile: the plain mean over bands
/// within 10° of the equator minus the mean over bands within 10° of either
/// pole. Without a band in a range, the band nearest the equator (or each
/// pole) stands in.
pub fn equator_pole_contrast(latitude: &[f64], values: &[f64]) -> f64 {
    assert_eq!(latitude.len(), values.len(), "one value per latitude");
    assert!(!latitude.is_empty(), "need at least one band");
    let mean =
        |picked: Vec<usize>| picked.iter().map(|&i| values[i]).sum::<f64>() / picked.len() as f64;
    let pick = |inside: &dyn Fn(f64) -> bool| -> Vec<usize> {
        (0..latitude.len())
            .filter(|&i| inside(latitude[i].abs()))
            .collect()
    };
    let mut tropics = pick(&|a| a <= 10.0);
    if tropics.is_empty() {
        // The nearest band, or the two either side of the equator.
        let nearest = latitude.iter().fold(f64::INFINITY, |m, l| m.min(l.abs()));
        tropics = pick(&|a| a == nearest);
    }
    let mut poles = pick(&|a| a >= 80.0);
    if poles.is_empty() {
        let south = (0..latitude.len())
            .min_by(|&a, &b| latitude[a].total_cmp(&latitude[b]))
            .expect("not empty");
        let north = (0..latitude.len())
            .max_by(|&a, &b| latitude[a].total_cmp(&latitude[b]))
            .expect("not empty");
        poles.push(south);
        if north != south {
            poles.push(north);
        }
    }
    mean(tropics) - mean(poles)
}

/// The general circulation of a run's time means, read from its overturning
/// and zonal winds.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Circulation {
    /// Overturning cells at mid-troposphere (the interface nearest σ = 0.5), south to north.
    pub cells: Vec<Cell>,
    /// Jets of the upper layer, south to north.
    pub jets: Vec<Jet>,
    /// Maxima of the lowest layer's westerlies, south to north.
    pub westerlies: Vec<Jet>,
    /// Where the lowest layer's trades end, [south, north], as positive latitudes.
    pub trades_edge: [Option<f64>; 2],
}

impl Circulation {
    /// The circulation of `means` (every layer) and `overturning`; `upper`
    /// and `lower` are the layer indices whose jets and westerlies (and
    /// trades) are read.
    pub fn of(means: &ZonalMeans, overturning: &Overturning, upper: usize, lower: usize) -> Self {
        let lat = &means.latitude;
        Self {
            cells: cells(&overturning.latitude, overturning.at(0.5)),
            jets: jets(lat, &means.wind[upper]),
            westerlies: jets(lat, &means.wind[lower]),
            trades_edge: [false, true].map(|north| trades_edge(lat, &means.wind[lower], north)),
        }
    }

    /// The thermally direct (Hadley) cell of one hemisphere: the strongest
    /// positive cell reaching north of the equator, or the strongest
    /// negative one reaching south of it. Its poleward edge is the Hadley
    /// cell's edge (±90 if it reaches the pole).
    pub fn hadley(&self, north: bool) -> Option<&Cell> {
        let candidates = self.cells.iter().filter(|c| {
            if north {
                c.strength > 0.0 && c.north > 0.0
            } else {
                c.strength < 0.0 && c.south < 0.0
            }
        });
        candidates.max_by(|a, b| a.strength.abs().total_cmp(&b.strength.abs()))
    }

    /// The poleward edge of one hemisphere's Hadley cell, as a positive
    /// latitude (degrees).
    pub fn hadley_edge(&self, north: bool) -> Option<f64> {
        self.hadley(north)
            .map(|c| if north { c.north } else { -c.south })
    }

    /// The number of cells whose extremum is in that hemisphere.
    pub fn cells_in(&self, north: bool) -> usize {
        self.cells
            .iter()
            .filter(|c| {
                if north {
                    c.latitude > 0.0
                } else {
                    c.latitude < 0.0
                }
            })
            .count()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `n + 1` boundary latitudes from pole to pole.
    fn boundaries(n: usize) -> Vec<f64> {
        (0..=n)
            .map(|b| -90.0 + b as f64 * 180.0 / n as f64)
            .collect()
    }

    /// `n` band-centre latitudes.
    fn centres(n: usize) -> Vec<f64> {
        (0..n)
            .map(|b| -90.0 + (b as f64 + 0.5) * 180.0 / n as f64)
            .collect()
    }

    fn profile(lat: &[f64], f: impl Fn(f64) -> f64) -> Vec<f64> {
        lat.iter().map(|l| f(l.to_radians())).collect()
    }

    #[test]
    fn cos_sin_6_has_six_alternating_cells() {
        let lat = boundaries(40);
        let psi = profile(&lat, |p| fms_math::cos(p) * fms_math::sin(6.0 * p));
        let c = cells(&lat, &psi);
        assert_eq!(c.len(), 6, "{c:?}");
        let edges = [-90.0, -60.0, -30.0, 0.0, 30.0, 60.0, 90.0];
        for (i, cell) in c.iter().enumerate() {
            assert!((cell.south - edges[i]).abs() < 0.5, "{cell:?}");
            assert!((cell.north - edges[i + 1]).abs() < 0.5, "{cell:?}");
            // Alternating, starting negative: sin 6φ < 0 just north of −90°.
            let sign = if i % 2 == 0 { -1.0 } else { 1.0 };
            assert_eq!(cell.strength.signum(), sign, "{cell:?}");
        }
        let circ = Circulation {
            cells: c,
            ..Circulation::default()
        };
        assert!((circ.hadley_edge(true).unwrap() - 30.0).abs() < 0.5);
        assert!((circ.hadley_edge(false).unwrap() - 30.0).abs() < 0.5);
        assert!(circ.hadley(true).unwrap().strength > 0.0);
        assert!(circ.hadley(false).unwrap().strength < 0.0);
        assert_eq!(circ.cells_in(true), 3);
        assert_eq!(circ.cells_in(false), 3);
    }

    #[test]
    fn weak_runs_are_dropped_and_their_neighbours_joined() {
        let lat = boundaries(8);
        // Poles, a strong positive run, a weak negative one, another strong
        // positive one, a strong negative one and a weak positive one at the
        // northern end.
        let psi = [0.0, 5.0, 10.0, -0.1, 8.0, -6.0, -6.0, 0.2, 0.1];
        let c = cells(&lat, &psi);
        assert_eq!(c.len(), 2, "{c:?}");
        // The joined cell runs from the pole to the zero between 8 and −6
        // and keeps the larger extremum.
        assert_eq!(c[0].south, -90.0);
        assert_eq!(c[0].strength, 10.0);
        assert_eq!(c[0].latitude, lat[2]);
        let zero = lat[4] + (lat[5] - lat[4]) * 8.0 / 14.0;
        assert!((c[0].north - zero).abs() < 1e-12);
        // The weak run at the end is just dropped.
        assert_eq!(c[1].strength, -6.0);
        assert!((c[1].south - zero).abs() < 1e-12);
        let end = lat[6] + (lat[7] - lat[6]) * 6.0 / 6.2;
        assert!((c[1].north - end).abs() < 1e-12);
        // All zero: no cells.
        assert!(cells(&lat, &[0.0; 9]).is_empty());
    }

    #[test]
    fn one_cell_per_hemisphere_reaches_the_pole() {
        let lat = boundaries(36);
        let psi = profile(&lat, |p| fms_math::cos(p) * fms_math::sin(2.0 * p));
        let circ = Circulation {
            cells: cells(&lat, &psi),
            ..Circulation::default()
        };
        assert_eq!(circ.cells.len(), 2, "{:?}", circ.cells);
        assert_eq!(circ.hadley_edge(true), Some(90.0));
        assert_eq!(circ.hadley_edge(false), Some(90.0));
    }

    #[test]
    fn the_trades_end_where_the_easterlies_do() {
        let lat = centres(36);
        let u = profile(&lat, |p| -fms_math::cos(3.0 * p));
        for north in [true, false] {
            let edge = trades_edge(&lat, &u, north).unwrap();
            assert!((edge - 30.0).abs() < 0.5, "{edge}");
        }
        // Westerlies at the equator, easterlies from 10° to 40°, westerlies
        // beyond.
        let u: Vec<f64> = lat
            .iter()
            .map(|l| {
                if (10.0..40.0).contains(&l.abs()) {
                    -2.0
                } else {
                    3.0
                }
            })
            .collect();
        for north in [true, false] {
            // Between the bands at 37.5° (−2) and 42.5° (3).
            let edge = trades_edge(&lat, &u, north).unwrap();
            assert!((edge - (37.5 + 5.0 * 2.0 / 5.0)).abs() < 1e-12, "{edge}");
        }
        // No easterlies; easterlies to the pole.
        assert_eq!(trades_edge(&lat, &[1.0; 36], true), None);
        assert_eq!(trades_edge(&lat, &[-1.0; 36], false), None);
    }

    #[test]
    fn jets_are_found_at_their_latitudes() {
        let lat = centres(90);
        // Per hemisphere: a 30 m/s jet at 48.7° and a 12 m/s one at 26.3°
        // (south: 25 m/s at 51.2°, 9 m/s at 23.8°), each 5° wide.
        let peaks = [(-51.2, 25.0), (-23.8, 9.0), (26.3, 12.0), (48.7, 30.0)];
        let gauss = |l: f64, (c, s): (f64, f64)| s * fms_math::exp(-((l - c) / 5.0).powi(2));
        let u: Vec<f64> = lat
            .iter()
            .map(|&l| peaks.iter().map(|&p| gauss(l, p)).sum::<f64>() - 1.0)
            .collect();
        let found = jets(&lat, &u);
        assert_eq!(found.len(), 4, "{found:?}");
        for (jet, (c, s)) in found.iter().zip(peaks) {
            assert!((jet.latitude - c).abs() < 0.5, "{jet:?} vs {c}");
            assert!((jet.speed - (s - 1.0)).abs() < 0.5, "{jet:?} vs {s}");
        }
    }

    #[test]
    fn ripples_and_plateaus_are_not_extra_jets() {
        let lat = centres(90);
        // A 30 m/s jet with a 1 m/s ripple on it.
        let u = profile(&lat, |p| {
            30.0 * fms_math::exp(-((p.to_degrees() - 45.0) / 15.0).powi(2))
                + fms_math::sin(p * 40.0)
        });
        let found = jets(&lat, &u);
        assert_eq!(found.len(), 1, "{found:?}");
        assert!((found[0].latitude - 45.0).abs() < 2.0, "{found:?}");
        // A flat top three bands wide counts once, at its middle.
        let mut u = vec![0.0; 9];
        u[3..6].fill(10.0);
        let lat = centres(9);
        let found = jets(&lat, &u);
        assert_eq!(
            found,
            vec![Jet {
                latitude: lat[4],
                speed: 10.0
            }]
        );
    }

    #[test]
    fn the_equator_pole_contrast_of_a_known_profile() {
        let lat = centres(18);
        // T = 300 − 40 sin²φ: the bands at ±5° and ±85°.
        let t = profile(&lat, |p| 300.0 - 40.0 * fms_math::sin(p).powi(2));
        let s = |d: f64| fms_math::sin(d.to_radians()).powi(2);
        let expected = 40.0 * (s(85.0) - s(5.0));
        assert!((equator_pole_contrast(&lat, &t) - expected).abs() < 1e-12);
        // Too few bands for the ranges: the nearest bands stand in.
        let lat = centres(4);
        let t = [1.0, 5.0, 7.0, 3.0];
        assert_eq!(equator_pole_contrast(&lat, &t), 6.0 - 2.0);
    }

    #[test]
    fn mirroring_twice_is_the_identity_and_at_picks_the_nearest_interface() {
        let o = Overturning {
            latitude: boundaries(4),
            sigma: vec![0.25, 0.5, 0.85, 1.0],
            streamfunction: vec![
                vec![0.0, 1.0, 2.0, 3.0, 0.0],
                vec![0.0, 4.0, 5.0, 6.0, 0.0],
                vec![0.0, 7.0, 8.0, 9.0, 0.0],
                vec![0.0, 0.1, -0.3, 0.2, 0.0],
            ],
        };
        assert_eq!(o.mirrored().mirrored(), o);
        assert_eq!(
            o.mirrored().streamfunction[0],
            [-0.0, -3.0, -2.0, -1.0, -0.0]
        );
        assert_eq!(o.at(0.5), &o.streamfunction[1][..]);
        assert_eq!(o.at(0.3), &o.streamfunction[0][..]);
        // Never the surface, however near.
        assert_eq!(o.at(0.99), &o.streamfunction[2][..]);
        assert!((o.residual() - 0.3).abs() < 1e-15);
        let one = Overturning {
            sigma: vec![1.0],
            streamfunction: vec![vec![0.0, 2.0, 0.0]],
            latitude: boundaries(2),
        };
        assert_eq!(one.at(0.5), &one.streamfunction[0][..]);
    }
}
