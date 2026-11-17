//! The circulation outputs shared by `fms-climate held-suarez` and
//! `fms-climate climate`: the overturning streamfunction as `overturning.csv`
//! and `overturning.png`, the summary's lines on the cells, trades and jets
//! (and, for `climate`, the rain belt), and the same measures as
//! `circulation.csv`.

use std::fmt::Write as _;
use std::path::Path;

use fms_atmos::{Cell, Circulation, Jet, Overturning, RAIN_BELT_LIMIT, belt_mean, rain_belt};

use crate::chart::{self, Series, TimeAxis};
use crate::render::write_png;

/// Streamfunctions are reported in this unit (kg/s).
const UNIT: f64 = 1e9;

/// Writes `overturning.csv` (kg/s) and `overturning.png` (10⁹ kg/s, one
/// panel per interior interface) to `dir`, the chart under `header`.
pub fn write_overturning(
    overturning: &Overturning,
    dir: &Path,
    width: usize,
    header: &str,
) -> Result<(), String> {
    let interior = overturning.sigma.len().saturating_sub(1);
    let mut csv = String::from("latitude_deg");
    for sigma in &overturning.sigma[..interior] {
        let _ = write!(csv, ",streamfunction_sigma_{sigma}_kg_s");
    }
    csv.push_str(",surface_residual_kg_s\n");
    for (b, latitude) in overturning.latitude.iter().enumerate() {
        let _ = write!(csv, "{latitude}");
        for row in &overturning.streamfunction {
            let _ = write!(csv, ",{}", row[b]);
        }
        csv.push('\n');
    }
    let path = dir.join("overturning.csv");
    std::fs::write(&path, csv).map_err(|e| format!("cannot write {}: {e}", path.display()))?;
    println!("overturning -> {}", path.display());

    let scaled: Vec<Vec<f64>> = overturning.streamfunction[..interior]
        .iter()
        .map(|row| row.iter().map(|v| v / UNIT).collect())
        .collect();
    let series: Vec<Series> = overturning.sigma[..interior]
        .iter()
        .zip(&scaled)
        .map(|(sigma, values)| Series {
            title: format!("OVERTURNING ABOVE SIGMA {sigma:.2} (1E9 KG/S)"),
            values,
        })
        .collect();
    let axis = TimeAxis {
        times: &overturning.latitude,
        caption: "LATITUDE (DEGREES)",
        hours: false,
    };
    let (pixels, height) = chart::render_chart(width, header, &axis, &series);
    let path = dir.join("overturning.png");
    write_png(&path, width, height, &pixels, &chart::palette())
        .map_err(|e| format!("cannot write {}: {e}", path.display()))?;
    println!("overturning chart -> {}", path.display());
    Ok(())
}

/// The summary's circulation block: per hemisphere, the Hadley cell, every
/// cell, where the trades end, the upper-layer jets and the lowest layer's
/// westerlies; the lowest layer's equator-to-pole temperature `contrast`
/// (K), and the rain belt and the seasonal range of surface temperature
/// (`seasons`), if given; and the overturning's mass-closure residual.
pub fn summary(
    circulation: &Circulation,
    overturning: &Overturning,
    contrast: Option<f64>,
    seasons: Option<(&RainBelt, &SeasonalRange)>,
) -> String {
    let mut out = String::new();
    let row = overturning.nearest_row(0.5);
    let _ = writeln!(
        out,
        "Overturning at σ {} (10⁹ kg/s, positive northward aloft):",
        overturning.sigma[row]
    );
    for (name, north) in [("North", true), ("South", false)] {
        let sign = if north { 1.0 } else { -1.0 };
        let hadley = match circulation.hadley(north) {
            Some(c) => {
                let edge = if north { c.north } else { c.south };
                let edge = if edge.abs() >= 90.0 {
                    "reaches the pole".into()
                } else {
                    format!("edge at {edge:.1}°")
                };
                format!(
                    "Hadley cell {:+.1} at {:.1}°, {edge}",
                    c.strength / UNIT,
                    c.latitude
                )
            }
            None => "no thermally direct cell".into(),
        };
        let cells = hemisphere_cells(circulation, north);
        let list = cells
            .iter()
            .map(|c| {
                format!(
                    "{:.1}° to {:.1}° ({:+.1})",
                    c.south,
                    c.north,
                    c.strength / UNIT
                )
            })
            .collect::<Vec<_>>()
            .join(", ");
        let trades = match circulation.trades_edge[usize::from(north)] {
            Some(edge) => format!("trades end at {:.1}°", sign * edge),
            None => "no trades' edge".into(),
        };
        let _ = writeln!(
            out,
            "{name}: {hadley}; {} cell{}{}{list}; {trades}",
            cells.len(),
            if cells.len() == 1 { "" } else { "s" },
            if cells.is_empty() { "" } else { ": " },
        );
        let _ = writeln!(
            out,
            "{name}: upper-layer jets {}; lowest-layer westerlies {}",
            jets(&circulation.jets, north),
            jets(&circulation.westerlies, north)
        );
    }
    if let Some(contrast) = contrast {
        let _ = writeln!(
            out,
            "Lowest-layer temperature equator − pole: {contrast:.1} K"
        );
    }
    if let Some((rain, range)) = seasons {
        out.push_str(&rain.summary());
        out.push_str(&range.summary());
    }
    let _ = writeln!(
        out,
        "Overturning mass-closure residual (largest at the surface): {:.3} × 10⁹ kg/s",
        overturning.residual() / UNIT
    );
    out
}

/// The cells whose extremum is in one hemisphere (none on the equator).
fn hemisphere_cells(circulation: &Circulation, north: bool) -> Vec<&Cell> {
    circulation
        .cells
        .iter()
        .filter(|c| (c.latitude > 0.0) == north && c.latitude != 0.0)
        .collect()
}

/// The jets in one hemisphere (an equatorial one in both).
fn hemisphere_jets(jets: &[Jet], north: bool) -> Vec<&Jet> {
    jets.iter()
        .filter(|j| j.latitude.abs() < 0.05 || (j.latitude > 0.0) == north)
        .collect()
}

/// The jets in one hemisphere (an equatorial one in both), as text.
fn jets(jets: &[Jet], north: bool) -> String {
    let list: Vec<String> = hemisphere_jets(jets, north)
        .iter()
        .map(|j| format!("{:.1} m/s at {:.1}°", j.speed, j.latitude))
        .collect();
    if list.is_empty() {
        "none".into()
    } else {
        list.join(", ")
    }
}

/// Writes `circulation.csv` to `dir`: the summary's circulation measures,
/// one `measure,value` row each at full precision (empty where there is
/// none), per hemisphere south then north, followed by `extra` rows.
pub fn write_measures(
    circulation: &Circulation,
    extra: &[(&str, Option<f64>)],
    dir: &Path,
) -> Result<(), String> {
    let mut csv = String::from("measure,value\n");
    let mut row = |name: &str, value: Option<f64>| {
        let _ = match value {
            Some(v) => writeln!(csv, "{name},{v}"),
            None => writeln!(csv, "{name},"),
        };
    };
    let strongest = |jets: &[Jet], north: bool| -> Option<Jet> {
        hemisphere_jets(jets, north)
            .into_iter()
            .copied()
            .max_by(|a, b| a.speed.total_cmp(&b.speed))
    };
    type Measure<'a> = (&'a str, &'a dyn Fn(bool) -> Option<f64>);
    let hadley = |north| circulation.hadley(north);
    let jet = |north| strongest(&circulation.jets, north);
    let westerlies = |north| strongest(&circulation.westerlies, north);
    let measures: [Measure; 11] = [
        ("hadley_strength_kg_s", &|n| hadley(n).map(|c| c.strength)),
        ("hadley_latitude_deg", &|n| hadley(n).map(|c| c.latitude)),
        ("hadley_edge_deg", &|n| {
            hadley(n).map(|c| if n { c.north } else { c.south })
        }),
        ("cells", &|n| {
            Some(hemisphere_cells(circulation, n).len() as f64)
        }),
        ("trades_edge_deg", &|n| {
            circulation.trades_edge[usize::from(n)].map(|e| if n { e } else { -e })
        }),
        ("jets", &|n| {
            Some(hemisphere_jets(&circulation.jets, n).len() as f64)
        }),
        ("jet_speed_m_s", &|n| jet(n).map(|j| j.speed)),
        ("jet_latitude_deg", &|n| jet(n).map(|j| j.latitude)),
        ("westerly_belts", &|n| {
            Some(hemisphere_jets(&circulation.westerlies, n).len() as f64)
        }),
        ("westerlies_speed_m_s", &|n| westerlies(n).map(|j| j.speed)),
        ("westerlies_latitude_deg", &|n| {
            westerlies(n).map(|j| j.latitude)
        }),
    ];
    for (name, value) in measures {
        for (side, north) in [("south", false), ("north", true)] {
            row(&format!("{name}_{side}"), value(north));
        }
    }
    for &(name, value) in extra {
        row(name, value);
    }
    let path = dir.join("circulation.csv");
    std::fs::write(&path, csv).map_err(|e| format!("cannot write {}: {e}", path.display()))?;
    println!("circulation measures -> {}", path.display());
    Ok(())
}

/// The months a rain belt is measured over.
pub const RAIN_BELT_MONTHS: usize = 12;

/// The rain belt over a run's last [`RAIN_BELT_MONTHS`] completed months,
/// from each month's zonal-mean precipitation.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct RainBelt {
    /// The months the run completed.
    pub completed: u64,
    /// The rain belt's latitude (degrees) in each of the last 12 months, by
    /// month number (1–12, from the northern spring equinox) in increasing
    /// order; empty if the run completed fewer.
    pub monthly: Vec<(u64, Option<f64>)>,
    /// The rain belt of the 12 months' mean precipitation (degrees).
    pub mean: Option<f64>,
    /// The 12 months' mean precipitation within 10° of the equator (mm/day).
    pub tropical: Option<f64>,
    /// The 12 months' mean precipitation from 20° to 30° (mm/day).
    pub subtropical: Option<f64>,
}

impl RainBelt {
    /// The rain belt of `months`, the run's last completed months (up to
    /// 12) as (month number, zonal-mean precipitation in mm/day on band
    /// centres at `latitude`, degrees), of `completed` in all.
    pub fn of(latitude: &[f64], months: &[(u64, Vec<f64>)], completed: u64) -> Self {
        if months.len() < RAIN_BELT_MONTHS {
            return Self {
                completed,
                ..Self::default()
            };
        }
        let months = &months[months.len() - RAIN_BELT_MONTHS..];
        let mut mean = vec![0.0; latitude.len()];
        for (_, profile) in months {
            for (m, p) in mean.iter_mut().zip(profile) {
                *m += p;
            }
        }
        for m in &mut mean {
            *m /= RAIN_BELT_MONTHS as f64;
        }
        let mut monthly: Vec<(u64, Option<f64>)> = months
            .iter()
            .map(|(month, profile)| (*month, rain_belt(latitude, profile)))
            .collect();
        monthly.sort_by_key(|&(month, _)| month);
        Self {
            completed,
            monthly,
            mean: rain_belt(latitude, &mean),
            tropical: belt_mean(latitude, &mean, 0.0, 10.0),
            subtropical: belt_mean(latitude, &mean, 20.0, 30.0),
        }
    }

    /// The month whose rain belt is furthest north (`north`) or south, and
    /// that latitude (degrees); the earliest month on a tie.
    pub fn extreme(&self, north: bool) -> Option<(u64, f64)> {
        let sign = if north { 1.0 } else { -1.0 };
        self.monthly
            .iter()
            .filter_map(|&(month, belt)| belt.map(|b| (month, b)))
            .fold(None, |best: Option<(u64, f64)>, (month, b)| match best {
                Some((_, a)) if sign * b <= sign * a => best,
                _ => Some((month, b)),
            })
    }

    /// How far the rain belt moves over the year (degrees): furthest north
    /// minus furthest south.
    pub fn swing(&self) -> Option<f64> {
        Some(self.extreme(true)?.1 - self.extreme(false)?.1)
    }

    /// The summary's lines on the rain belt.
    pub fn summary(&self) -> String {
        let degrees = |v: Option<f64>| v.map_or("none".into(), |v| format!("{v:+.1}°"));
        let rate = |v: Option<f64>| v.map_or("none".into(), |v| format!("{v:.1} mm/day"));
        let head = format!("Rain belt (centroid of rain within {RAIN_BELT_LIMIT}° of the equator)");
        if self.monthly.is_empty() {
            return format!(
                "{head}: not measured; the run completed {} month{}, and it needs {RAIN_BELT_MONTHS}\n",
                self.completed,
                if self.completed == 1 { "" } else { "s" }
            );
        }
        let extreme = |north: bool| match self.extreme(north) {
            Some((month, b)) => format!("{} in month {month}", degrees(Some(b))),
            None => "none".into(),
        };
        format!(
            "{head}: annual mean {}; furthest north {}, furthest south {} (months from the \
             northern spring equinox); swing {}\n\
             Rain within 10° of the equator {}; at 20° to 30° {}\n",
            degrees(self.mean),
            extreme(true),
            extreme(false),
            self.swing().map_or("none".into(), |s| format!("{s:.1}°")),
            rate(self.tropical),
            rate(self.subtropical),
        )
    }

    /// The rows of `circulation.csv` on the rain belt.
    pub fn measures(&self) -> Vec<(&'static str, Option<f64>)> {
        let (north, south) = (self.extreme(true), self.extreme(false));
        vec![
            ("rain_belt_mean_deg", self.mean),
            ("rain_belt_north_deg", north.map(|e| e.1)),
            ("rain_belt_north_month", north.map(|e| e.0 as f64)),
            ("rain_belt_south_deg", south.map(|e| e.1)),
            ("rain_belt_south_month", south.map(|e| e.0 as f64)),
            ("rain_belt_swing_deg", self.swing()),
            ("tropical_rain_mm_day", self.tropical),
            ("subtropical_rain_mm_day", self.subtropical),
        ]
    }

    /// Writes `rain_belt.png`: the rain belt's latitude in each of the 12
    /// months, under `header`. Nothing if it was not measured.
    pub fn write_chart(&self, dir: &Path, width: usize, header: &str) -> Result<(), String> {
        if self.monthly.is_empty() {
            return Ok(());
        }
        let months: Vec<f64> = self.monthly.iter().map(|m| m.0 as f64).collect();
        let belt: Vec<f64> = self
            .monthly
            .iter()
            .map(|m| m.1.unwrap_or(f64::NAN))
            .collect();
        let series = [Series {
            title: "LATITUDE OF THE RAIN BELT (DEGREES)".into(),
            values: &belt,
        }];
        let axis = TimeAxis {
            times: &months,
            caption: "MONTH (FROM THE NORTHERN SPRING EQUINOX)",
            hours: false,
        };
        let (pixels, height) = chart::render_chart(width, header, &axis, &series);
        let path = dir.join("rain_belt.png");
        write_png(&path, width, height, &pixels, &chart::palette())
            .map_err(|e| format!("cannot write {}: {e}", path.display()))?;
        println!("rain-belt chart -> {}", path.display());
        Ok(())
    }
}

/// How far from the equator (degrees) the seasonal range of surface
/// temperature is measured: band centres within it count.
pub const SEASONAL_RANGE_LIMIT: f64 = 60.0;

/// The largest seasonal range of surface temperature near the equator: over
/// a run's last [`RAIN_BELT_MONTHS`] completed months (those the rain belt
/// uses), for each latitude band whose centre is within
/// [`SEASONAL_RANGE_LIMIT`] of the equator, the warmest month's zonal-mean
/// surface temperature minus the coldest's; the largest of those, and where.
/// Without seasons it is near zero.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct SeasonalRange {
    /// The months the run completed.
    pub completed: u64,
    /// The largest range (K) and its band's latitude (degrees); `None` if
    /// the run completed fewer than 12 months.
    pub largest: Option<(f64, f64)>,
}

impl SeasonalRange {
    /// The range from `months`, the run's last completed months (up to 12;
    /// the last 12 are used) as zonal-mean surface temperature (K) on band
    /// centres at `latitude` (degrees), of `completed` in all. The first
    /// band (from the south) wins a tie.
    pub fn of(latitude: &[f64], months: &[Vec<f64>], completed: u64) -> Self {
        if months.len() < RAIN_BELT_MONTHS {
            return Self {
                completed,
                largest: None,
            };
        }
        let months = &months[months.len() - RAIN_BELT_MONTHS..];
        let mut largest: Option<(f64, f64)> = None;
        for (b, &lat) in latitude.iter().enumerate() {
            if lat.abs() > SEASONAL_RANGE_LIMIT {
                continue;
            }
            let (lo, hi) = months
                .iter()
                .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), m| {
                    (lo.min(m[b]), hi.max(m[b]))
                });
            let range = hi - lo;
            if largest.is_none_or(|(r, _)| range > r) {
                largest = Some((range, lat));
            }
        }
        Self { completed, largest }
    }

    /// The summary's line.
    pub fn summary(&self) -> String {
        let head = format!(
            "Largest seasonal range of surface temperature within {SEASONAL_RANGE_LIMIT}° of the equator"
        );
        match self.largest {
            Some((range, lat)) => format!("{head}: {range:.1} K at {lat:.1}°\n"),
            None => format!(
                "{head}: not measured; the run completed {} month{}, and it needs {RAIN_BELT_MONTHS}\n",
                self.completed,
                if self.completed == 1 { "" } else { "s" }
            ),
        }
    }

    /// The rows of `circulation.csv`.
    pub fn measures(&self) -> Vec<(&'static str, Option<f64>)> {
        vec![
            ("seasonal_temperature_range_k", self.largest.map(|l| l.0)),
            (
                "seasonal_temperature_range_latitude_deg",
                self.largest.map(|l| l.1),
            ),
        ]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// 2.5° bands from pole to pole.
    fn centres() -> Vec<f64> {
        (0..72).map(|b| -90.0 + (b as f64 + 0.5) * 2.5).collect()
    }

    #[test]
    fn the_rain_belt_follows_the_months_and_says_when_it_was_not_measured() {
        let lat = centres();
        // Rain in one band each month, moving north to +6.25° in month 7
        // and south to −8.75° in month 1 (and month 12, a tie).
        let band = |month: u64| match month {
            1 | 12 => -8.75,
            7 => 6.25,
            _ => 1.25,
        };
        let profile = |month: u64| -> Vec<f64> {
            lat.iter()
                .map(|&l| if l == band(month) { 3.0 } else { 0.0 })
                .collect()
        };
        // Fourteen months, from month 11: the last twelve are months 1 to 12.
        let months: Vec<(u64, Vec<f64>)> = (0..14)
            .map(|k| (10 + k) % 12 + 1)
            .map(|m| (m, profile(m)))
            .collect();
        let belt = RainBelt::of(&lat, &months, 14);
        assert_eq!(belt.monthly.len(), 12);
        assert_eq!(belt.monthly[0], (1, Some(-8.75)));
        assert_eq!(belt.extreme(true), Some((7, 6.25)));
        assert_eq!(belt.extreme(false), Some((1, -8.75)));
        assert_eq!(belt.swing(), Some(15.0));
        let measures = belt.measures();
        assert_eq!(measures[2], ("rain_belt_north_month", Some(7.0)));
        assert_eq!(measures[4], ("rain_belt_south_month", Some(1.0)));
        // Rain at ±10° counts as tropical; none falls from 20° to 30°.
        assert_eq!(belt.subtropical, Some(0.0));
        let summary = belt.summary();
        assert!(
            summary.contains("furthest north +6.2° in month 7, furthest south -8.8° in month 1"),
            "{summary}"
        );
        assert_eq!(summary.lines().count(), 2);

        // Eleven months are too few.
        let short = RainBelt::of(&lat, &months[..11], 11);
        assert!(short.monthly.is_empty() && short.mean.is_none());
        assert!(short.measures().iter().all(|(_, v)| v.is_none()));
        assert_eq!(
            short.summary(),
            "Rain belt (centroid of rain within 20° of the equator): not measured; the run \
             completed 11 months, and it needs 12\n"
        );
    }

    #[test]
    fn the_seasonal_range_is_the_largest_within_60_degrees() {
        let lat = centres();
        // Twelve months of a profile whose range grows away from the
        // equator to 8 K at ±58.75°, and is larger still beyond 60°; the
        // southern band comes first and wins the tie.
        let months: Vec<Vec<f64>> = (0..13)
            .map(|m| {
                let phase = fms_math::sin(m as f64 / 12.0 * std::f64::consts::TAU);
                lat.iter()
                    .map(|&l| {
                        280.0
                            + 0.5
                                * (8.0 * l.abs() / 58.75)
                                * phase
                                * if l.abs() > 60.0 { 3.0 } else { 1.0 }
                    })
                    .collect()
            })
            .collect();
        let range = SeasonalRange::of(&lat, &months, 13);
        let (r, at) = range.largest.unwrap();
        assert_eq!(at, -58.75);
        // The sine reaches ±1 at months 3 and 9 of the last twelve.
        assert!((r - 8.0).abs() < 1e-9, "{r}");
        assert_eq!(
            range.summary(),
            "Largest seasonal range of surface temperature within 60° of the equator: 8.0 K at -58.8°\n"
        );
        assert_eq!(
            range.measures()[0],
            ("seasonal_temperature_range_k", Some(r))
        );
        let short = SeasonalRange::of(&lat, &months[..11], 11);
        assert!(short.measures().iter().all(|(_, v)| v.is_none()));
        assert!(short.summary().contains("not measured"));
    }
}
