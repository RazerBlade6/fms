//! The circulation outputs shared by `fms-climate held-suarez` and
//! `fms-climate climate`: the overturning streamfunction as `overturning.csv`
//! and `overturning.png`, and the summary's lines on the cells, trades and
//! jets.

use std::fmt::Write as _;
use std::path::Path;

use fms_atmos::{Cell, Circulation, Jet, Overturning};

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
/// (K) if given; and the overturning's mass-closure residual.
pub fn summary(
    circulation: &Circulation,
    overturning: &Overturning,
    contrast: Option<f64>,
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
        let cells: Vec<&Cell> = circulation
            .cells
            .iter()
            .filter(|c| (c.latitude > 0.0) == north && c.latitude != 0.0)
            .collect();
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
    let _ = writeln!(
        out,
        "Overturning mass-closure residual (largest at the surface): {:.3} × 10⁹ kg/s",
        overturning.residual() / UNIT
    );
    out
}

/// The jets in one hemisphere (an equatorial one in both), as text.
fn jets(jets: &[Jet], north: bool) -> String {
    let list: Vec<String> = jets
        .iter()
        .filter(|j| j.latitude.abs() < 0.05 || (j.latitude > 0.0) == north)
        .map(|j| format!("{:.1} m/s at {:.1}°", j.speed, j.latitude))
        .collect();
    if list.is_empty() {
        "none".into()
    } else {
        list.join(", ")
    }
}
