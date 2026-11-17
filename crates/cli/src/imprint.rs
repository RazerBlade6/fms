//! The grid-imprint outputs of `fms-climate held-suarez` and of
//! `fms-climate run --imprint-block`: `imprint.csv`, a
//! map of each field's time mean minus its zonal mean (`imprint.<field>.png`)
//! and of the part of it with the grid's symmetry
//! (`imprint.<field>.locked.png`), and the summary's lines on what the two
//! halves of the run share.

use std::fmt::Write as _;
use std::path::Path;

use fms_atmos::{GridImprint, Imprint};

use crate::ImageArgs;
use crate::basins::{DivergingMap, write_diverging_map};

/// The share of cells whose |residual| an imprint map's colours span; the
/// rest saturate, so a few strong cells don't wash out the pattern.
const MAP_QUANTILE: f64 = 0.99;

/// How a field is shown: its name in the summary, the factor from its
/// units to the shown ones, and the shown units.
fn shown(field: &str) -> (&'static str, f64, &'static str) {
    match field {
        "surface_pressure" => ("surface pressure", 0.01, "hPa"),
        "temperature_lower" => ("lowest-layer temperature", 1.0, "K"),
        "east_wind_upper" => ("upper-layer eastward wind", 1.0, "m/s"),
        "north_wind_lower" => ("lowest-layer northward wind", 1.0, "m/s"),
        "atm.precipitation" => ("precipitation", 1.0, "mm/day"),
        "atm.evaporation" => ("evaporation", 1.0, "mm/day"),
        "rad.surface_temperature" => ("surface temperature", 1.0, "K"),
        _ => ("field", 1.0, ""),
    }
}

/// Writes `imprint.csv` (each field's measures, in its own units), and an
/// `imprint.<field>.png` map of each field's residual and an
/// `imprint.<field>.locked.png` map of its symmetric part (on the same
/// colour range) to `dir`, the maps stamped `stamp`.
pub fn write(
    imprint: &GridImprint,
    results: &[Imprint],
    grid: &fms_sim::Grid,
    dir: &Path,
    image: &ImageArgs,
    stamp: &str,
) -> Result<(), String> {
    let mut csv = String::from(
        "field,units,samples_a,samples_b,rms,zonal,noise,covariance,persistent,correlation,\
         locked,covariance_even_even,covariance_even_odd,covariance_odd_even,covariance_odd_odd,\
         transient,pentagon,ring",
    );
    for m in 1..=fms_atmos::imprint::WAVENUMBERS {
        let _ = write!(csv, ",wave_{m}");
    }
    csv.push('\n');
    let [a, b] = imprint.counts();
    for r in results {
        let _ = write!(
            csv,
            "{},{},{a},{b},{},{},{},{},{},{},{},{},{},{},{},{},{},{}",
            r.name,
            r.units,
            r.rms,
            r.zonal,
            r.noise,
            r.covariance,
            r.persistent,
            r.correlation,
            r.locked,
            r.symmetric[0],
            r.symmetric[1],
            r.symmetric[2],
            r.symmetric[3],
            r.transient,
            r.pentagon,
            r.ring
        );
        for p in r.waves {
            let _ = write!(csv, ",{p}");
        }
        csv.push('\n');
    }
    let path = dir.join("imprint.csv");
    std::fs::write(&path, csv).map_err(|e| format!("cannot write {}: {e}", path.display()))?;
    println!("grid imprint -> {}", path.display());
    for r in results {
        let (_, scale, units) = shown(&r.name);
        let values: Vec<f64> = r.residual.iter().map(|v| v * scale).collect();
        let name = format!("imprint.{}", r.name);
        let map = DivergingMap {
            name: &name,
            units,
            note: ", time mean minus zonal mean",
            stamp,
            quantile: MAP_QUANTILE,
            cap: None,
        };
        let cap = write_diverging_map(grid, &values, &map, image, dir)?;
        let values: Vec<f64> = r.locked_residual.iter().map(|v| v * scale).collect();
        let name = format!("imprint.{}.locked", r.name);
        let map = DivergingMap {
            name: &name,
            units,
            // Short, so the title clears the stamp at the default width.
            note: ", grid-symmetric part",
            stamp,
            quantile: MAP_QUANTILE,
            cap: Some(cap),
        };
        write_diverging_map(grid, &values, &map, image, dir)?;
    }
    Ok(())
}

/// The summary's grid-imprint block for blocks of `block_days`: per field,
/// the shared pattern with the grid's symmetry ("locked") and the
/// cross-half covariance of each symmetry part it comes from, all the
/// halves share, the noise and the weather's spread, and the pentagons.
/// Says so instead if there are fewer than two complete blocks (`results`
/// is then not read).
pub fn summary(imprint: &GridImprint, results: &[Imprint], block_days: f64) -> String {
    let mut out = String::new();
    if imprint.complete_blocks() < 2 {
        let _ = writeln!(
            out,
            "\nGrid imprint: fewer than two full {block_days}-day blocks after spin-up \
             ({} samples, {} per block); not measured.",
            imprint.samples(),
            imprint.block()
        );
        return out;
    }
    let [a, b] = imprint.counts();
    let _ = writeln!(
        out,
        "\nGrid imprint (time mean minus its local zonal mean), halves of alternating \
         {block_days}-day blocks ({} complete; {a} and {b} samples). The grid repeats under a \
         half-turn about the axis and a mirror in the equator (to {:.1e} rad); \"locked\" is the \
         shared pattern with both symmetries beyond what the other three symmetry parts share:",
        imprint.complete_blocks(),
        imprint.analysis().symmetry_error()
    );
    for r in results {
        let (name, scale, units) = shown(&r.name);
        let v = |x: f64| format!("{} {units}", significant(x * scale));
        let c = |x: f64| format!("{:+.2e}", x * scale * scale);
        let squared = if units.contains('/') {
            format!("({units})²")
        } else {
            format!("{units}²")
        };
        let _ = writeln!(
            out,
            "  {name}: locked {}; covariance with both symmetries {} {squared} against {}, {}, {} \
             in the other parts; persistent {}, correlation {:+.2} between halves; noise {}; \
             transient s.d. {}; pentagons {}, their ring {}",
            v(r.locked),
            c(r.symmetric[0]),
            c(r.symmetric[1]),
            c(r.symmetric[2]),
            c(r.symmetric[3]),
            v(r.persistent),
            r.correlation,
            v(r.noise),
            v(r.transient),
            v(r.pentagon),
            v(r.ring)
        );
    }
    out
}

/// `x` to three significant figures, without an exponent.
fn significant(x: f64) -> String {
    if x == 0.0 || !x.is_finite() {
        return format!("{x}");
    }
    let digits = 2 - fms_math::log10(x.abs()).floor() as i32;
    format!("{x:.*}", digits.max(0) as usize)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn three_significant_figures() {
        assert_eq!(significant(0.012345), "0.0123");
        assert_eq!(significant(1.2345), "1.23");
        assert_eq!(significant(123.45), "123");
        assert_eq!(significant(12345.0), "12345");
        assert_eq!(significant(0.0), "0");
    }
}
