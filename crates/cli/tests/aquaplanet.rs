//! The reference aquaplanet: `fms-climate climate` on a blank world (all
//! ocean, with the coupled ocean), run as shipped as the real binary at 642
//! cells for nine years, has an Earth-like general circulation and a rain
//! belt that follows the sun.

use std::collections::HashMap;
use std::path::PathBuf;
use std::process::Command;
use std::time::{SystemTime, UNIX_EPOCH};

/// A fresh directory under the system's temporary one, removed when
/// dropped.
struct TempDir(PathBuf);

impl TempDir {
    fn new(name: &str) -> Self {
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let dir = std::env::temp_dir().join(format!(
            "fms-climate-aquaplanet-{name}-{}-{nanos}",
            std::process::id()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        Self(dir)
    }
}

impl Drop for TempDir {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

/// `circulation.csv` as a map from measure to value (`None` where empty).
fn measures(csv: &str) -> HashMap<String, Option<f64>> {
    csv.lines()
        .skip(1)
        .map(|line| {
            let (name, value) = line.split_once(',').expect("measure,value");
            let value = (!value.is_empty()).then(|| value.parse().expect("a number"));
            (name.to_string(), value)
        })
        .collect()
}

#[test]
#[ignore = "slow: nine model years at 642 cells (3–7 minutes); run with --ignored"]
fn the_aquaplanet_has_an_earth_like_circulation_and_rain_belt() {
    let dir = TempDir::new("reference");
    let output = Command::new(env!("CARGO_BIN_EXE_fms-climate"))
        .args([
            "climate",
            "--terrain-frequency",
            "16",
            "--climate-frequency",
            "8",
            "--bands",
            "24",
            "--width",
            "256",
            "--frames",
            "2",
            "--movie-days",
            "2",
            "--out-dir",
        ])
        .arg(&dir.0)
        .output()
        .expect("the binary runs");
    assert!(
        output.status.success(),
        "fms-climate climate failed:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let path = dir.0.join("climate/circulation.csv");
    let csv = std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    let m = measures(&csv);
    let mut names: Vec<&String> = m.keys().collect();
    names.sort();
    for name in names {
        eprintln!("{name} = {:?}", m[name]);
    }
    let get = |name: &str| -> f64 {
        m.get(name)
            .copied()
            .flatten()
            .unwrap_or_else(|| panic!("no {name}"))
    };
    let within = |name: &str, v: f64, lo: f64, hi: f64| {
        assert!((lo..=hi).contains(&v), "{name} {v} is outside {lo} to {hi}");
    };

    // Thermally direct Hadley cells: measured +38.0e9 north, −40.6e9 south.
    let (north, south) = (
        get("hadley_strength_kg_s_north"),
        get("hadley_strength_kg_s_south"),
    );
    assert!(north > 0.0 && south < 0.0, "{north} {south}");
    within("northern Hadley cell", north.abs(), 20e9, 80e9);
    within("southern Hadley cell", south.abs(), 20e9, 80e9);

    // Hadley edges: measured +22.6° and −24.1°, 1.5° from mirroring.
    let (north, south) = (get("hadley_edge_deg_north"), get("hadley_edge_deg_south"));
    within("northern Hadley edge", north.abs(), 15.0, 32.0);
    within("southern Hadley edge", south.abs(), 15.0, 32.0);
    within("Hadley edges' asymmetry", (north + south).abs(), 0.0, 6.0);

    // Trades' edges: measured +25.4° and −24.3°.
    for side in ["north", "south"] {
        let name = format!("trades_edge_deg_{side}");
        within(&name, get(&name).abs(), 18.0, 35.0);
    }

    // Upper-layer jets: measured one each side, 25.6 m/s at +30.4° and
    // 25.5 m/s at −30.2°.
    for side in ["north", "south"] {
        assert!(get(&format!("jets_{side}")) >= 1.0, "no {side}ern jet");
        let latitude = get(&format!("jet_latitude_deg_{side}"));
        within(
            &format!("{side}ern jet latitude"),
            latitude.abs(),
            22.0,
            45.0,
        );
        let speed = get(&format!("jet_speed_m_s_{side}"));
        within(&format!("{side}ern jet speed"), speed, 15.0, 45.0);
    }

    // Equator-to-pole contrast: measured 16.3 K.
    within(
        "equator-pole contrast",
        get("equator_pole_contrast_k"),
        10.0,
        30.0,
    );

    // Rain belt: measured annual mean +0.12°, furthest north +2.26° in
    // month 5, furthest south −1.40° in month 2, swing 3.66°. Months count
    // from the northern spring equinox, so the northern summer solstice ends
    // month 3 and the southern one month 9: the belt is furthest into each
    // hemisphere in the half year after its solstice.
    within(
        "annual-mean rain belt",
        get("rain_belt_mean_deg").abs(),
        0.0,
        1.5,
    );
    let (north, north_month) = (get("rain_belt_north_deg"), get("rain_belt_north_month"));
    assert!(north >= 1.0, "the rain belt reaches only {north}° north");
    assert!(
        (4.0..=9.0).contains(&north_month),
        "furthest north in month {north_month}"
    );
    let (south, south_month) = (get("rain_belt_south_deg"), get("rain_belt_south_month"));
    assert!(south <= -1.0, "the rain belt reaches only {south}° south");
    assert!(
        [10.0, 11.0, 12.0, 1.0, 2.0, 3.0].contains(&south_month),
        "furthest south in month {south_month}"
    );
    assert!(get("rain_belt_swing_deg") >= 2.0);

    // Wetter tropics than subtropics: measured 4.14 against 3.69 mm/day
    // (1.12 times).
    let (tropical, subtropical) = (get("tropical_rain_mm_day"), get("subtropical_rain_mm_day"));
    assert!(
        tropical > 1.05 * subtropical,
        "tropical rain {tropical} against subtropical {subtropical}"
    );
}
