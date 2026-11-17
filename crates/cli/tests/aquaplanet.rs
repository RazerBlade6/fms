//! The reference aquaplanet: `fms-climate run` on a blank world (all
//! ocean, with the coupled ocean), run as shipped as the real binary at 642
//! cells for nine years, has an Earth-like general circulation and a rain
//! belt that follows the sun. Without tilt, on a circular orbit, it has no
//! seasons, mirrored hemispheres and its rain belt on the equator. With a
//! day twice as long, its Hadley cells are wider and stronger, its jets
//! further poleward and its equator-to-pole contrast smaller. With a day half
//! as long (at 2,562 cells, which its smaller storms need), its Hadley cells
//! are narrower and weaker, its contrast larger, and each hemisphere has a
//! second belt of westerlies.

use std::collections::HashMap;
use std::path::PathBuf;
use std::process::Command;
use std::sync::OnceLock;
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

/// The grid of a run: 642 cells with 24 latitude bands, or 2,562 cells with
/// the default 32.
#[derive(Clone, Copy)]
enum Grid {
    Cells642,
    Cells2562,
}

/// Runs `fms-climate run` as shipped at 642 cells, with 24 latitude
/// bands and small images, and `extra` options, into a fresh directory
/// named after `name`; prints every measure of `circulation.csv` and
/// returns the directory and the measures.
fn run_aquaplanet(name: &str, extra: &[&str]) -> (TempDir, HashMap<String, Option<f64>>) {
    run_on(Grid::Cells642, name, extra)
}

/// As [`run_aquaplanet`], on `grid`.
fn run_on(grid: Grid, name: &str, extra: &[&str]) -> (TempDir, HashMap<String, Option<f64>>) {
    let dir = TempDir::new(name);
    let grid: &[&str] = match grid {
        Grid::Cells642 => &[
            "--terrain-frequency",
            "16",
            "--climate-frequency",
            "8",
            "--bands",
            "24",
        ],
        Grid::Cells2562 => &["--terrain-frequency", "32", "--climate-frequency", "16"],
    };
    let output = Command::new(env!("CARGO_BIN_EXE_fms-climate"))
        .arg("run")
        .args(grid)
        .args(extra)
        .args([
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
        "fms-climate run failed:\n{}",
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
    (dir, m)
}

/// The value of measure `name`, which must be there.
#[track_caller]
fn get(m: &HashMap<String, Option<f64>>, name: &str) -> f64 {
    m.get(name)
        .copied()
        .flatten()
        .unwrap_or_else(|| panic!("no {name}"))
}

/// Asserts that measure `name`'s value `v` is within `lo..=hi`.
#[track_caller]
fn within(name: &str, v: f64, lo: f64, hi: f64) {
    assert!((lo..=hi).contains(&v), "{name} {v} is outside {lo} to {hi}");
}

/// The reference aquaplanet's measures: one run, shared by the tests here.
fn reference() -> &'static HashMap<String, Option<f64>> {
    static REFERENCE: OnceLock<HashMap<String, Option<f64>>> = OnceLock::new();
    REFERENCE.get_or_init(|| run_aquaplanet("reference", &[]).1)
}

#[test]
#[ignore = "slow: nine model years at 642 cells (3–7 minutes); run with --ignored"]
fn the_aquaplanet_has_an_earth_like_circulation_and_rain_belt() {
    let m = reference();
    let get = |name: &str| get(m, name);

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

#[test]
#[ignore = "slow: nine model years at 642 cells (3–7 minutes); run with --ignored"]
fn without_tilt_the_aquaplanet_has_no_seasons_and_mirrored_hemispheres() {
    let (dir, m) = run_aquaplanet("no-tilt", &["--obliquity", "0", "--eccentricity", "0"]);
    let get = |name: &str| get(&m, name);

    // The final year's zonal-mean surface temperature against its mirror
    // image, printed before any check fails.
    let path = dir.0.join("climate/zonal_mean.csv");
    let csv = std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    let mut lines = csv.lines();
    let header: Vec<&str> = lines.next().expect("a header").split(',').collect();
    let column = |name: &str| {
        header
            .iter()
            .position(|h| *h == name)
            .unwrap_or_else(|| panic!("no column {name}"))
    };
    let (lat, ts) = (column("latitude_deg"), column("rad.surface_temperature"));
    let rows: Vec<(f64, f64)> = lines
        .map(|line| {
            let cells: Vec<f64> = line.split(',').map(|c| c.parse().unwrap()).collect();
            (cells[lat], cells[ts])
        })
        .collect();
    let mut largest = 0.0f64;
    for (i, &(latitude, t)) in rows.iter().enumerate() {
        let (mirror, t_mirror) = rows[rows.len() - 1 - i];
        assert!(
            (latitude + mirror).abs() < 1e-6,
            "bands {latitude} and {mirror}"
        );
        if latitude.abs() <= 60.0 {
            eprintln!("surface temperature at {latitude:+.2}°: {t:.3} K, mirrored {t_mirror:.3} K");
            largest = largest.max((t - t_mirror).abs());
        }
    }
    eprintln!("surface temperature's largest north-south difference within 60°: {largest:.3} K");

    // No seasons: measured a largest range of 0.37 K (at −3.75°) and a
    // swing of 1.06° (+0.29° in month 4 to −0.78° in month 6).
    within(
        "seasonal range of surface temperature",
        get("seasonal_temperature_range_k"),
        0.0,
        1.0,
    );
    within("rain belt's swing", get("rain_belt_swing_deg"), 0.0, 1.5);

    // Mirrored Hadley cells: measured edges +23.3° and −26.6°, 3.25° from
    // mirroring (the bands are 7.5° wide; at 2,562 cells it is 1.1°);
    // strengths +52.8e9 and −49.6e9 kg/s, 6.1% apart.
    let (north, south) = (get("hadley_edge_deg_north"), get("hadley_edge_deg_south"));
    within("Hadley edges' asymmetry", (north + south).abs(), 0.0, 5.0);
    let (north, south) = (
        get("hadley_strength_kg_s_north").abs(),
        get("hadley_strength_kg_s_south").abs(),
    );
    within(
        "Hadley strengths' difference, of the stronger",
        (north - south).abs() / north.max(south),
        0.0,
        0.25,
    );

    // Mirrored jets: measured 30.7 m/s at +31.6° and 30.0 m/s at −32.0°,
    // 0.41° from mirroring and 0.68 m/s apart.
    let (north, south) = (get("jet_latitude_deg_north"), get("jet_latitude_deg_south"));
    within("jets' asymmetry", (north + south).abs(), 0.0, 4.0);
    let (north, south) = (get("jet_speed_m_s_north"), get("jet_speed_m_s_south"));
    within("jets' speed difference", (north - south).abs(), 0.0, 4.0);

    // Mirrored surface temperature within 60° of the equator, from
    // `zonal_mean.csv`'s `rad.surface_temperature` (the final year's zonal
    // mean): measured at most 0.17 K apart (at ±41.25°).
    within("surface temperature's asymmetry", largest, 0.0, 1.0);

    // The rain belt on the equator: measured an annual mean of −0.20°, and
    // 4.38 against 3.70 mm/day (1.18 times).
    within(
        "annual-mean rain belt",
        get("rain_belt_mean_deg").abs(),
        0.0,
        1.0,
    );
    let (tropical, subtropical) = (get("tropical_rain_mm_day"), get("subtropical_rain_mm_day"));
    assert!(
        tropical > 1.05 * subtropical,
        "tropical rain {tropical} against subtropical {subtropical}"
    );

    // The circulation is still there: measured Hadley edges +23.3° and
    // −26.6°, one jet each side.
    for side in ["north", "south"] {
        let edge = get(&format!("hadley_edge_deg_{side}"));
        within(&format!("{side}ern Hadley edge"), edge.abs(), 15.0, 32.0);
        assert!(get(&format!("jets_{side}")) >= 1.0, "no {side}ern jet");
    }
}

#[test]
#[ignore = "slow: nine model years at 642 cells, and the reference run if no other test here has made it (3–14 minutes); run with --ignored"]
fn a_longer_day_widens_the_hadley_cells_and_moves_the_jets_poleward() {
    // A 48-hour day in a year of the same length in seconds.
    let (_dir, slow) = run_aquaplanet("slow", &["--solar-day", "48", "--year-days", "182.6211"]);
    let reference = reference();
    for side in ["north", "south"] {
        let both = |name: &str| {
            let name = format!("{name}_{side}");
            (get(&slow, &name).abs(), get(reference, &name).abs())
        };

        // Hadley edges: measured +35.7° and −35.8° against +22.6° and −24.1°.
        let (edge, was) = both("hadley_edge_deg");
        assert!(
            edge >= was + 5.0,
            "{side}ern Hadley edge {edge}° against {was}°"
        );

        // Hadley strengths: measured +83.8e9 and −112.3e9 kg/s
        // against +38.0e9 and −40.6e9.
        let (strength, was) = both("hadley_strength_kg_s");
        assert!(
            strength >= 1.5 * was,
            "{side}ern Hadley cell {strength} against {was} kg/s"
        );

        // Upper-layer jets: measured +41.0° and −41.1° against +30.4° and −30.2°.
        let (jet, was) = both("jet_latitude_deg");
        assert!(jet >= was + 5.0, "{side}ern jet at {jet}° against {was}°");
    }

    // Equator-to-pole contrast: measured 15.8 K against 16.3 K
    // (16.5 against 18.5 K at 2,562 cells).
    let (contrast, was) = (
        get(&slow, "equator_pole_contrast_k"),
        get(reference, "equator_pole_contrast_k"),
    );
    assert!(
        contrast < was,
        "equator-pole contrast {contrast} K against {was} K"
    );
}

#[test]
#[ignore = "slow: two runs of nine model years at 2,562 cells, together (20–60 minutes); run with --ignored"]
fn a_shorter_day_narrows_the_hadley_cells_and_adds_a_belt_of_westerlies() {
    // A 12-hour day in a year of the same length in seconds, and the
    // reference on the same grid. 642 cells don't resolve the storms of the
    // shorter day, which are half the size.
    let (fast, reference) = std::thread::scope(|scope| {
        let fast = scope.spawn(|| {
            run_on(
                Grid::Cells2562,
                "fast",
                &["--solar-day", "12", "--year-days", "730.4844"],
            )
            .1
        });
        let reference = run_on(Grid::Cells2562, "reference-2562", &[]).1;
        (fast.join().expect("the fast run"), reference)
    });
    for side in ["north", "south"] {
        let both = |name: &str| {
            let name = format!("{name}_{side}");
            (get(&fast, &name).abs(), get(&reference, &name).abs())
        };

        // Hadley edges: measured +21.2° and −20.7° against +26.2° and −25.7°.
        let (edge, was) = both("hadley_edge_deg");
        assert!(
            edge <= was - 3.0,
            "{side}ern Hadley edge {edge}° against {was}°"
        );

        // Hadley strengths: measured +17.3e9 and −16.9e9 kg/s against
        // +35.1e9 and −40.1e9.
        let (strength, was) = both("hadley_strength_kg_s");
        assert!(
            strength <= 0.7 * was,
            "{side}ern Hadley cell {strength} against {was} kg/s"
        );

        // Upper-layer jets: measured ±25.5° against +32.3° and −32.1°.
        let (jet, was) = both("jet_latitude_deg");
        assert!(jet <= was - 3.0, "{side}ern jet at {jet}° against {was}°");

        // Belts of lowest-layer westerlies: measured two each side against
        // one.
        let (belts, was) = both("westerly_belts");
        assert!(
            belts > was,
            "{belts} {side}ern belts of westerlies against {was}"
        );
    }

    // Equator-to-pole contrast: measured 24.1 K against 18.0 K.
    let (contrast, was) = (
        get(&fast, "equator_pole_contrast_k"),
        get(&reference, "equator_pole_contrast_k"),
    );
    assert!(
        contrast >= was + 3.0,
        "equator-pole contrast {contrast} K against {was} K"
    );
}
