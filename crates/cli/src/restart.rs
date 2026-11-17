//! `fms-climate run --save-restart` and `--restart`: the run's final
//! state saved to a file, and a run started from one.
//!
//! A restarted run is built from its command line exactly as a fresh one,
//! then the file is loaded over the fresh state. The file never supplies
//! settings; the CLI adds a fingerprint of its run's (`run.settings`) so a
//! restart can say whether it continues the saved run exactly.

use std::fmt::Debug;
use std::path::Path;

use fms_atmos::{HyperdiffusionOptions, LayeredSemiImplicit, PhysicsParams, Restart};
use fms_planet::PlanetParams;

/// The CLI's array: the settings' fingerprint, its high and low 32 bits.
const SETTINGS: &str = "run.settings";

/// A 64-bit FNV-1a hash.
struct Fnv(u64);

impl Fnv {
    fn new() -> Self {
        Self(0xcbf2_9ce4_8422_2325)
    }

    fn bytes(&mut self, bytes: &[u8]) {
        for &b in bytes {
            self.0 = (self.0 ^ u64::from(b)).wrapping_mul(0x0000_0100_0000_01b3);
        }
    }

    /// One part, its length first so that parts can't run into each other.
    fn part(&mut self, bytes: &[u8]) {
        self.bytes(&(bytes.len() as u64).to_le_bytes());
        self.bytes(bytes);
    }

    fn text(&mut self, value: &impl Debug) {
        self.part(format!("{value:?}").as_bytes());
    }

    fn values(&mut self, values: &[f64]) {
        let bytes: Vec<u8> = values.iter().flat_map(|v| v.to_le_bytes()).collect();
        self.part(&bytes);
    }
}

/// What shapes a climate run's computation, to fingerprint.
#[derive(Clone, Copy)]
pub struct Settings<'a> {
    pub planet: &'a PlanetParams,
    /// The physics as the run uses them (after the step and the ocean's
    /// interval are set).
    pub physics: &'a PhysicsParams,
    /// The σ interfaces.
    pub sigma: &'a [f64],
    pub dt: f64,
    /// The hyperdiffusion's options, or `None` when inviscid.
    pub diffusion: Option<HyperdiffusionOptions>,
    /// The terrain the model uses: the smoothed surface geopotential, the
    /// relief's spread, the surface height and the land fraction.
    pub terrain: [&'a [f64]; 4],
}

impl Settings<'_> {
    /// The 64-bit FNV-1a fingerprint of the settings: their `Debug` text,
    /// and the bits of the numbers and arrays.
    pub fn fingerprint(&self) -> u64 {
        let mut h = Fnv::new();
        h.text(self.planet);
        h.text(self.physics);
        h.values(self.sigma);
        h.values(&[self.dt]);
        h.text(&self.diffusion);
        for values in self.terrain {
            h.values(values);
        }
        h.0
    }
}

/// Write the stepper's state, with the `settings` fingerprint, to `path`,
/// creating its directory if needed.
pub fn save(path: &Path, si: &LayeredSemiImplicit, settings: u64) -> Result<(), String> {
    if let Some(dir) = path.parent().filter(|d| !d.as_os_str().is_empty()) {
        std::fs::create_dir_all(dir)
            .map_err(|e| format!("cannot create {}: {e}", dir.display()))?;
    }
    let mut restart = si.save_restart();
    restart.put(
        SETTINGS,
        &[(settings >> 32) as f64, (settings & 0xffff_ffff) as f64],
    );
    restart
        .save(path)
        .map_err(|e| format!("cannot write the restart {}: {e}", path.display()))?;
    println!("restart -> {}", path.display());
    Ok(())
}

/// Load the state in `path` into the freshly built stepper, whose run has
/// the `settings` fingerprint. Returns how the state was adapted to this
/// run, its settings first; none means it continues the saved run exactly.
pub fn load(
    path: &Path,
    si: &mut LayeredSemiImplicit,
    settings: u64,
) -> Result<Vec<String>, String> {
    let refused =
        |e: &dyn std::fmt::Display| format!("cannot restart from {}: {e}", path.display());
    let restart = Restart::load(path).map_err(|e| refused(&e))?;
    let loaded = si.load_restart(&restart).map_err(|e| refused(&e))?;
    let mut notes = Vec::new();
    match restart.get_array::<2>(SETTINGS) {
        Ok([high, low]) => {
            let same = high.to_bits() == ((settings >> 32) as f64).to_bits()
                && low.to_bits() == ((settings & 0xffff_ffff) as f64).to_bits();
            if !same {
                notes.push("settings differ from the saved run's".to_owned());
            }
        }
        Err(_) => notes.push(
            "the file has no record of its run's settings (it was not saved by fms-climate), \
             so they could not be compared"
                .to_owned(),
        ),
    }
    notes.extend(loaded);
    Ok(notes)
}

/// What a restart from `path` is, for stderr and the run's description:
/// an exact continuation, or a warm start with each of its `notes` on a
/// line of its own. The stepper has just loaded it; a year is `per_year`
/// steps and a day `day` seconds.
pub fn report(
    path: &Path,
    si: &LayeredSemiImplicit,
    notes: &[String],
    per_year: u64,
    day: f64,
) -> String {
    let step = si.steps();
    let when = format!(
        "step {step} (day {:.1}, in year {})",
        step as f64 * si.dt() / day,
        step / per_year.max(1) + 1
    );
    if notes.is_empty() {
        return format!(
            "restart: an exact continuation of {} from {when}",
            path.display()
        );
    }
    let mut report = format!(
        "restart: a warm start from {} at {when}, not an exact continuation:",
        path.display()
    );
    for note in notes {
        report.push_str("\n  - ");
        report.push_str(note);
    }
    report
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fnv_1a_matches_its_reference_values() {
        let hash = |s: &str| {
            let mut h = Fnv::new();
            h.bytes(s.as_bytes());
            h.0
        };
        assert_eq!(hash(""), 0xcbf2_9ce4_8422_2325);
        assert_eq!(hash("a"), 0xaf63_dc4c_8601_ec8c);
        assert_eq!(hash("foobar"), 0x8594_4171_f739_67e8);
    }

    #[test]
    fn the_fingerprint_sees_every_setting() {
        let planet = PlanetParams::earth();
        let physics = PhysicsParams::default();
        let terrain = [vec![0.0; 4], vec![1.0; 4], vec![2.0; 4], vec![0.5; 4]];
        let base = Settings {
            planet: &planet,
            physics: &physics,
            sigma: &[0.0, 0.5, 1.0],
            dt: 1800.0,
            diffusion: None,
            terrain: terrain.each_ref().map(Vec::as_slice),
        };
        let f = base.fingerprint();
        assert_eq!(f, base.fingerprint());
        let co2 = PhysicsParams {
            co2: 840.0,
            ..physics
        };
        let mut land = terrain.clone();
        land[3][2] = 0.6;
        let others = [
            Settings {
                physics: &co2,
                ..base
            },
            Settings { dt: 1200.0, ..base },
            Settings {
                sigma: &[0.0, 0.4, 1.0],
                ..base
            },
            Settings {
                diffusion: Some(HyperdiffusionOptions::default()),
                ..base
            },
            Settings {
                terrain: land.each_ref().map(Vec::as_slice),
                ..base
            },
        ];
        for other in &others {
            assert_ne!(other.fingerprint(), f);
        }
    }
}
