//! Global diagnostics of a dynamics run, sampled through it: the mass,
//! potential-temperature and energy budgets, and the state's extremes. Written
//! as a CSV after every run, and quoted in the report when a run blows up.

use std::fmt::Write as _;
use std::path::Path;

use fms_atmos::{EnergyBudget, Health, TwoLayerDiagnostics};

/// One sample.
#[derive(Clone, Copy, Debug)]
pub struct Record {
    pub step: u64,
    pub day: f64,
    pub diagnostics: TwoLayerDiagnostics,
    /// Energy added by forcing and diffusion since the start of the run (J).
    pub budget: EnergyBudget,
    pub health: Health,
}

/// All samples of a run, with the first as the reference for changes.
#[derive(Default)]
pub struct RunLog {
    records: Vec<Record>,
}

impl RunLog {
    pub fn push(&mut self, record: Record) {
        self.records.push(record);
    }

    /// Relative changes since the first sample: mass, θ, and the
    /// energy-budget residual `(E − E₀ − forcing − diffusion) / E₀` (the
    /// time-stepping error).
    pub fn drift(&self, r: &Record) -> (f64, f64, f64) {
        let Some(s) = self.records.first() else {
            return (0.0, 0.0, 0.0);
        };
        let (a, b) = (&s.diagnostics, &r.diagnostics);
        let added = |x: &Record| x.budget.forcing.total() + x.budget.diffusion.total();
        (
            b.mass / a.mass - 1.0,
            b.theta / a.theta - 1.0,
            (b.energy - a.energy - (added(r) - added(s))) / a.energy,
        )
    }

    pub fn first(&self) -> Option<&Record> {
        self.records.first()
    }

    pub fn last(&self) -> Option<&Record> {
        self.records.last()
    }

    /// Write every sample as CSV.
    pub fn write_csv(&self, path: &Path) -> Result<(), String> {
        let mut out = String::from(
            "step,day,mass_kg,theta_kg_k,kinetic_energy_j,energy_j,\
             forcing_kinetic_j,forcing_thermal_j,diffusion_kinetic_j,diffusion_thermal_j,\
             mass_drift,theta_drift,energy_residual,max_wind_m_s,ps_min_pa,ps_max_pa,\
             theta_min_k,theta_max_k\n",
        );
        for r in &self.records {
            let (dm, dt, de) = self.drift(r);
            let (d, h, b) = (&r.diagnostics, &r.health, &r.budget);
            let _ = writeln!(
                out,
                "{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}",
                r.step,
                r.day,
                d.mass,
                d.theta,
                d.kinetic_energy,
                d.energy,
                b.forcing.kinetic,
                b.forcing.thermal,
                b.diffusion.kinetic,
                b.diffusion.thermal,
                dm,
                dt,
                de,
                h.max_wind,
                h.surface_pressure.0,
                h.surface_pressure.1,
                h.theta.0,
                h.theta.1
            );
        }
        std::fs::write(path, out).map_err(|e| format!("cannot write {}: {e}", path.display()))
    }

    /// The last `n` samples as an aligned text table.
    pub fn tail_table(&self, n: usize) -> String {
        let mut out = format!(
            "{:>7} {:>10} {:>11} {:>11} {:>11} {:>9} {:>15} {:>13}\n",
            "day", "KE (J)", "mass drift", "θ drift", "E residual", "max |u|", "p_s (hPa)", "θ (K)"
        );
        let start = self.records.len().saturating_sub(n);
        for r in &self.records[start..] {
            let (dm, dt, de) = self.drift(r);
            let h = &r.health;
            let _ = writeln!(
                out,
                "{:>7.2} {:>10.3e} {:>+11.1e} {:>+11.1e} {:>+11.1e} {:>9.1} {:>7.0}–{:<7.0} {:>5.0}–{:<5.0}",
                r.day,
                r.diagnostics.kinetic_energy,
                dm,
                dt,
                de,
                h.max_wind,
                h.surface_pressure.0 / 100.0,
                h.surface_pressure.1 / 100.0,
                h.theta.0,
                h.theta.1,
            );
        }
        out
    }
}
