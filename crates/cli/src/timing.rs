//! The time loop's timing report: where a run's wall-clock time went, phase
//! by phase ([`fms_atmos::timing`]), and how the Helmholtz solver did in
//! each vertical mode. Printed to stderr after the run, and written as
//! `timing.csv` and `solver.csv`.

use std::fmt::Write as _;
use std::path::Path;
use std::time::Duration;

use fms_atmos::{Phase, SolveReport, Timings};

/// One vertical mode's solves in this run.
#[derive(Clone, Copy, Debug, Default)]
pub struct ModeSolves {
    pub solves: u64,
    pub iterations: u64,
    /// Most iterations any one solve took.
    pub most: u32,
    /// Worst final relative residual.
    pub worst_residual: f64,
}

impl ModeSolves {
    pub fn add(&mut self, r: &SolveReport) {
        self.solves += 1;
        self.iterations += u64::from(r.iterations);
        self.most = self.most.max(r.iterations);
        self.worst_residual = self.worst_residual.max(r.relative_residual);
    }

    fn per_solve(&self) -> f64 {
        self.iterations as f64 / self.solves.max(1) as f64
    }
}

/// Everything the report needs about one run's time loop.
pub struct LoopTiming {
    pub steps: u64,
    /// The loop's wall-clock time.
    pub elapsed: Duration,
    pub timings: Timings,
    pub threads: usize,
    pub cells: usize,
    /// Per vertical mode, external first: wave speed (m s⁻¹), solves and
    /// solve time (the modes are solved side by side, so these overlap).
    pub wave_speeds: Vec<f64>,
    pub modes: Vec<ModeSolves>,
    pub solve_time: Vec<Duration>,
}

impl LoopTiming {
    /// The loop time not covered by its top-level phases: the loop's own
    /// overhead, and the step's hand-off to the thread pool.
    pub fn unaccounted(&self) -> Duration {
        let covered: Duration = Phase::ALL
            .iter()
            .filter(|p| p.level() == 0)
            .map(|&p| self.timings.time(p))
            .sum();
        self.elapsed.saturating_sub(covered)
    }

    fn ms_per_step(&self, d: Duration) -> f64 {
        1e3 * d.as_secs_f64() / self.steps.max(1) as f64
    }

    fn share(&self, d: Duration) -> f64 {
        100.0 * d.as_secs_f64() / self.elapsed.as_secs_f64().max(f64::MIN_POSITIVE)
    }

    /// The report, as printed to stderr.
    pub fn report(&self) -> String {
        let mut out = String::new();
        let _ = writeln!(
            out,
            "timing: {} steps in {:.2} s ({:.3} ms a step), {} threads, {} cells",
            self.steps,
            self.elapsed.as_secs_f64(),
            self.ms_per_step(self.elapsed),
            self.threads,
            self.cells
        );
        let _ = writeln!(
            out,
            "  {:<44} {:>9} {:>7} {:>8} {:>9}",
            "phase", "ms/step", "share", "calls", "ms/call"
        );
        for p in Phase::ALL {
            let calls = self.timings.calls(p);
            if calls == 0 {
                continue;
            }
            let time = self.timings.time(p);
            let name = format!(
                "{}{}{}",
                "  ".repeat(usize::from(p.level())),
                if p.side_by_side() { "∥ " } else { "" },
                p.name()
            );
            let _ = writeln!(
                out,
                "  {name:<44} {:>9.3} {:>6.1}% {calls:>8} {:>9.3}",
                self.ms_per_step(time),
                self.share(time),
                1e3 * time.as_secs_f64() / calls as f64,
            );
        }
        let rest = self.unaccounted();
        let _ = writeln!(
            out,
            "  {:<44} {:>9.3} {:>6.1}%",
            "unaccounted",
            self.ms_per_step(rest),
            self.share(rest)
        );
        let _ = writeln!(
            out,
            "  (∥ runs side by side with its siblings, inside a rayon::join: their times overlap)"
        );
        let _ = writeln!(
            out,
            "solver: conjugate gradient by vertical mode (solved side by side: their times overlap)"
        );
        let _ = writeln!(
            out,
            "  {:<12} {:>11} {:>8} {:>11} {:>6} {:>10} {:>9}",
            "mode", "speed (m/s)", "solves", "iter/solve", "most", "worst res", "ms/solve"
        );
        let mut total = ModeSolves::default();
        for (k, m) in self.modes.iter().enumerate() {
            let time = self.solve_time.get(k).copied().unwrap_or_default();
            let _ = writeln!(
                out,
                "  {:<12} {:>11.1} {:>8} {:>11.1} {:>6} {:>10.1e} {:>9.3}",
                mode_name(k),
                self.wave_speeds.get(k).copied().unwrap_or(f64::NAN),
                m.solves,
                m.per_solve(),
                m.most,
                m.worst_residual,
                1e3 * time.as_secs_f64() / m.solves.max(1) as f64,
            );
            total.solves += m.solves;
            total.iterations += m.iterations;
            total.most = total.most.max(m.most);
            total.worst_residual = total.worst_residual.max(m.worst_residual);
        }
        let _ = writeln!(
            out,
            "  {:<12} {:>11} {:>8} {:>11.1} {:>6} {:>10.1e}",
            "total",
            "",
            total.solves,
            total.per_solve(),
            total.most,
            total.worst_residual,
        );
        out
    }

    /// `timing.csv` (the loop, then every phase beneath it, then what the
    /// phases leave unaccounted) and `solver.csv` (a row per vertical mode)
    /// in `dir`.
    pub fn write_csv(&self, dir: &Path) -> Result<(), String> {
        let mut csv = String::from("phase,level,side_by_side,seconds,calls\n");
        let _ = writeln!(
            csv,
            "loop,0,false,{:.9},{}",
            self.elapsed.as_secs_f64(),
            self.steps
        );
        for p in Phase::ALL {
            let _ = writeln!(
                csv,
                "{},{},{},{:.9},{}",
                p.name(),
                p.level() + 1,
                p.side_by_side(),
                self.timings.time(p).as_secs_f64(),
                self.timings.calls(p)
            );
        }
        let _ = writeln!(
            csv,
            "unaccounted,1,false,{:.9},{}",
            self.unaccounted().as_secs_f64(),
            self.steps
        );
        let path = dir.join("timing.csv");
        std::fs::write(&path, csv).map_err(|e| format!("cannot write {}: {e}", path.display()))?;
        println!("timing -> {}", path.display());

        let mut csv = String::from(
            "mode,wave_speed_m_s,solves,iterations,most_iterations,worst_residual,seconds\n",
        );
        for (k, m) in self.modes.iter().enumerate() {
            let _ = writeln!(
                csv,
                "{},{},{},{},{},{:e},{:.9}",
                mode_name(k),
                self.wave_speeds.get(k).copied().unwrap_or(f64::NAN),
                m.solves,
                m.iterations,
                m.most,
                m.worst_residual,
                self.solve_time
                    .get(k)
                    .copied()
                    .unwrap_or_default()
                    .as_secs_f64()
            );
        }
        let path = dir.join("solver.csv");
        std::fs::write(&path, csv).map_err(|e| format!("cannot write {}: {e}", path.display()))?;
        println!("solver -> {}", path.display());
        Ok(())
    }
}

/// "external", then "internal 1", "internal 2" …
fn mode_name(k: usize) -> String {
    if k == 0 {
        "external".into()
    } else {
        format!("internal {k}")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn unaccounted_is_the_loop_less_its_top_level_phases() {
        let mut timings = Timings::default();
        let ms = Duration::from_millis(1);
        timings.add(Phase::Step, 6 * ms);
        // Nested: already inside the step.
        timings.add(Phase::Radiation, 4 * ms);
        timings.add(Phase::HealthCheck, ms);
        let t = LoopTiming {
            steps: 2,
            elapsed: 10 * ms,
            timings,
            threads: 1,
            cells: 642,
            wave_speeds: vec![300.0],
            modes: vec![ModeSolves::default()],
            solve_time: vec![Duration::ZERO],
        };
        assert_eq!(t.unaccounted(), 3 * ms);
        let report = t.report();
        assert!(report.contains("radiation"));
        assert!(
            !report.contains("dynamics"),
            "phases never called are left out"
        );
    }
}
