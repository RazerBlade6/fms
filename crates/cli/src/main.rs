//! Headless command-line runner for the Fantasy Map Simulator.
//!
//! See `docs/PLAN.md` for the design.

// Linked so the dependency graph is enforced from the start.
use fms_sim as _;

fn main() {
    println!(
        "fms {}: headless runner (not implemented yet)",
        env!("CARGO_PKG_VERSION")
    );
}
