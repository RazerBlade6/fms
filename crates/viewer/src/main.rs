//! Desktop viewer for the fms-climate simulator.
//!
//! See `docs/climate/PLAN.md` for the design.

// Linked so the dependency graph is enforced from the start.
use fms_sim as _;

fn main() {
    println!(
        "fms-viewer {}: GUI (not implemented yet)",
        env!("CARGO_PKG_VERSION")
    );
}
