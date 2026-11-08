//! Single-column radiative–convective equilibrium: Earth calibration, CO₂
//! doubling, and sweeps over CO₂ and surface pressure.
//!
//! ```sh
//! cargo run --release -p fms-atmos --example rcm
//! ```
//!
//! Parameters can be overridden with environment variables named after the
//! fields of `RadiationParams` (e.g. `co2_width=150`), for calibration.

use fms_atmos::{Equilibrium, SingleColumn};

fn env(name: &str, default: f64) -> f64 {
    std::env::var(name)
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(default)
}

fn column() -> SingleColumn {
    let mut c = SingleColumn::earth(1361.0);
    let r = &mut c.params.radiation;
    r.co2_width = env("co2_width", r.co2_width);
    r.co2_wing = env("co2_wing", r.co2_wing);
    r.water_opacity = env("water_opacity", r.water_opacity);
    r.water_edge = env("water_edge", r.water_edge);
    r.continuum_opacity = env("continuum_opacity", r.continuum_opacity);
    r.water_absorption = env("water_absorption", r.water_absorption);
    r.ozone_absorption = env("ozone_absorption", r.ozone_absorption);
    r.albedo = env("albedo", r.albedo);
    r.cloud_albedo = env("cloud_albedo", r.cloud_albedo);
    r.surface_humidity = env("surface_humidity", r.surface_humidity);
    c.params.lapse_rate = env("lapse_rate", c.params.lapse_rate);
    c.params.transfer = env("transfer", c.params.transfer);
    c
}

fn show(label: &str, e: &Equilibrium, insolation: f64) {
    let f = &e.fluxes;
    println!(
        "{label:<22} T_s {:>6.1} K  T {:>5.1}/{:>5.1}  OLR {:>6.1}  ASR {:>6.1}  albedo {:.3}  SH {:>5.1} LH {:>5.1}  water {:>5.1} kg/m²  band {:>5.0} cm⁻¹  surface SW {:>5.1} LW down {:>5.1} up {:>5.1} net {:>5.1}  ({} h)",
        e.surface_temperature,
        e.temperature[0],
        e.temperature[1],
        f.outgoing,
        f.absorbed,
        1.0 - f.absorbed / insolation,
        e.surface_fluxes.0,
        e.surface_fluxes.1,
        f.water[0] + f.water[1],
        f.co2_width,
        f.surface_shortwave,
        f.surface_down,
        f.surface_up,
        f.surface_up - f.surface_down,
        e.steps
    );
}

fn main() {
    // Fluxes of an observed global-mean Earth column, before any equilibrium.
    let mut fixed = column();
    fixed.surface_temperature = 288.0;
    fixed.set_temperatures([226.0, 272.0]);
    let f = fixed.fluxes();
    let (sh, lh) = fixed.surface_fluxes();
    println!(
        "observed column: OLR {:.1}  ASR {:.1} (atmosphere {:.1}, surface {:.1})  albedo {:.3}  surface LW down {:.1} up {:.1}  SH {sh:.1} LH {lh:.1}  water {:.1}/{:.1} kg/m²  heating {:.1}/{:.1}",
        f.outgoing,
        f.absorbed,
        f.absorbed - f.surface_shortwave,
        f.surface_shortwave,
        1.0 - f.absorbed / fixed.insolation,
        f.surface_down,
        f.surface_up,
        f.water[0],
        f.water[1],
        f.layers[0],
        f.layers[1] + sh + lh,
    );
    let mut doubled = fixed;
    doubled.params.co2 *= 2.0;
    let mut moist = fixed;
    moist.set_temperatures([227.0, 273.0]);
    moist.surface_temperature = 289.0;
    println!(
        "  fixed-column forcing per CO₂ doubling {:.2} W/m²; OLR change for +1 K everywhere (fixed RH) {:.2} W/m²",
        f.outgoing - doubled.fluxes().outgoing,
        moist.fluxes().outgoing - f.outgoing
    );
    if std::env::var("FIXED_ONLY").is_ok() {
        return;
    }
    let base = column();
    let e = base.equilibrium();
    show("Earth (420 ppm)", &e, base.insolation);

    // Instantaneous doubling forcing at the equilibrium's temperatures.
    let mut eq_col = base;
    eq_col.surface_temperature = e.surface_temperature;
    eq_col.set_temperatures(e.temperature);
    let mut at = eq_col;
    let olr1 = at.fluxes().outgoing - at.fluxes().absorbed;
    at.params.co2 *= 2.0;
    let olr2 = at.fluxes().outgoing - at.fluxes().absorbed;
    println!(
        "doubling forcing at fixed temperatures: {:.2} W/m²",
        olr1 - olr2
    );
    // Adjusted forcing: surface held, atmosphere relaxes (two years).
    let adjusted = |co2: f64| {
        let mut c = eq_col;
        c.params.co2 = co2;
        c.params.slab_depth = 1e12;
        for _ in 0..(2 * 365 * 24) {
            c.step(3600.0);
        }
        let f = c.fluxes();
        f.absorbed - f.outgoing
    };
    println!(
        "adjusted forcing per doubling (surface fixed): {:.2} W/m²",
        adjusted(2.0 * base.params.co2) - adjusted(base.params.co2)
    );
    let mut doubled = base;
    doubled.params.co2 *= 2.0;
    let e2 = doubled.equilibrium();
    show("2× CO₂", &e2, base.insolation);
    println!(
        "equilibrium warming per doubling: {:.2} K",
        e2.surface_temperature - e.surface_temperature
    );

    println!("CO₂ sweep:");
    for co2 in [0.0, 1.0, 10.0, 100.0, 280.0, 420.0, 1000.0, 1e4, 1e5] {
        let mut c = base;
        c.params.co2 = co2;
        show(&format!("  {co2} ppm"), &c.equilibrium(), c.insolation);
    }
    println!("pressure sweep:");
    for bar in [0.1, 0.3, 1.0, 3.0, 10.0] {
        let mut c = base;
        c.ps = bar * 1.0e5;
        show(&format!("  {bar} bar"), &c.equilibrium(), c.insolation);
    }
    println!("sunlight sweep:");
    for s in [0.7, 0.85, 1.0, 1.1, 1.2] {
        let mut c = base;
        c.insolation *= s;
        show(&format!("  {s}× S"), &c.equilibrium(), c.insolation);
    }
}
