//! Radiation (M2 step 8), single column: the same physics the 3-D model runs,
//! in radiative–convective equilibrium. Earth's global-mean climate, the
//! response to CO₂, and sanity across the planet parameters' full ranges.
//! The 3-D tests are in `climate.rs`.

use fms_atmos::{Equilibrium, SingleColumn};

/// Earth's mean stellar flux (W m⁻²).
const S0: f64 = 1361.0;

fn earth() -> SingleColumn {
    SingleColumn::earth(S0)
}

fn albedo(e: &Equilibrium, column: &SingleColumn) -> f64 {
    1.0 - e.fluxes.absorbed / column.insolation
}

#[test]
fn earths_column_is_earthlike() {
    let c = earth();
    let e = c.equilibrium();
    let f = &e.fluxes;
    let water = f.water[0] + f.water[1];
    eprintln!(
        "T_s {:.1} K, layers {:.1}/{:.1} K, OLR {:.1}, albedo {:.3}, water {water:.1} kg/m², SH {:.1} LH {:.1}",
        e.surface_temperature,
        e.temperature[0],
        e.temperature[1],
        f.outgoing,
        albedo(&e, &c),
        e.surface_fluxes.0,
        e.surface_fluxes.1
    );
    // Measured: 286.7 K, 232/274 K, 247 W/m², 0.274, 21.5 kg/m², 22 and 77 W/m².
    // Earth: 288 K, ≈226/272 K, 240 W/m², 0.30, 25 kg/m², 20 and 85 W/m².
    assert!((283.0..291.0).contains(&e.surface_temperature));
    assert!((222.0..240.0).contains(&e.temperature[0]));
    assert!((265.0..280.0).contains(&e.temperature[1]));
    assert!((235.0..255.0).contains(&f.outgoing));
    assert!((0.24..0.32).contains(&albedo(&e, &c)));
    assert!((15.0..30.0).contains(&water));
    assert!((10.0..35.0).contains(&e.surface_fluxes.0));
    assert!((60.0..100.0).contains(&e.surface_fluxes.1));
    // In equilibrium the top of the atmosphere balances.
    assert!((f.absorbed - f.outgoing).abs() < 0.01);
}

#[test]
fn doubling_co2_forces_and_warms_like_earth() {
    let c = earth();
    let e = c.equilibrium();
    // Forcing: the change in net downward flux at the top when CO₂ doubles,
    // at the equilibrium column's temperatures.
    let mut at = c;
    at.surface_temperature = e.surface_temperature;
    at.set_temperatures(e.temperature);
    let before = at.fluxes();
    at.params.co2 *= 2.0;
    let after = at.fluxes();
    let forcing = before.outgoing - after.outgoing;
    let mut doubled = c;
    doubled.params.co2 *= 2.0;
    let warming = doubled.equilibrium().surface_temperature - e.surface_temperature;
    eprintln!("doubling: forcing {forcing:.2} W/m², warming {warming:.2} K");
    // Myhre et al. (1998): 3.7 W/m² (adjusted). Measured: 3.6 W/m².
    assert!((3.2..4.2).contains(&forcing), "forcing {forcing}");
    // Measured 2.4 K: clear sky, fixed relative humidity, no ice or cloud
    // feedbacks yet (Earth: ~3 K with them).
    assert!((2.0..4.5).contains(&warming), "warming {warming}");
}

/// Equilibrium surface temperatures across a sweep of one parameter.
fn sweep(values: &[f64], set: impl Fn(&mut SingleColumn, f64)) -> Vec<f64> {
    values
        .iter()
        .map(|&v| {
            let mut c = earth();
            set(&mut c, v);
            let e = c.equilibrium();
            let f = &e.fluxes;
            assert!(
                e.surface_temperature.is_finite() && (f.absorbed - f.outgoing).abs() < 0.01,
                "no equilibrium at {v}: {e:?}"
            );
            e.surface_temperature
        })
        .collect()
}

fn increasing(name: &str, t: &[f64]) {
    eprintln!("{name}: {t:.1?}");
    assert!(
        t.windows(2).all(|w| w[1] > w[0]),
        "{name} not increasing: {t:?}"
    );
}

#[test]
fn co2_across_its_whole_range() {
    // The planet parameter's range: 0–100,000 ppm.
    let t = sweep(&[0.0, 1.0, 10.0, 100.0, 420.0, 1e3, 1e4, 1e5], |c, v| {
        c.params.co2 = v
    });
    increasing("CO₂", &t);
    // Without CO₂ the planet is still kept warm by water vapour (measured
    // 278 K); 10% CO₂ warms it by about 20 K, not a runaway.
    assert!(t[0] > 270.0 && t[7] < t[4] + 30.0);
}

#[test]
fn pressure_across_its_whole_range() {
    // 0.1–10 bar: thicker air holds more greenhouse gas (and broadens its
    // lines), so the surface is warmer.
    let t = sweep(&[1e4, 3e4, 1e5, 3e5, 1e6], |c, v| c.ps = v);
    increasing("pressure", &t);
    assert!(t[4] < 330.0);
}

#[test]
fn sunlight_and_gravity() {
    let t = sweep(&[0.7, 0.85, 1.0, 1.1], |c, v| c.insolation *= v);
    increasing("sunlight", &t);
    // Weaker gravity: more air mass per pascal, so more CO₂ and water in the
    // column at the same surface pressure.
    let t = sweep(&[20.0, 9.81, 3.0], |c, v| c.gravity = v);
    increasing("1/gravity", &t);
}
