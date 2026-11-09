//! Single-column cases for the boundary layer and shallow convection: the
//! column physics as a per-column function, forced by a case's prescribed
//! large-scale tendencies, surface fluxes and radiation, and run for hours in
//! milliseconds.
//!
//! - [`bomex`]: trade cumulus over the tropical Atlantic (Siebesma et al.
//!   2003).
//! - [`dycoms`]: nocturnal marine stratocumulus off California, DYCOMS-II
//!   RF01 (Stevens et al. 2005).
//!
//! The cases' setups are written from memory of those papers, which aren't in
//! `docs/sources`; check them against the papers before trusting a close
//! comparison. See `docs/OBSERVATIONS.md` O-35.
//!
//! These cases resolve boundary-layer structure the model's few layers can't:
//! a case's bottom layer holds a whole cloud-topped boundary layer, so only
//! bulk quantities can be compared.

use crate::column_physics::{BoundaryColumn, PhysicsParams, boundary_layer_top, theta_tendency};
use crate::layered::Gas;
use crate::levels::{Levels, MAX_LAYERS};
use crate::moisture::{LATENT_HEAT, MoistColumn, Setting, moist_physics, stratiform_cover};
use crate::radiation::saturation_humidity;
use crate::williamson::GRAVITY;

/// Sub-intervals of a layer in its height averages.
const LAYER_SAMPLES: usize = 32;
/// The case's time step (s).
const STEP: f64 = 300.0;
const DAY: f64 = 86_400.0;

/// A column case.
#[derive(Clone, Copy, Debug)]
pub struct ColumnCase {
    /// Name, for messages.
    pub name: &'static str,
    /// Surface pressure (Pa).
    pub ps: f64,
    /// Length of the run (hours); the last half is averaged.
    pub hours: f64,
    /// Temperature (K) for converting between σ and height.
    pub mean_temperature: f64,
    /// Interior layer interfaces (σ).
    pub interior: [f64; 3],
    /// Surface sensible and latent heat flux (W m⁻²).
    pub sensible: f64,
    /// See `sensible`.
    pub latent: f64,
    /// Longwave cooling of the bottom layer by cloud at its top (W m⁻²).
    pub cloud_top_cooling: f64,
    /// Initial potential temperature (K) with height (m).
    pub theta: fn(f64) -> f64,
    /// Initial specific humidity (kg kg⁻¹) with height.
    pub humidity: fn(f64) -> f64,
    /// Prescribed cooling (K s⁻¹; negative cools) with height.
    pub cooling: fn(f64) -> f64,
    /// Prescribed drying (s⁻¹; negative dries) with height.
    pub drying: fn(f64) -> f64,
    /// Large-scale vertical velocity (m s⁻¹; negative sinks) with height.
    pub subsidence: fn(f64) -> f64,
}

/// Results averaged over the last half of a case, and the bottom layer's
/// change over the whole run.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct CaseResult {
    /// Change of the bottom layer's potential temperature (K).
    pub theta_change: f64,
    /// Change of the bottom layer's specific humidity (kg kg⁻¹).
    pub humidity_change: f64,
    /// Bottom layer's relative humidity.
    pub relative_humidity: f64,
    /// Rain reaching the surface (mm day⁻¹).
    pub rain: f64,
    /// Convective part of the rain (mm day⁻¹).
    pub convective_rain: f64,
    /// Stratiform cloud cover of the bottom layer.
    pub cloud: f64,
    /// Convective velocity `w*` (m s⁻¹).
    pub convective_velocity: f64,
    /// Entrainment velocity (m s⁻¹).
    pub entrainment_velocity: f64,
    /// Share of time with shallow cumulus at the top of the boundary layer.
    pub cumulus: f64,
    /// Column water at the end minus what its sources and sinks account for
    /// (kg m⁻²).
    pub water_residual: f64,
}

/// Run `case` with `params`.
pub fn run_case(case: &ColumnCase, params: &PhysicsParams) -> CaseResult {
    let gas = Gas::EARTH_AIR;
    let g = GRAVITY;
    let levels = Levels::new(&case.interior);
    let n = levels.count();
    let bottom = n - 1;
    let above = bottom - 1;
    let ds = levels.dsigma();
    let ps = case.ps;
    let scale_height = gas.gas_constant * case.mean_temperature / g;
    let layer_mean = |k: usize, f: &dyn Fn(f64) -> f64| {
        let (top, base) = levels.bounds(k);
        (0..LAYER_SAMPLES)
            .map(|j| {
                let sigma = top + (base - top) * (j as f64 + 0.5) / LAYER_SAMPLES as f64;
                f(-scale_height * fms_math::ln(sigma))
            })
            .sum::<f64>()
            / LAYER_SAMPLES as f64
    };
    let advection = |profile: fn(f64) -> f64, z: f64| {
        let slope = 0.5 * (profile(z + 1.0) - profile(z - 1.0));
        -(case.subsidence)(z) * slope
    };
    let theta_forcing: Vec<f64> = (0..n)
        .map(|k| layer_mean(k, &|z| (case.cooling)(z) + advection(case.theta, z)))
        .collect();
    let humidity_forcing: Vec<f64> = (0..n)
        .map(|k| layer_mean(k, &|z| (case.drying)(z) + advection(case.humidity, z)))
        .collect();
    let mut theta: Vec<f64> = (0..n).map(|k| layer_mean(k, &case.theta)).collect();
    let mut humidity: Vec<f64> = (0..n).map(|k| layer_mean(k, &case.humidity)).collect();
    let (theta_start, humidity_start) = (theta[bottom], humidity[bottom]);

    let setting = Setting {
        params: &params.moisture,
        radiation: &params.radiation,
        gas,
        levels: &levels,
        gravity: g,
    };
    let sl = &params.surface_layer;
    let cp = gas.heat_capacity;
    let rd = gas.gas_constant;
    let top = levels.interfaces()[bottom];
    let evaporation = case.latent / LATENT_HEAT;
    let surface_exner = gas.exner(ps);
    let mut exner = [0.0; MAX_LAYERS];
    levels.layer_exner_into(gas, ps, &mut exner[..n]);
    let exner_top = gas.exner(top * ps);

    let masses: Vec<f64> = ds.iter().map(|d| d * ps / g).collect();
    let layer_pressure = levels.mean_pressures(gas, ps)[bottom];
    let water = |q: &[f64]| (0..n).map(|k| q[k] * masses[k]).sum::<f64>();
    let start_water = water(&humidity);
    let mut expected = 0.0;
    let steps = (case.hours * 3600.0 / STEP).round() as usize;
    let averaged_from = steps / 2;
    let mut result = CaseResult::default();
    for step in 0..steps {
        let temperature: Vec<f64> = (0..n).map(|k| theta[k] * exner[k]).collect();
        let air = theta[bottom] * surface_exner;
        let moist = moist_physics(
            &setting,
            &MoistColumn {
                ps,
                temperature: temperature.clone(),
                surface_humidity: humidity[bottom].clamp(0.0, saturation_humidity(air, ps)),
                humidity: humidity.clone(),
                surface_air: air,
                stratosphere: 210.0,
                orographic_ascent: 0.0,
            },
        );
        let t = temperature[bottom];
        let depth = -(rd * t / g) * fms_math::ln(top);
        let layer_top = boundary_layer_top(
            &BoundaryColumn {
                gas,
                gravity: g,
                lapse_rate: params.radiation.lapse_rate,
                top,
                depth,
                ps,
                theta: theta[bottom],
                humidity: humidity[bottom].max(0.0),
                layer_humidity: humidity[bottom],
                air,
                density: ps / (rd * air),
                theta_above: theta[above],
                humidity_above: humidity[above],
                exner_above: exner[above],
                sensible: case.sensible,
                evaporation,
                cloud_top_cooling: case.cloud_top_cooling,
            },
            sl,
        );
        let entrainment = layer_top.density * layer_top.entrainment_velocity;
        let cumulus = if layer_top.cumulus {
            sl.cumulus * layer_top.density * layer_top.convective_velocity
        } else {
            0.0
        };
        let heat = entrainment * cp * exner_top * (theta[above] - theta[bottom]);
        let moisture = g * (entrainment + cumulus) * (humidity[above] - humidity[bottom]);

        expected += STEP
            * (evaporation - moist.rain - moist.snow
                + (0..n).map(|k| humidity_forcing[k] * masses[k]).sum::<f64>());
        if step >= averaged_from {
            let rh = humidity[bottom] / saturation_humidity(temperature[bottom], layer_pressure);
            result.relative_humidity += rh;
            result.rain += (moist.rain + moist.snow) * DAY;
            result.convective_rain += moist.convective * DAY;
            let clouds = &params.moisture.clouds;
            result.cloud +=
                stratiform_cover(rh, clouds.critical_humidity[1], clouds.overcast_humidity[1]);
            result.convective_velocity += layer_top.convective_velocity;
            result.entrainment_velocity += layer_top.entrainment_velocity;
            result.cumulus += f64::from(u8::from(layer_top.cumulus));
        }

        for k in 0..n {
            let mut heating = 0.0;
            let mut water = ps * moist.humidity[k];
            if k == bottom {
                heating += case.sensible - case.cloud_top_cooling + heat;
                water += (g * evaporation + moisture) / ds[k];
            } else if k == above {
                heating -= heat;
                water -= moisture / ds[k];
            }
            let d_theta = theta_tendency(gas, g, ds[k], exner[k], heating)
                + ps * moist.temperature[k] / exner[k]
                + ps * theta_forcing[k];
            theta[k] += STEP * d_theta / ps;
            humidity[k] += STEP * (water + ps * humidity_forcing[k]) / ps;
        }
    }
    let samples = (steps - averaged_from) as f64;
    for x in [
        &mut result.relative_humidity,
        &mut result.rain,
        &mut result.convective_rain,
        &mut result.cloud,
        &mut result.convective_velocity,
        &mut result.entrainment_velocity,
        &mut result.cumulus,
    ] {
        *x /= samples;
    }
    result.theta_change = theta[bottom] - theta_start;
    result.humidity_change = humidity[bottom] - humidity_start;
    result.water_residual = water(&humidity) - start_water - expected;
    result
}

/// Linear from `a` at `z0` to `b` at `z1`, constant outside.
fn ramp(z: f64, z0: f64, z1: f64, a: f64, b: f64) -> f64 {
    a + (b - a) * ((z - z0) / (z1 - z0)).clamp(0.0, 1.0)
}

/// Piecewise-linear through `(z, value)` points; beyond the last, `above`.
fn through(z: f64, points: &[(f64, f64)], above: impl Fn(f64, f64) -> f64) -> f64 {
    let &(top, last) = points.last().unwrap();
    if z >= top {
        return above(z - top, last);
    }
    for pair in points.windows(2) {
        let ((z0, a), (z1, b)) = (pair[0], pair[1]);
        if z < z1 {
            return ramp(z, z0, z1, a, b);
        }
    }
    last
}

/// Trade cumulus: BOMEX (Siebesma et al. 2003). The forcing was designed
/// to hold the observed state steady, with cloud base near 500 m, top near
/// 1500 m and cover of about a tenth.
pub fn bomex() -> ColumnCase {
    const DENSITY: f64 = 1.16;
    ColumnCase {
        name: "BOMEX",
        ps: 101_500.0,
        hours: 6.0,
        mean_temperature: 295.0,
        interior: [0.25, 0.5, 0.85],
        sensible: DENSITY * Gas::EARTH_AIR.heat_capacity * 8e-3,
        latent: DENSITY * LATENT_HEAT * 5.2e-5,
        cloud_top_cooling: 0.0,
        theta: |z| {
            through(
                z,
                &[
                    (0.0, 298.7),
                    (520.0, 298.7),
                    (1480.0, 302.4),
                    (2000.0, 308.2),
                    (3000.0, 311.85),
                ],
                |dz, top| top + 3.65e-3 * dz,
            )
        },
        humidity: |z| {
            1e-3 * through(
                z,
                &[
                    (0.0, 17.0),
                    (520.0, 16.3),
                    (1480.0, 10.7),
                    (2000.0, 4.2),
                    (3000.0, 3.0),
                ],
                |dz, top| top * fms_math::exp(-dz / 2000.0),
            )
        },
        cooling: |z| ramp(z, 1500.0, 2500.0, -2.0 / DAY, 0.0),
        drying: |z| ramp(z, 300.0, 500.0, -1.2e-8, 0.0),
        subsidence: |z| {
            if z < 1500.0 {
                -6.5e-3 * z / 1500.0
            } else {
                ramp(z, 1500.0, 2100.0, -6.5e-3, 0.0)
            }
        },
    }
}

/// Inversion height of [`dycoms`] (m).
const DYCOMS_INVERSION: f64 = 840.0;

/// Nocturnal stratocumulus: DYCOMS-II RF01 (Stevens et al. 2005). The bottom
/// layer's top sits at the inversion, 840 m, instead of the climate's σ
/// 0.85, or the layer would average the free troposphere into the cloud.
pub fn dycoms() -> ColumnCase {
    let mean_temperature = 289.0;
    let scale = Gas::EARTH_AIR.gas_constant * mean_temperature / GRAVITY;
    ColumnCase {
        name: "DYCOMS-II RF01",
        ps: 101_780.0,
        hours: 4.0,
        mean_temperature,
        interior: [0.25, 0.5, fms_math::exp(-DYCOMS_INVERSION / scale)],
        sensible: 15.0,
        latent: 115.0,
        cloud_top_cooling: 70.0,
        theta: |z| {
            if z < DYCOMS_INVERSION {
                289.0
            } else {
                297.5 + fms_math::pow(z - DYCOMS_INVERSION, 1.0 / 3.0)
            }
        },
        humidity: |z| {
            if z < DYCOMS_INVERSION {
                9e-3
            } else {
                1.5e-3 * fms_math::exp(-(z - DYCOMS_INVERSION) / 3000.0)
            }
        },
        cooling: |_| 0.0,
        drying: |_| 0.0,
        subsidence: |z| -3.75e-6 * z,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::surface_layer::{SurfaceLayerParams, entrainment_velocity};

    /// [`PhysicsParams`] with deep convection switched off.
    fn without_convection() -> PhysicsParams {
        let mut params = PhysicsParams::default();
        params.moisture.convection_time = 1e30;
        params
    }

    #[test]
    fn a_case_conserves_its_column_water() {
        for case in [bomex(), dycoms()] {
            let r = run_case(&case, &PhysicsParams::default());
            assert!(r.water_residual.abs() < 1e-10, "{}: {r:?}", case.name);
        }
    }

    #[test]
    fn the_boundary_layer_holds_bomex_steady_without_deep_convection() {
        // BOMEX's forcing was designed to hold the trade-wind state steady
        // over the run, with no rain: measured −0.004 K, +0.29 g/kg and
        // 0.0006 mm/day. Its turbulence, `w*` 0.92 m/s, is near the LES's.
        let r = run_case(&bomex(), &without_convection());
        assert!(r.theta_change.abs() < 0.5, "{r:?}");
        assert!(r.humidity_change.abs() < 1e-3, "{r:?}");
        assert!(r.rain < 0.1, "{r:?}");
        assert!((0.5..1.5).contains(&r.convective_velocity), "{r:?}");
    }

    #[test]
    fn entrainment_gets_stratocumulus_right_given_the_real_inversion() {
        // DYCOMS-II's observed entrainment is about 0.4 cm/s. The closure,
        // with the case's own `w*` and its 8.5 K jump at the inversion, gives
        // 0.19 cm/s: the right order, so the closure isn't what's wrong.
        let r = run_case(&dycoms(), &PhysicsParams::default());
        let params = SurfaceLayerParams::default();
        let w_e = entrainment_velocity(&params, GRAVITY, 800.0, r.convective_velocity, 290.0, 8.5);
        assert!((1e-3..8e-3).contains(&w_e), "w_e {w_e}");
    }

    #[test]
    #[ignore = "known gap (O-35): deep convection rains under the trade inversion"]
    fn bomex_trade_cumulus_do_not_rain() {
        let r = run_case(&bomex(), &PhysicsParams::default());
        assert!(r.rain < 0.1, "{r:?}");
        assert!(r.cumulus > 0.5, "{r:?}");
        assert!(r.cloud < 0.3, "{r:?}");
    }

    #[test]
    #[ignore = "known gap (O-35): the model's layers can't hold an inversion"]
    fn dycoms_stratocumulus_entrains_slowly_and_stays_overcast() {
        let r = run_case(&dycoms(), &PhysicsParams::default());
        assert!((2e-3..8e-3).contains(&r.entrainment_velocity), "{r:?}");
        assert!(r.cloud > 0.8, "{r:?}");
    }
}
