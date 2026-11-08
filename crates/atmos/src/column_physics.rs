//! Column physics for a radiatively driven planet: [`crate::radiation`], a
//! temporary slab surface with bulk sensible and latent heat fluxes, and
//! convective adjustment. Used by [`ColumnPhysics`] (the 3-D forcing) and
//! [`SingleColumn`] (a 1-D radiative–convective model for calibration and
//! quick experiments), so both run exactly the same physics.
//!
//! # Pieces
//!
//! - **Surface**: one slab of water of fixed depth everywhere (an
//!   aquaplanet; land and ocean come in M2 steps 10–11). It absorbs sunlight
//!   and longwave, emits as a black body, and loses heat to the lower layer
//!   by bulk formulas: sensible `ρ c_p C V (T_s − T_a)` and latent
//!   `ρ L C V (q*(T_s) − RH_s q*(T_a))`, where `T_a` is the lower layer's
//!   temperature extrapolated to the surface along the convective lapse rate
//!   and `V` includes a gust speed.
//! - **Latent heat**: with relative humidity fixed, evaporated water is
//!   condensed at once. Its latent heat warms the whole column by the same
//!   temperature change per unit mass, as deep convection does, so the
//!   lapse rate is kept. Moisture becomes a model variable in M2 step 9.
//! - **Convection**: when the layers are less stable than a critical lapse
//!   rate (Manabe & Strickler 1964: the moist adiabat's effect until the
//!   model has moisture), both layers relax towards it, conserving
//!   enthalpy. Measured between the layer means, Earth's troposphere sits at
//!   4.4 K/km (see [`PhysicsParams::lapse_rate`]).
//! - **Friction**: the Held–Suarez boundary-layer drag on the lower layer,
//!   until M2 step 10 adds surface drag.
//!
//! # Time stepping
//!
//! Radiation runs every [`PhysicsParams::radiation_interval`] with sunlight
//! averaged over the coming interval; surface fluxes and the slab update
//! every step ([`Forcing::advance`]). The atmosphere and the slab use the same
//! stored fluxes, so the planet's energy budget closes exactly: what the
//! atmosphere and slab gain is what was absorbed from the Sun minus what was
//! emitted to space (and the kinetic energy friction and diffusion remove).

use fms_geo::Operators;
use fms_geo::operators::MIN_CHUNK;
use fms_geo::zonal::east_north;
use fms_planet::Planet;
use rayon::prelude::*;

use crate::forcing::{CellField, Forcing, ForcingBudget};
use crate::radiation::{
    Column, Fluxes, RadiationParams, column_fluxes, saturation_pressure, specific_humidity,
};
use crate::two_layer::{Gas, REFERENCE_PRESSURE, TwoLayer, TwoLayerState};

const DAY: f64 = 86_400.0;
/// Latent heat of vaporization (J kg⁻¹).
pub const LATENT_HEAT: f64 = 2.5e6;
/// Volumetric heat capacity of water (J m⁻³ K⁻¹).
const WATER_HEAT_CAPACITY: f64 = 1000.0 * 4186.0;

/// Column-physics parameters. [`Default`] is Earth-like.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PhysicsParams {
    /// Radiation.
    pub radiation: RadiationParams,
    /// CO₂ (ppm).
    pub co2: f64,
    /// Slab depth (m of water).
    pub slab_depth: f64,
    /// Bulk transfer coefficient for heat and moisture.
    pub transfer: f64,
    /// Gust speed added in quadrature to the wind (m s⁻¹).
    pub gust: f64,
    /// Critical lapse rate for convection between the layer means (K m⁻¹).
    /// Earth's observed layer means (272 and 226 K, 10.5 km apart) sit at
    /// 4.4 K/km: less than the tropospheric 6.5 K/km because the upper layer
    /// includes the stratosphere.
    pub lapse_rate: f64,
    /// Convective relaxation time (s).
    pub convection_time: f64,
    /// Lower-layer Rayleigh friction (s⁻¹).
    pub friction: f64,
    /// Time between radiation calls (s).
    pub radiation_interval: f64,
}

impl Default for PhysicsParams {
    fn default() -> Self {
        Self {
            radiation: RadiationParams::default(),
            co2: 420.0,
            slab_depth: 10.0,
            transfer: 1.3e-3,
            gust: 3.0,
            lapse_rate: 4.4e-3,
            convection_time: 2.0 * 3600.0,
            friction: 0.3 / DAY,
            radiation_interval: 3.0 * 3600.0,
        }
    }
}

impl PhysicsParams {
    /// Slab heat capacity (J m⁻² K⁻¹).
    pub fn heat_capacity(&self) -> f64 {
        WATER_HEAT_CAPACITY * self.slab_depth
    }
}

/// Layer-mean Exner functions `(p/p₀)^κ` for surface pressure `ps`.
pub fn layer_exner(gas: Gas, sigma: f64, ps: f64) -> [f64; 2] {
    let kappa = gas.kappa();
    let pi = sigma * ps;
    let (ei, es) = (gas.exner(pi), gas.exner(ps));
    [
        ei / (1.0 + kappa),
        (ps * es - pi * ei) / ((1.0 + kappa) * (1.0 - sigma) * ps),
    ]
}

/// Air temperature at the surface (K): the lower layer's mean temperature
/// extrapolated down the convective lapse rate.
pub fn surface_air_temperature(
    params: &PhysicsParams,
    gas: Gas,
    sigma: f64,
    gravity: f64,
    ps: f64,
    theta_lower: f64,
) -> f64 {
    let exner = layer_exner(gas, sigma, ps)[1];
    let p_mean = REFERENCE_PRESSURE * fms_math::pow(exner, 1.0 / gas.kappa());
    let beta = gas.gas_constant * params.lapse_rate / gravity;
    theta_lower * exner * fms_math::pow(ps / p_mean, beta)
}

/// Sensible and latent heat flux from the surface (W m⁻², positive upward)
/// for surface air at temperature `air`.
pub fn surface_fluxes(
    params: &PhysicsParams,
    gas: Gas,
    ps: f64,
    surface_temperature: f64,
    air: f64,
    wind: f64,
) -> (f64, f64) {
    let density = ps / (gas.gas_constant * air);
    let speed = (wind * wind + params.gust * params.gust).sqrt();
    let k = density * params.transfer * speed;
    let sensible = k * gas.heat_capacity * (surface_temperature - air);
    let q_surface = specific_humidity(saturation_pressure(surface_temperature), ps);
    let q_air = specific_humidity(
        params.radiation.surface_humidity * saturation_pressure(air),
        ps,
    );
    (sensible, k * LATENT_HEAT * (q_surface - q_air))
}

/// Potential temperatures the column relaxes to if it is less stable than
/// the critical lapse rate, else `None`. The adjusted column keeps its
/// enthalpy `Σ δσₖ Πₖ θₖ`.
pub fn convective_target(
    params: &PhysicsParams,
    gas: Gas,
    sigma: f64,
    gravity: f64,
    ps: f64,
    theta: [f64; 2],
) -> Option<[f64; 2]> {
    let exner = layer_exner(gas, sigma, ps);
    let t = [theta[0] * exner[0], theta[1] * exner[1]];
    // Thickness between the levels whose Exner functions are the layer means.
    let p = exner.map(|e| REFERENCE_PRESSURE * fms_math::pow(e, 1.0 / gas.kappa()));
    let dz = gas.gas_constant * 0.5 * (t[0] + t[1]) / gravity * fms_math::ln(p[1] / p[0]);
    let drop = params.lapse_rate * dz;
    if t[1] - t[0] <= drop {
        return None;
    }
    let ds = [sigma, 1.0 - sigma];
    let enthalpy = ds[0] * t[0] + ds[1] * t[1];
    // T₀ = T₁ − drop with the same enthalpy.
    let t1 = enthalpy + ds[0] * drop;
    Some([(t1 - drop) / exner[0], t1 / exner[1]])
}

/// Heating rate of `Θₖ = p_s θₖ` (Pa K s⁻¹) from `q` W m⁻² absorbed by a
/// layer.
fn theta_tendency(gas: Gas, gravity: f64, dsigma: f64, exner: f64, q: f64) -> f64 {
    gravity * q / (gas.heat_capacity * dsigma * exner)
}

/// The radiative–convective forcing of the two-layer model: radiation, a
/// slab surface, surface heat fluxes, convection and boundary-layer friction.
pub struct ColumnPhysics<'a> {
    params: PhysicsParams,
    gas: Gas,
    sigma: f64,
    gravity: f64,
    ops: Operators<'a>,
    planet: &'a Planet,
    /// Model time (s since t = 0).
    time: f64,
    next_radiation: f64,
    surface_temperature: Vec<f64>,
    /// Last radiation call.
    radiation: Vec<Fluxes>,
    /// Sensible and latent heat flux per cell, from the last step.
    surface: Vec<(f64, f64)>,
    absorbed: f64,
    outgoing: f64,
}

impl<'a> ColumnPhysics<'a> {
    /// Physics for `model` on `planet`, starting at time `time` (s) from
    /// `state` with surface temperatures `surface_temperature` (K).
    pub fn new(
        model: &TwoLayer<'a>,
        planet: &'a Planet,
        params: PhysicsParams,
        time: f64,
        state: &TwoLayerState,
        surface_temperature: Vec<f64>,
    ) -> Self {
        let n = state.ps.len();
        assert_eq!(surface_temperature.len(), n);
        let mut physics = Self {
            params,
            gas: model.gas(),
            sigma: model.sigma(),
            gravity: model.gravity(),
            ops: *model.operators(),
            planet,
            time,
            next_radiation: time,
            surface_temperature,
            radiation: vec![Fluxes::default(); n],
            surface: vec![(0.0, 0.0); n],
            absorbed: 0.0,
            outgoing: 0.0,
        };
        physics.update_surface_fluxes(state);
        physics.radiate(state);
        physics
    }

    /// The parameters.
    pub fn params(&self) -> &PhysicsParams {
        &self.params
    }

    /// Surface temperature per cell (K).
    pub fn surface_temperature(&self) -> &[f64] {
        &self.surface_temperature
    }

    /// The last radiation call's fluxes per cell.
    pub fn fluxes(&self) -> &[Fluxes] {
        &self.radiation
    }

    /// Model time (s).
    pub fn time(&self) -> f64 {
        self.time
    }

    fn update_surface_fluxes(&mut self, state: &TwoLayerState) {
        let (mesh, geometry) = (self.ops.mesh(), self.ops.geometry());
        let winds = geometry.reconstruct_vectors(mesh, &state.u[1]);
        let centers = mesh.cell_centers();
        let (params, gas, sigma, gravity) = (&self.params, self.gas, self.sigma, self.gravity);
        let ts = &self.surface_temperature;
        self.surface
            .par_iter_mut()
            .with_min_len(MIN_CHUNK)
            .enumerate()
            .for_each(|(i, f)| {
                let (e, n) = east_north(centers[i], winds[i]);
                let theta = state.theta_ps[1][i] / state.ps[i];
                let air = surface_air_temperature(params, gas, sigma, gravity, state.ps[i], theta);
                *f = surface_fluxes(params, gas, state.ps[i], ts[i], air, (e * e + n * n).sqrt());
            });
    }

    /// Radiation for the interval starting now, with sunlight averaged over
    /// it.
    fn radiate(&mut self, state: &TwoLayerState) {
        const SAMPLES: usize = 6;
        let centers = self.ops.mesh().cell_centers();
        let n = centers.len();
        let mut insolation = vec![0.0; n];
        let mut sample = vec![0.0; n];
        let h = self.params.radiation_interval / SAMPLES as f64;
        for j in 0..SAMPLES {
            self.planet
                .insolation_into(centers, self.time + (j as f64 + 0.5) * h, &mut sample);
            for (a, s) in insolation.iter_mut().zip(&sample) {
                *a += s / SAMPLES as f64;
            }
        }
        let (params, gas, sigma, gravity) = (&self.params, self.gas, self.sigma, self.gravity);
        let ts = &self.surface_temperature;
        self.radiation
            .par_iter_mut()
            .with_min_len(MIN_CHUNK)
            .enumerate()
            .for_each(|(i, f)| {
                *f = column_fluxes(
                    &params.radiation,
                    gas,
                    &Column {
                        ps: state.ps[i],
                        surface_temperature: ts[i],
                        theta: [
                            state.theta_ps[0][i] / state.ps[i],
                            state.theta_ps[1][i] / state.ps[i],
                        ],
                        sigma,
                        gravity,
                        insolation: insolation[i],
                        co2: params.co2,
                    },
                );
            });
        self.next_radiation = self.time + self.params.radiation_interval;
    }
}

impl Forcing for ColumnPhysics<'_> {
    fn add_tendency(&self, state: &TwoLayerState, out: &mut TwoLayerState) {
        let (params, gas, sigma, gravity) = (&self.params, self.gas, self.sigma, self.gravity);
        let ds = [sigma, 1.0 - sigma];
        let [t0, t1] = &mut out.theta_ps;
        t0.par_iter_mut()
            .with_min_len(MIN_CHUNK)
            .zip(t1.par_iter_mut())
            .enumerate()
            .for_each(|(i, (d0, d1))| {
                let ps = state.ps[i];
                let exner = layer_exner(gas, sigma, ps);
                let f = &self.radiation[i];
                let (sensible, latent) = self.surface[i];
                let heating = [
                    f.layers[0] + ds[0] * latent,
                    f.layers[1] + sensible + ds[1] * latent,
                ];
                *d0 += theta_tendency(gas, gravity, ds[0], exner[0], heating[0]);
                *d1 += theta_tendency(gas, gravity, ds[1], exner[1], heating[1]);
                let theta = [state.theta_ps[0][i] / ps, state.theta_ps[1][i] / ps];
                if let Some(target) = convective_target(params, gas, sigma, gravity, ps, theta) {
                    *d0 += ps * (target[0] - theta[0]) / params.convection_time;
                    *d1 += ps * (target[1] - theta[1]) / params.convection_time;
                }
            });
        let rate = self.params.friction;
        if rate > 0.0 {
            out.u[1]
                .par_iter_mut()
                .with_min_len(MIN_CHUNK)
                .zip(state.u[1].par_iter())
                .for_each(|(d, &u)| *d -= rate * u);
        }
    }

    fn advance(&mut self, state: &TwoLayerState, dt: f64) {
        // The slab and the budget use the same stored fluxes the atmosphere
        // was heated with.
        let areas = &self.ops.geometry().cell_areas;
        let r2 = self.ops.radius() * self.ops.radius();
        let capacity = self.params.heat_capacity();
        let (radiation, surface) = (&self.radiation, &self.surface);
        self.surface_temperature
            .par_iter_mut()
            .with_min_len(MIN_CHUNK)
            .enumerate()
            .for_each(|(i, t)| {
                let (sensible, latent) = surface[i];
                *t += dt * (radiation[i].surface() - sensible - latent) / capacity;
            });
        let (mut absorbed, mut outgoing) = (0.0, 0.0);
        for (f, a) in radiation.iter().zip(areas) {
            absorbed += a * f.absorbed;
            outgoing += a * f.outgoing;
        }
        self.absorbed += dt * r2 * absorbed;
        self.outgoing += dt * r2 * outgoing;
        self.time += dt;
        self.update_surface_fluxes(state);
        if self.time >= self.next_radiation - 1e-6 {
            self.radiate(state);
        }
    }

    fn budget(&self) -> Option<ForcingBudget> {
        let areas = &self.ops.geometry().cell_areas;
        let r2 = self.ops.radius() * self.ops.radius();
        let stored: f64 = self
            .surface_temperature
            .iter()
            .zip(areas)
            .map(|(t, a)| a * r2 * self.params.heat_capacity() * t)
            .sum();
        Some(ForcingBudget {
            absorbed: self.absorbed,
            outgoing: self.outgoing,
            stored,
        })
    }

    fn cell_fields(&self) -> Vec<CellField> {
        let pick = |f: fn(&Fluxes) -> f64| self.radiation.iter().map(f).collect();
        vec![
            CellField {
                name: "rad.surface_temperature",
                description: "Surface temperature",
                units: "K",
                values: self.surface_temperature.clone(),
            },
            CellField {
                name: "rad.outgoing_longwave",
                description: "Outgoing longwave radiation",
                units: "W m-2",
                values: pick(|f| f.outgoing),
            },
            CellField {
                name: "rad.absorbed_shortwave",
                description: "Absorbed sunlight",
                units: "W m-2",
                values: pick(|f| f.absorbed),
            },
            CellField {
                name: "rad.water_vapour",
                description: "Water-vapour column (fixed relative humidity)",
                units: "kg m-2",
                values: pick(|f| f.water[0] + f.water[1]),
            },
            CellField {
                name: "rad.latent_heat_flux",
                description: "Surface latent heat flux",
                units: "W m-2",
                values: self.surface.iter().map(|s| s.1).collect(),
            },
        ]
    }
}

/// A single-column radiative–convective model: the same physics as
/// [`ColumnPhysics`] without dynamics, under a fixed insolation. For
/// calibration and quick experiments (CO₂, pressure, gravity, sunlight).
#[derive(Clone, Copy, Debug)]
pub struct SingleColumn {
    /// Parameters (the slab depth only sets how fast it converges).
    pub params: PhysicsParams,
    /// Gas constants.
    pub gas: Gas,
    /// Interface σ.
    pub sigma: f64,
    /// Gravity (m s⁻²).
    pub gravity: f64,
    /// Surface pressure (Pa).
    pub ps: f64,
    /// Mean insolation (W m⁻²).
    pub insolation: f64,
    /// Fixed wind for the surface fluxes (m s⁻¹).
    pub wind: f64,
    /// Layer potential temperatures (K).
    pub theta: [f64; 2],
    /// Surface temperature (K).
    pub surface_temperature: f64,
}

/// A single column's equilibrium.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Equilibrium {
    /// Radiative fluxes.
    pub fluxes: Fluxes,
    /// Surface temperature (K).
    pub surface_temperature: f64,
    /// Layer temperatures, upper first (K).
    pub temperature: [f64; 2],
    /// Sensible and latent heat flux (W m⁻²).
    pub surface_fluxes: (f64, f64),
    /// Steps taken.
    pub steps: u64,
}

impl SingleColumn {
    /// Cloud albedo for a single global-mean column: lower than the 3-D
    /// model's (0.33) because one column under the planet's mean sunlight
    /// isn't the planet's mean. The physics is nonlinear in temperature
    /// (water vapour, `σT⁴`), and the 3-D model's warm tropics hold more
    /// water. With it, the column reaches Earth's global mean as the 3-D
    /// model does.
    pub const SINGLE_COLUMN_CLOUD_ALBEDO: f64 = 0.30;

    /// An Earth-like column under the global-mean sunlight of a planet with
    /// mean stellar flux `s0`, with [`Self::SINGLE_COLUMN_CLOUD_ALBEDO`].
    pub fn earth(s0: f64) -> Self {
        let mut params = PhysicsParams {
            slab_depth: 1.0,
            ..PhysicsParams::default()
        };
        params.radiation.cloud_albedo = Self::SINGLE_COLUMN_CLOUD_ALBEDO;
        Self {
            params,
            gas: Gas::EARTH_AIR,
            sigma: 0.5,
            gravity: 9.81,
            ps: 101_325.0,
            insolation: s0 / 4.0,
            wind: 5.0,
            theta: [330.0, 290.0],
            surface_temperature: 288.0,
        }
    }

    /// Fluxes of the current state.
    pub fn fluxes(&self) -> Fluxes {
        column_fluxes(
            &self.params.radiation,
            self.gas,
            &Column {
                ps: self.ps,
                surface_temperature: self.surface_temperature,
                theta: self.theta,
                sigma: self.sigma,
                gravity: self.gravity,
                insolation: self.insolation,
                co2: self.params.co2,
            },
        )
    }

    /// Sensible and latent heat flux (W m⁻²).
    pub fn surface_fluxes(&self) -> (f64, f64) {
        let air = surface_air_temperature(
            &self.params,
            self.gas,
            self.sigma,
            self.gravity,
            self.ps,
            self.theta[1],
        );
        surface_fluxes(
            &self.params,
            self.gas,
            self.ps,
            self.surface_temperature,
            air,
            self.wind,
        )
    }

    /// Set the layer temperatures (K, upper first).
    pub fn set_temperatures(&mut self, t: [f64; 2]) {
        let exner = layer_exner(self.gas, self.sigma, self.ps);
        self.theta = [t[0] / exner[0], t[1] / exner[1]];
    }

    /// One step of `dt` seconds (forward Euler, radiation every step).
    pub fn step(&mut self, dt: f64) {
        let f = self.fluxes();
        let (sensible, latent) = self.surface_fluxes();
        let exner = layer_exner(self.gas, self.sigma, self.ps);
        let ds = [self.sigma, 1.0 - self.sigma];
        let heating = [
            f.layers[0] + ds[0] * latent,
            f.layers[1] + sensible + ds[1] * latent,
        ];
        for k in 0..2 {
            self.theta[k] +=
                dt * theta_tendency(self.gas, self.gravity, ds[k], exner[k], heating[k]) / self.ps;
        }
        if let Some(target) = convective_target(
            &self.params,
            self.gas,
            self.sigma,
            self.gravity,
            self.ps,
            self.theta,
        ) {
            let a = (dt / self.params.convection_time).min(1.0);
            for (theta, target) in self.theta.iter_mut().zip(target) {
                *theta += a * (target - *theta);
            }
        }
        self.surface_temperature +=
            dt * (f.surface() - sensible - latent) / self.params.heat_capacity();
    }

    /// Step until the top-of-atmosphere imbalance and the surface
    /// temperature settle (or 50 model years pass).
    pub fn equilibrium(mut self) -> Equilibrium {
        let dt = 3600.0;
        let mut steps = 0;
        let mut last = f64::INFINITY;
        while steps < (50.0 * 365.0 * 24.0) as u64 {
            for _ in 0..240 {
                self.step(dt);
            }
            steps += 240;
            let f = self.fluxes();
            let change = (self.surface_temperature - last).abs();
            last = self.surface_temperature;
            if !self.surface_temperature.is_finite() {
                break;
            }
            if change < 1e-6 && (f.absorbed - f.outgoing).abs() < 1e-3 {
                break;
            }
        }
        let exner = layer_exner(self.gas, self.sigma, self.ps);
        Equilibrium {
            fluxes: self.fluxes(),
            surface_temperature: self.surface_temperature,
            temperature: [self.theta[0] * exner[0], self.theta[1] * exner[1]],
            surface_fluxes: self.surface_fluxes(),
            steps,
        }
    }
}
