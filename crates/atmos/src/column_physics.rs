//! Column physics for a radiatively driven planet: [`crate::radiation`],
//! [`crate::moisture`], and a temporary slab surface with bulk sensible and
//! latent heat fluxes. [`ColumnPhysics`] is the 3-D forcing; [`SingleColumn`]
//! is a 1-D radiative–convective model with fixed relative humidity, kept
//! to check the radiation against Earth's global-mean column.
//!
//! # Pieces (3-D)
//!
//! - **Surface**: one slab of water of fixed depth everywhere (an
//!   aquaplanet; land and ocean come in M2 steps 10–11). It absorbs sunlight
//!   and longwave, emits as a black body, and exchanges heat and water with
//!   the lower layer by bulk formulas: sensible `ρ c_p C V (T_s − T_a)` and
//!   evaporation `ρ C V (q*(T_s) − q_a)`, where `T_a` is the lower layer's
//!   temperature extrapolated to the surface along
//!   [`PhysicsParams::lapse_rate`], `q_a` the lower layer's humidity profile
//!   at the surface ([`crate::moisture`]), and `V` includes a gust speed.
//!   The slab pays the latent heat of what evaporates, and of fusion for the
//!   snow that falls on it and melts.
//! - **Moisture**: convection, large-scale condensation, re-evaporation,
//!   rain or snow, and diagnosed clouds ([`crate::moisture`]).
//! - **Friction**: the Held–Suarez boundary-layer drag on the lower layer,
//!   until M2 step 10 adds surface drag.
//!
//! # Time stepping
//!
//! Radiation runs every [`PhysicsParams::radiation_interval`] with sunlight
//! averaged over the coming interval and the clouds of the latest moist
//! physics; surface fluxes and the slab update every step
//! ([`Forcing::advance`]). The atmosphere and the slab use the same stored
//! fluxes and the same precipitation, so the planet's energy and water
//! budgets close exactly: what the atmosphere and slab gain is what was
//! absorbed from the Sun minus what was emitted to space (and the kinetic
//! energy friction and diffusion remove), and the atmosphere's water
//! changes by evaporation minus precipitation.

use fms_geo::Operators;
use fms_geo::operators::MIN_CHUNK;
use fms_geo::zonal::east_north;
use fms_planet::Planet;
use rayon::prelude::*;

use crate::forcing::{CellField, Forcing, ForcingBudget};
use crate::moisture::{
    FUSION_HEAT, MoistColumn, MoistTendency, MoistureParams, Setting, moist_physics,
};
use crate::radiation::{
    Cloud, Column, Fluxes, Humidity, RadiationParams, STEFAN_BOLTZMANN, column_fluxes,
    saturation_humidity, saturation_pressure, specific_humidity,
};
use crate::two_layer::{Gas, REFERENCE_PRESSURE, TwoLayer, TwoLayerState};

pub use crate::moisture::LATENT_HEAT;

const DAY: f64 = 86_400.0;
/// Volumetric heat capacity of water (J m⁻³ K⁻¹).
const WATER_HEAT_CAPACITY: f64 = 1000.0 * 4186.0;

/// Column-physics parameters. [`Default`] is Earth-like.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PhysicsParams {
    /// Radiation.
    pub radiation: RadiationParams,
    /// Moisture (3-D only).
    pub moisture: MoistureParams,
    /// CO₂ (ppm).
    pub co2: f64,
    /// Slab depth (m of water).
    pub slab_depth: f64,
    /// Bulk transfer coefficient for heat and moisture.
    pub transfer: f64,
    /// Gust speed added in quadrature to the wind (m s⁻¹).
    pub gust: f64,
    /// Lapse rate (K m⁻¹) from the lower layer's mean down to the surface
    /// air, and the single column's critical lapse rate for convection
    /// between the layer means. Earth's observed layer means (272 and 226 K,
    /// 10.5 km apart) sit at 4.4 K/km: less than the tropospheric 6.5 K/km
    /// because the upper layer includes the stratosphere.
    pub lapse_rate: f64,
    /// The single column's convective relaxation time (s).
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
            moisture: MoistureParams::default(),
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

/// Bulk exchange `ρ C V` (kg m⁻² s⁻¹) over the surface for surface air at
/// temperature `air` and wind speed `wind`.
fn exchange(params: &PhysicsParams, gas: Gas, ps: f64, air: f64, wind: f64) -> f64 {
    let density = ps / (gas.gas_constant * air);
    let speed = (wind * wind + params.gust * params.gust).sqrt();
    density * params.transfer * speed
}

/// Fluxes between the surface and the lower layer in one cell.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
struct SurfaceFlux {
    /// Sensible heat (W m⁻², upward).
    sensible: f64,
    /// Evaporation (kg m⁻² s⁻¹).
    evaporation: f64,
}

/// The radiative–convective forcing of the two-layer model: radiation, moist
/// physics, a slab surface, surface fluxes and boundary-layer friction. Needs
/// a moist state.
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
    /// Last radiation call, all-sky and clear-sky.
    radiation: Vec<Fluxes>,
    clear: Vec<Fluxes>,
    /// Surface fluxes, from the last step.
    surface: Vec<SurfaceFlux>,
    /// The moist physics of the last tendency evaluation.
    moist: Vec<MoistTendency>,
    absorbed: f64,
    outgoing: f64,
    evaporated: f64,
    precipitated: f64,
}

impl<'a> ColumnPhysics<'a> {
    /// Physics for `model` on `planet`, starting at time `time` (s) from
    /// `state` with surface temperatures `surface_temperature` (K).
    ///
    /// # Panics
    /// If `state` is dry.
    pub fn new(
        model: &TwoLayer<'a>,
        planet: &'a Planet,
        params: PhysicsParams,
        time: f64,
        state: &TwoLayerState,
        surface_temperature: Vec<f64>,
    ) -> Self {
        assert!(state.is_moist(), "column physics needs a moist state");
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
            clear: vec![Fluxes::default(); n],
            surface: vec![SurfaceFlux::default(); n],
            moist: vec![MoistTendency::default(); n],
            absorbed: 0.0,
            outgoing: 0.0,
            evaporated: 0.0,
            precipitated: 0.0,
        };
        // Clouds for the first radiation call; radiation (clear sky) for the
        // first convection's stratosphere temperature.
        physics.radiate(state);
        physics.update_moisture(state);
        physics.radiate(state);
        physics.update_surface_fluxes(state);
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

    /// The last moist physics per cell.
    pub fn moisture(&self) -> &[MoistTendency] {
        &self.moist
    }

    /// Model time (s).
    pub fn time(&self) -> f64 {
        self.time
    }

    fn setting(&self) -> Setting<'_> {
        Setting {
            params: &self.params.moisture,
            radiation: &self.params.radiation,
            gas: self.gas,
            sigma: self.sigma,
            gravity: self.gravity,
        }
    }

    /// One column's view for the moist physics.
    fn moist_column(&self, state: &TwoLayerState, i: usize) -> MoistColumn {
        let ps = state.ps[i];
        let exner = layer_exner(self.gas, self.sigma, ps);
        let theta = [state.theta_ps[0][i] / ps, state.theta_ps[1][i] / ps];
        // The stratosphere sits near the skin temperature of the column's
        // outgoing longwave.
        let olr = self.radiation[i].outgoing.max(1.0);
        MoistColumn {
            ps,
            temperature: [theta[0] * exner[0], theta[1] * exner[1]],
            humidity: [state.q_ps[0][i] / ps, state.q_ps[1][i] / ps],
            surface_air: surface_air_temperature(
                &self.params,
                self.gas,
                self.sigma,
                self.gravity,
                ps,
                theta[1],
            ),
            stratosphere: fms_math::pow(olr / (2.0 * STEFAN_BOLTZMANN), 0.25),
        }
    }

    fn update_moisture(&mut self, state: &TwoLayerState) {
        let moist: Vec<MoistTendency> = (0..state.ps.len())
            .into_par_iter()
            .with_min_len(MIN_CHUNK / 16)
            .map(|i| moist_physics(&self.setting(), &self.moist_column(state, i)))
            .collect();
        self.moist = moist;
    }

    fn update_surface_fluxes(&mut self, state: &TwoLayerState) {
        let (mesh, geometry) = (self.ops.mesh(), self.ops.geometry());
        let winds = geometry.reconstruct_vectors(mesh, &state.u[1]);
        let centers = mesh.cell_centers();
        let surface: Vec<SurfaceFlux> = (0..state.ps.len())
            .into_par_iter()
            .with_min_len(MIN_CHUNK / 16)
            .map(|i| {
                let (e, n) = east_north(centers[i], winds[i]);
                let column = self.moist_column(state, i);
                let ps = column.ps;
                let k = exchange(
                    &self.params,
                    self.gas,
                    ps,
                    column.surface_air,
                    (e * e + n * n).sqrt(),
                );
                let q_air = crate::moisture::surface_humidity(
                    &self.setting(),
                    ps,
                    column.temperature,
                    column.humidity[1],
                    column.surface_air,
                );
                let ts = self.surface_temperature[i];
                SurfaceFlux {
                    sensible: k * self.gas.heat_capacity * (ts - column.surface_air),
                    evaporation: k * (saturation_humidity(ts, ps) - q_air),
                }
            })
            .collect();
        self.surface = surface;
    }

    /// Radiation for the interval starting now, with sunlight averaged over
    /// it, all-sky and clear-sky.
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
        let clouds = &params.moisture.clouds;
        let moist = &self.moist;
        let both: Vec<(Fluxes, Fluxes)> = (0..n)
            .into_par_iter()
            .with_min_len(MIN_CHUNK / 16)
            .map(|i| {
                let ps = state.ps[i];
                let cover = moist[i].clouds;
                let mut column = Column {
                    ps,
                    surface_temperature: self.surface_temperature[i],
                    theta: [state.theta_ps[0][i] / ps, state.theta_ps[1][i] / ps],
                    sigma,
                    gravity,
                    insolation: insolation[i],
                    co2: params.co2,
                    humidity: Humidity::Layers([state.q_ps[0][i] / ps, state.q_ps[1][i] / ps]),
                    clouds: std::array::from_fn(|c| Cloud {
                        sigma: clouds.sigma[c],
                        reflectance: cover[c] * clouds.albedo[c],
                        emissivity: cover[c] * clouds.emissivity[c],
                    }),
                };
                let all = column_fluxes(&params.radiation, gas, &column);
                column.clouds = [Cloud::NONE; 2];
                (all, column_fluxes(&params.radiation, gas, &column))
            })
            .collect();
        (self.radiation, self.clear) = both.into_iter().unzip();
        self.next_radiation = self.time + self.params.radiation_interval;
    }
}

impl Forcing for ColumnPhysics<'_> {
    fn add_tendency(&mut self, state: &TwoLayerState, out: &mut TwoLayerState) {
        assert!(
            state.is_moist() && out.is_moist(),
            "column physics needs a moist state"
        );
        self.update_moisture(state);
        let (gas, sigma, gravity) = (self.gas, self.sigma, self.gravity);
        let ds = [sigma, 1.0 - sigma];
        let (radiation, surface, moist) = (&self.radiation, &self.surface, &self.moist);
        let [t0, t1] = &mut out.theta_ps;
        let [q0, q1] = &mut out.q_ps;
        t0.par_iter_mut()
            .with_min_len(MIN_CHUNK)
            .zip(t1.par_iter_mut())
            .zip(q0.par_iter_mut())
            .zip(q1.par_iter_mut())
            .enumerate()
            .for_each(|(i, (((d0, d1), dq0), dq1))| {
                let ps = state.ps[i];
                let exner = layer_exner(gas, sigma, ps);
                let f = &radiation[i];
                let m = &moist[i];
                let heating = [f.layers[0], f.layers[1] + surface[i].sensible];
                *d0 += theta_tendency(gas, gravity, ds[0], exner[0], heating[0])
                    + ps * m.temperature[0] / exner[0];
                *d1 += theta_tendency(gas, gravity, ds[1], exner[1], heating[1])
                    + ps * m.temperature[1] / exner[1];
                *dq0 += ps * m.humidity[0];
                *dq1 += ps * m.humidity[1] + gravity * surface[i].evaporation / ds[1];
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
        // The slab and the budgets use the same stored fluxes and
        // precipitation the atmosphere was given.
        let areas = &self.ops.geometry().cell_areas;
        let r2 = self.ops.radius() * self.ops.radius();
        let capacity = self.params.heat_capacity();
        let (radiation, surface, moist) = (&self.radiation, &self.surface, &self.moist);
        self.surface_temperature
            .par_iter_mut()
            .with_min_len(MIN_CHUNK)
            .enumerate()
            .for_each(|(i, t)| {
                let s = surface[i];
                let loss = s.sensible + LATENT_HEAT * s.evaporation + FUSION_HEAT * moist[i].snow;
                *t += dt * (radiation[i].surface() - loss) / capacity;
            });
        let (mut absorbed, mut outgoing, mut evaporated, mut precipitated) = (0.0, 0.0, 0.0, 0.0);
        for i in 0..areas.len() {
            let a = areas[i];
            absorbed += a * radiation[i].absorbed;
            outgoing += a * radiation[i].outgoing;
            evaporated += a * surface[i].evaporation;
            precipitated += a * (moist[i].rain + moist[i].snow);
        }
        self.absorbed += dt * r2 * absorbed;
        self.outgoing += dt * r2 * outgoing;
        self.evaporated += dt * r2 * evaporated;
        self.precipitated += dt * r2 * precipitated;
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
            evaporated: self.evaporated,
            precipitated: self.precipitated,
        })
    }

    fn cell_fields(&self) -> Vec<CellField> {
        let rad = |f: fn(&Fluxes) -> f64| self.radiation.iter().map(f).collect();
        let moist = |f: fn(&MoistTendency) -> f64| self.moist.iter().map(f).collect();
        let field = |name, description, units, values| CellField {
            name,
            description,
            units,
            values,
        };
        vec![
            field(
                "rad.surface_temperature",
                "Surface temperature",
                "K",
                self.surface_temperature.clone(),
            ),
            field(
                "rad.outgoing_longwave",
                "Outgoing longwave radiation",
                "W m-2",
                rad(|f| f.outgoing),
            ),
            field(
                "rad.absorbed_shortwave",
                "Absorbed sunlight",
                "W m-2",
                rad(|f| f.absorbed),
            ),
            field(
                "rad.surface_shortwave",
                "Sunlight absorbed by the surface",
                "W m-2",
                rad(|f| f.surface_shortwave),
            ),
            field(
                "rad.surface_longwave",
                "Net longwave loss of the surface",
                "W m-2",
                rad(|f| f.surface_up - f.surface_down),
            ),
            field(
                "rad.clear_surface_longwave",
                "Net longwave loss of the surface under clear skies",
                "W m-2",
                self.clear
                    .iter()
                    .map(|f| f.surface_up - f.surface_down)
                    .collect(),
            ),
            field(
                "rad.shortwave_cloud_effect",
                "Shortwave cloud radiative effect",
                "W m-2",
                self.radiation
                    .iter()
                    .zip(&self.clear)
                    .map(|(a, c)| a.absorbed - c.absorbed)
                    .collect(),
            ),
            field(
                "rad.longwave_cloud_effect",
                "Longwave cloud radiative effect",
                "W m-2",
                self.radiation
                    .iter()
                    .zip(&self.clear)
                    .map(|(a, c)| c.outgoing - a.outgoing)
                    .collect(),
            ),
            field(
                "atm.water_vapour",
                "Water-vapour column",
                "kg m-2",
                rad(|f| f.water[0] + f.water[1]),
            ),
            field(
                "atm.precipitation",
                "Precipitation",
                "mm day-1",
                moist(|m| (m.rain + m.snow) * DAY),
            ),
            field(
                "atm.snowfall",
                "Snowfall (water equivalent)",
                "mm day-1",
                moist(|m| m.snow * DAY),
            ),
            field(
                "atm.convective_precipitation",
                "Convective precipitation",
                "mm day-1",
                moist(|m| m.convective * DAY),
            ),
            field(
                "atm.evaporation",
                "Evaporation",
                "mm day-1",
                self.surface.iter().map(|s| s.evaporation * DAY).collect(),
            ),
            field(
                "atm.latent_heat_flux",
                "Surface latent heat flux",
                "W m-2",
                self.surface
                    .iter()
                    .map(|s| s.evaporation * LATENT_HEAT)
                    .collect(),
            ),
            field(
                "atm.sensible_heat_flux",
                "Surface sensible heat flux",
                "W m-2",
                self.surface.iter().map(|s| s.sensible).collect(),
            ),
            field(
                "atm.high_cloud",
                "High cloud cover",
                "1",
                moist(|m| m.clouds[0]),
            ),
            field(
                "atm.low_cloud",
                "Low cloud cover",
                "1",
                moist(|m| m.clouds[1]),
            ),
            field(
                "atm.cloud_cover",
                "Total cloud cover (random overlap)",
                "1",
                moist(|m| 1.0 - (1.0 - m.clouds[0]) * (1.0 - m.clouds[1])),
            ),
            field(
                "atm.relative_humidity_upper",
                "Upper-layer relative humidity",
                "1",
                moist(|m| m.relative_humidity[0]),
            ),
            field(
                "atm.relative_humidity_lower",
                "Lower-layer relative humidity",
                "1",
                moist(|m| m.relative_humidity[1]),
            ),
        ]
    }
}

/// A moist start for the 3-D model on `planet`: everywhere the single
/// column's equilibrium (with `params`' CO₂) at rest, with its fixed-profile
/// humidity, plus 0.1 K of seeded noise in θ to break the symmetry. Returns
/// the state and the surface temperature (K).
pub fn initial_state(
    model: &TwoLayer,
    planet: &Planet,
    params: &PhysicsParams,
) -> (TwoLayerState, f64) {
    let mut column = SingleColumn::earth(planet.mean_flux());
    column.params.co2 = params.co2;
    column.ps = planet.params().surface_pressure;
    column.gravity = model.gravity();
    let eq = column.equilibrium();
    let (gas, sigma, ps) = (model.gas(), model.sigma(), column.ps);
    let exner = layer_exner(gas, sigma, ps);
    let theta0 = [eq.temperature[0] / exner[0], eq.temperature[1] / exner[1]];
    let q = fixed_humidity(params, gas, sigma, model.gravity(), ps, theta0);
    let mesh = model.operators().mesh();
    let (nc, ne) = (mesh.n_cells(), mesh.n_edges());
    let mut x = 0x2545_F491_4F6C_DD1Du64;
    let theta: [Vec<f64>; 2] = std::array::from_fn(|k| {
        (0..nc)
            .map(|_| {
                x ^= x << 13;
                x ^= x >> 7;
                x ^= x << 17;
                theta0[k] + 0.1 * ((x >> 11) as f64 / (1u64 << 52) as f64 - 1.0)
            })
            .collect()
    });
    let state = TwoLayerState::new(vec![ps; nc], [vec![0.0; ne], vec![0.0; ne]], theta)
        .with_humidity([vec![q[0]; nc], vec![q[1]; nc]]);
    (state, eq.surface_temperature)
}

/// Initial humidity (kg kg⁻¹, per layer) for a column at fixed relative
/// humidity: the Manabe–Wetherald profile of the single column.
pub fn fixed_humidity(
    params: &PhysicsParams,
    gas: Gas,
    sigma: f64,
    gravity: f64,
    ps: f64,
    theta: [f64; 2],
) -> [f64; 2] {
    let f = column_fluxes(
        &params.radiation,
        gas,
        &Column {
            ps,
            surface_temperature: 288.0,
            theta,
            sigma,
            gravity,
            insolation: 0.0,
            co2: params.co2,
            humidity: Humidity::Fixed,
            clouds: [Cloud::NONE; 2],
        },
    );
    let ds = [sigma, 1.0 - sigma];
    [0, 1].map(|k| f.water[k] * gravity / (ds[k] * ps))
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
    /// Albedo of the single column's fixed cloud deck, calibrated so the
    /// column reaches Earth's global-mean climate (with the clear air's
    /// Rayleigh scattering, [`RadiationParams::rayleigh`], on top). One
    /// column under the planet's mean sunlight isn't the planet's mean (the
    /// physics is nonlinear in temperature), so this is not the 3-D model's
    /// cloud albedo; the 3-D model diagnoses its clouds
    /// ([`crate::moisture`]).
    pub const SINGLE_COLUMN_CLOUD_ALBEDO: f64 = 0.27;

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
                humidity: Humidity::Fixed,
                clouds: [
                    Cloud::deck(self.sigma, self.params.radiation.cloud_albedo),
                    Cloud::NONE,
                ],
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
