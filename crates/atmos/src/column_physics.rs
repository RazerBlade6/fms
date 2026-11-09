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
//!   and longwave, emits as a black body, and pays the latent heat of what
//!   evaporates, and of fusion for the snow that falls on it and melts.
//! - **Boundary layer**: the bottom layer is a well-mixed boundary layer
//!   ([`crate::surface_layer`]): surface air has its potential temperature
//!   and specific humidity. Bulk formulas with stability-dependent
//!   coefficients and an ocean roughness give the surface's sensible heat,
//!   evaporation and drag on the bottom layer's wind; the layer above is
//!   entrained at mixed-layer theory's rate, exchanging dry static energy,
//!   water and momentum, and shallow cumulus carry its water and momentum
//!   up where its air saturates below its top.
//! - **Moisture**: convection, large-scale condensation, re-evaporation,
//!   rain or snow, and diagnosed clouds ([`crate::moisture`]).
//! - **Dissipation**: the kinetic energy drag and entrainment remove heats
//!   the boundary layer ([`crate::dissipation`]).
//!
//! # Time stepping
//!
//! Radiation runs every [`PhysicsParams::radiation_interval`] with sunlight
//! averaged over the coming interval and the clouds of the latest moist
//! physics; surface fluxes and the slab update every step
//! ([`Forcing::advance`]). The atmosphere and the slab use the same stored
//! fluxes and the same precipitation, so the planet's energy and water
//! budgets close exactly: what the atmosphere and slab gain is what was
//! absorbed from the Sun minus what was emitted to space (friction's lost
//! kinetic energy comes back as heat), and the atmosphere's water changes
//! by evaporation minus precipitation.

use fms_geo::Operators;
use fms_geo::operators::MIN_CHUNK;
use fms_geo::zonal::east_north;
use fms_planet::Planet;
use rayon::prelude::*;

use crate::dissipation::KineticHeating;
use crate::forcing::{CellField, Forcing, ForcingBudget};
use crate::layered::{Gas, LayeredModel, LayeredState, REFERENCE_PRESSURE};
use crate::levels::{Levels, MAX_LAYERS};
use crate::moisture::{
    FUSION_HEAT, MoistColumn, MoistTendency, MoistureParams, Setting, lift, moist_physics,
};
use crate::radiation::{
    Cloud, Column, Fluxes, Humidity, Profile, RadiationParams, STEFAN_BOLTZMANN, column_fluxes,
    saturation_humidity, saturation_pressure, specific_humidity,
};
use crate::surface_layer::{
    Surface, SurfaceAir, SurfaceLayerParams, convective_velocity, entrainment_velocity, exchange,
};

pub use crate::moisture::LATENT_HEAT;

const DAY: f64 = 86_400.0;

/// The climate's interior layer interfaces (σ): four layers. The bottom one
/// (~1.4 km) is the boundary layer. The top is split at 250 hPa so that the
/// Hadley outflow (σ 0.25–0.5) is not averaged with the stratosphere, which
/// inflates its moist static energy and suppresses the subtropical eddy
/// momentum flux.
pub const CLIMATE_LEVELS: [f64; 3] = [0.25, 0.5, 0.85];
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
    /// The surface layer and boundary-layer top (3-D).
    pub surface_layer: SurfaceLayerParams,
    /// The single column's bulk transfer coefficient for heat and moisture.
    pub transfer: f64,
    /// The single column's gust speed, added in quadrature to the wind
    /// (m s⁻¹).
    pub gust: f64,
    /// The single column's lapse rate (K m⁻¹) from its lower layer's mean
    /// down to the surface air, and its critical lapse rate for convection
    /// between its two layer means. Earth's observed layer means (272 and
    /// 226 K, 10.5 km apart) sit at 4.4 K/km: less than the tropospheric
    /// 6.5 K/km because the upper layer includes the stratosphere.
    pub lapse_rate: f64,
    /// The single column's convective relaxation time (s).
    pub convection_time: f64,
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
            surface_layer: SurfaceLayerParams::default(),
            transfer: 1.3e-3,
            gust: 3.0,
            lapse_rate: 4.4e-3,
            convection_time: 2.0 * 3600.0,
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

/// Sea-surface temperature profiles of the Aqua-Planet Experiment (Neale &
/// Hoskins 2000; Blackburn et al. 2013), zonally and hemispherically
/// symmetric, 0 °C poleward of 60°.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ApeSst {
    /// `27 (1 − sin²(3φ/2))` °C: APE's main experiment.
    Control,
    /// The mean of CONTROL and FLAT (`27 (1 − sin⁴(3φ/2))` °C): closest to
    /// the observed Pacific profile.
    Qobs,
}

impl ApeSst {
    /// Sea-surface temperature (K) at latitude `lat` (radians).
    pub fn at(self, lat: f64) -> f64 {
        let celsius = if lat.abs() < std::f64::consts::FRAC_PI_3 {
            let s2 = fms_math::sin(1.5 * lat).powi(2);
            match self {
                Self::Control => 27.0 * (1.0 - s2),
                Self::Qobs => 27.0 * (1.0 - 0.5 * (s2 + s2 * s2)),
            }
        } else {
            0.0
        };
        273.15 + celsius
    }
}

/// Air temperature at the surface (K): the lowest layer's mean temperature
/// extrapolated down the convective lapse rate, given that layer's Exner
/// function `exner` and potential temperature `theta_lowest`.
pub fn surface_air_temperature(
    params: &PhysicsParams,
    gas: Gas,
    gravity: f64,
    ps: f64,
    exner: f64,
    theta_lowest: f64,
) -> f64 {
    let p_mean = REFERENCE_PRESSURE * fms_math::pow(exner, 1.0 / gas.kappa());
    let beta = gas.gas_constant * params.lapse_rate / gravity;
    theta_lowest * exner * fms_math::pow(ps / p_mean, beta)
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

/// Potential temperatures a two-layer column (interface `sigma`) relaxes to
/// if it is less stable than the critical lapse rate, else `None`. The
/// adjusted column keeps its enthalpy `Σ δσₖ Πₖ θₖ`.
pub fn convective_target(
    params: &PhysicsParams,
    gas: Gas,
    sigma: f64,
    gravity: f64,
    ps: f64,
    theta: [f64; 2],
) -> Option<[f64; 2]> {
    let exner = Levels::two(sigma).layer_exner(gas, ps);
    let exner = [exner[0], exner[1]];
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

/// Exchange between the surface, the boundary layer and the layer above in
/// one cell.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
struct SurfaceFlux {
    /// Sensible heat (W m⁻², upward).
    sensible: f64,
    /// Evaporation (kg m⁻² s⁻¹).
    evaporation: f64,
    /// Drag on the bottom layer's wind (s⁻¹).
    drag: f64,
    /// Air entrained across the boundary layer's top (kg m⁻² s⁻¹).
    entrainment: f64,
    /// Shallow cumulus mass flux out of the boundary layer (kg m⁻² s⁻¹).
    cumulus: f64,
    /// Entrainment velocity (m s⁻¹).
    entrainment_velocity: f64,
    /// The boundary layer's convective velocity `w*` (m s⁻¹).
    convective_velocity: f64,
    /// 10 m wind speed (m s⁻¹) and surface stress (N m⁻²).
    wind: f64,
    stress: f64,
    /// Water carried out of the boundary layer by entrainment and shallow
    /// cumulus (kg m⁻² s⁻¹).
    export: f64,
    /// Eastward and northward stress of the surface on the air (N m⁻²).
    stress_east: f64,
    stress_north: f64,
}

/// The radiative–convective forcing of the layered model: radiation, moist
/// physics, a slab surface, surface fluxes and boundary-layer friction. Needs
/// a moist state.
pub struct ColumnPhysics<'a> {
    params: PhysicsParams,
    gas: Gas,
    levels: Levels,
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
    heating: KineticHeating,
    /// Whether the surface temperature is held fixed (an infinite heat
    /// reservoir) instead of evolving as a slab.
    prescribed: bool,
    /// Heat the surface has taken up since the start (J).
    uptake: f64,
    /// Kinetic energy turned to heat by the last tendency evaluation (W),
    /// and since the start (J).
    dissipation: f64,
    dissipated: f64,
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
        model: &LayeredModel<'a>,
        planet: &'a Planet,
        params: PhysicsParams,
        time: f64,
        state: &LayeredState,
        surface_temperature: Vec<f64>,
    ) -> Self {
        assert!(state.is_moist(), "column physics needs a moist state");
        let n = state.ps.len();
        assert_eq!(surface_temperature.len(), n);
        let mut physics = Self {
            params,
            gas: model.gas(),
            levels: model.levels().clone(),
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
            heating: model.kinetic_heating(),
            prescribed: false,
            uptake: 0.0,
            dissipation: 0.0,
            dissipated: 0.0,
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

    /// The physics with the surface temperature held at its starting values
    /// (fixed sea-surface temperature, as models are tuned; e.g. [`ApeSst`]):
    /// the surface becomes an infinite heat reservoir whose uptake the
    /// budget counts as stored heat. The atmosphere settles in about a month
    /// instead of the slab's years.
    pub fn with_prescribed_surface(mut self) -> Self {
        self.prescribed = true;
        self
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
            levels: &self.levels,
            gravity: self.gravity,
        }
    }

    /// One column's view for the moist physics.
    fn moist_column(&self, state: &LayeredState, i: usize) -> MoistColumn {
        let ps = state.ps[i];
        let n = self.levels.count();
        let mut exner = [0.0; MAX_LAYERS];
        self.levels.layer_exner_into(self.gas, ps, &mut exner[..n]);
        let theta: Vec<f64> = (0..n).map(|k| state.theta_ps[k][i] / ps).collect();
        // The stratosphere sits near the skin temperature of the column's
        // outgoing longwave.
        let olr = self.radiation[i].outgoing.max(1.0);
        let humidity: Vec<f64> = (0..n).map(|k| state.q_ps[k][i] / ps).collect();
        // The well-mixed boundary layer reaches the surface.
        let surface_air = theta[n - 1] * self.gas.exner(ps);
        MoistColumn {
            ps,
            temperature: (0..n).map(|k| theta[k] * exner[k]).collect(),
            surface_humidity: humidity[n - 1].clamp(0.0, saturation_humidity(surface_air, ps)),
            humidity,
            surface_air,
            stratosphere: fms_math::pow(olr / (2.0 * STEFAN_BOLTZMANN), 0.25),
        }
    }

    fn update_moisture(&mut self, state: &LayeredState) {
        let moist: Vec<MoistTendency> = (0..state.ps.len())
            .into_par_iter()
            .with_min_len(MIN_CHUNK / 16)
            .map(|i| moist_physics(&self.setting(), &self.moist_column(state, i)))
            .collect();
        self.moist = moist;
    }

    fn update_surface_fluxes(&mut self, state: &LayeredState) {
        let (mesh, geometry) = (self.ops.mesh(), self.ops.geometry());
        let n = self.levels.count();
        let bottom = n - 1;
        let winds = geometry.reconstruct_vectors(mesh, &state.u[bottom]);
        let centers = mesh.cell_centers();
        let (gas, g) = (self.gas, self.gravity);
        let (rd, cp) = (gas.gas_constant, gas.heat_capacity);
        let sl = &self.params.surface_layer;
        let ds = self.levels.dsigma();
        let top = self.levels.interfaces()[bottom];
        let old = &self.surface;
        let surface: Vec<SurfaceFlux> = (0..state.ps.len())
            .into_par_iter()
            .with_min_len(MIN_CHUNK / 16)
            .map(|i| {
                let (e, nn) = east_north(centers[i], winds[i]);
                let ps = state.ps[i];
                let mut exner = [0.0; MAX_LAYERS];
                self.levels.layer_exner_into(gas, ps, &mut exner[..n]);
                let theta = state.theta_ps[bottom][i] / ps;
                let q = (state.q_ps[bottom][i] / ps).max(0.0);
                let surface_exner = gas.exner(ps);
                let air = theta * surface_exner;
                // Heights of the layer's mean level and top (hypsometric).
                let t = theta * exner[bottom];
                let scale = rd * t / g;
                let p_mean = REFERENCE_PRESSURE * fms_math::pow(exner[bottom], 1.0 / gas.kappa());
                let height = scale * fms_math::ln(ps / p_mean);
                let depth = -scale * fms_math::ln(top);
                let ts = self.surface_temperature[i];
                let x = exchange(
                    sl,
                    Surface::Ocean,
                    &SurfaceAir {
                        wind: (e * e + nn * nn).sqrt(),
                        height,
                        theta,
                        humidity: q,
                        surface_temperature: ts,
                        ps,
                        surface_exner,
                        convective_velocity: old[i].convective_velocity,
                        gravity: g,
                    },
                );
                let density = ps / (rd * air);
                let sensible = density * cp * x.heat * x.speed * (ts - air);
                let evaporation =
                    density * x.moisture * x.speed * (saturation_humidity(ts, ps) - q);
                // Clouds' longwave cooling of the layer drives turbulence
                // from its top.
                let cloud_top_cooling =
                    self.clear[i].layers[bottom] - self.radiation[i].layers[bottom];
                let w_star = convective_velocity(
                    g,
                    depth,
                    theta,
                    density,
                    cp,
                    sensible,
                    evaporation,
                    cloud_top_cooling,
                );
                // Entrainment of the layer above, by the jump in virtual θ
                // at the boundary layer's top: the layer above extrapolated
                // down to it along the profile's lapse rate.
                let above = bottom - 1;
                let p_above = REFERENCE_PRESSURE * fms_math::pow(exner[above], 1.0 / gas.kappa());
                let beta = rd * sl.top_lapse_factor * self.params.radiation.lapse_rate / g;
                let t_top = state.theta_ps[above][i] / ps
                    * exner[above]
                    * fms_math::pow(top * ps / p_above, beta);
                let virtual_factor = |k: usize| 1.0 + 0.608 * (state.q_ps[k][i] / ps).max(0.0);
                let below = theta * virtual_factor(bottom);
                let jump = t_top / gas.exner(top * ps) * virtual_factor(above) - below;
                let w_e = entrainment_velocity(sl, g, depth, w_star, below, jump);
                // Shallow cumulus: surface air saturated by the top and still
                // buoyant there.
                let p_top = top * ps;
                let (t_parcel, q_parcel) = lift(gas, air, q, ps, p_top);
                let cumulus = q_parcel < q
                    && t_parcel * (1.0 + 0.608 * q_parcel) > t_top * virtual_factor(above);
                let density_top = p_top / (rd * t_top);
                let momentum = x.momentum(density);
                SurfaceFlux {
                    sensible,
                    evaporation,
                    drag: g * momentum * x.wind_ratio / (ps * ds[bottom]),
                    entrainment: density_top * w_e,
                    cumulus: if cumulus {
                        sl.cumulus * density_top * w_star
                    } else {
                        0.0
                    },
                    entrainment_velocity: w_e,
                    convective_velocity: w_star,
                    wind: x.wind,
                    stress: momentum * x.wind,
                    export: (density_top * w_e
                        + if cumulus {
                            sl.cumulus * density_top * w_star
                        } else {
                            0.0
                        })
                        * (state.q_ps[bottom][i] - state.q_ps[above][i])
                        / ps,
                    stress_east: -momentum * x.wind_ratio * e,
                    stress_north: -momentum * x.wind_ratio * nn,
                }
            })
            .collect();
        self.surface = surface;
    }

    /// Radiation for the interval starting now, with sunlight averaged over
    /// it, all-sky and clear-sky.
    fn radiate(&mut self, state: &LayeredState) {
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
        let (params, gas, gravity) = (&self.params, self.gas, self.gravity);
        let levels = &self.levels;
        let layers = levels.count();
        let clouds = &params.moisture.clouds;
        let moist = &self.moist;
        let both: Vec<(Fluxes, Fluxes)> = (0..n)
            .into_par_iter()
            .with_min_len(MIN_CHUNK / 16)
            .map(|i| {
                let ps = state.ps[i];
                let cover = moist[i].clouds;
                let theta: Vec<f64> = (0..layers).map(|k| state.theta_ps[k][i] / ps).collect();
                let q: Vec<f64> = (0..layers).map(|k| state.q_ps[k][i] / ps).collect();
                let sheets: [Cloud; 2] = std::array::from_fn(|c| Cloud {
                    sigma: clouds.sigma[c],
                    reflectance: cover[c] * clouds.albedo[c],
                    emissivity: cover[c] * clouds.emissivity[c],
                });
                let mut column = Column {
                    ps,
                    surface_temperature: self.surface_temperature[i],
                    theta: &theta,
                    levels,
                    gravity,
                    insolation: insolation[i],
                    co2: params.co2,
                    humidity: Humidity::Layers(&q),
                    clouds: &sheets,
                };
                let all = column_fluxes(&params.radiation, gas, &column);
                column.clouds = &[Cloud::NONE; 2];
                (all, column_fluxes(&params.radiation, gas, &column))
            })
            .collect();
        (self.radiation, self.clear) = both.into_iter().unzip();
        self.next_radiation = self.time + self.params.radiation_interval;
    }
}

impl Forcing for ColumnPhysics<'_> {
    fn add_tendency(&mut self, state: &LayeredState, out: &mut LayeredState) {
        assert!(
            state.is_moist() && out.is_moist(),
            "column physics needs a moist state"
        );
        self.update_moisture(state);
        let (gas, gravity) = (self.gas, self.gravity);
        let levels = &self.levels;
        let n = levels.count();
        let lowest = n - 1;
        let ds = levels.dsigma();
        let mut exner = vec![0.0; n * state.ps.len()];
        exner
            .par_chunks_mut(n)
            .with_min_len(MIN_CHUNK / 16)
            .zip(state.ps.par_iter())
            .for_each(|(e, &ps)| levels.layer_exner_into(gas, ps, e));
        let (radiation, surface, moist) = (&self.radiation, &self.surface, &self.moist);
        let above = lowest - 1;
        let cp = gas.heat_capacity;
        let top = levels.interfaces()[lowest];
        for (k, (dt, dq)) in out.theta_ps.iter_mut().zip(&mut out.q_ps).enumerate() {
            let exner = &exner;
            dt.par_iter_mut()
                .with_min_len(MIN_CHUNK)
                .zip(dq.par_iter_mut())
                .enumerate()
                .for_each(|(i, (dt, dq))| {
                    let ps = state.ps[i];
                    let e = exner[i * n + k];
                    let m = &moist[i];
                    let s = &surface[i];
                    let mut heating = radiation[i].layers[k];
                    let mut water = ps * m.humidity[k];
                    if k == lowest || k == above {
                        // Entrainment: dry static energy `c_p Π θ + gz`
                        // and water mixed across the boundary layer's top.
                        let theta = |k: usize| state.theta_ps[k][i] / ps;
                        let q = |k: usize| state.q_ps[k][i] / ps;
                        let exner_top = gas.exner(top * ps);
                        let heat = s.entrainment * cp * exner_top * (theta(above) - theta(lowest));
                        let moisture =
                            gravity * (s.entrainment + s.cumulus) * (q(above) - q(lowest));
                        if k == lowest {
                            heating += s.sensible + heat;
                            water += gravity * s.evaporation / ds[k] + moisture / ds[k];
                        } else {
                            heating -= heat;
                            water -= moisture / ds[k];
                        }
                    }
                    *dt +=
                        theta_tendency(gas, gravity, ds[k], e, heating) + ps * m.temperature[k] / e;
                    *dq += water;
                });
        }
        // Surface drag on the boundary layer's wind and entrainment of the
        // layer above's momentum; the kinetic energy they remove heats the
        // boundary layer.
        let edges = self.ops.mesh().edges();
        let (u_low, u_up) = (&state.u[lowest], &state.u[above]);
        let mut du_low = vec![0.0; edges.len()];
        let mut du_up = vec![0.0; edges.len()];
        du_low
            .par_iter_mut()
            .with_min_len(MIN_CHUNK)
            .zip(du_up.par_iter_mut())
            .enumerate()
            .for_each(|(e, (dl, du))| {
                let [a, b] = edges[e].cells.map(|c| c as usize);
                let drag = 0.5 * (surface[a].drag + surface[b].drag);
                let exchange = |s: &SurfaceFlux| s.entrainment + s.cumulus;
                let mixing = 0.5 * (exchange(&surface[a]) + exchange(&surface[b])) * gravity
                    / (0.5 * (state.ps[a] + state.ps[b]));
                let shear = u_up[e] - u_low[e];
                *dl = -drag * u_low[e] + mixing * shear / ds[lowest];
                *du = -mixing * shear / ds[above];
            });
        let heating = &mut self.heating;
        heating.prepare(&state.ps);
        let ops = &self.ops;
        heating.add(
            ops,
            &state.ps,
            lowest,
            lowest,
            u_low,
            &du_low,
            &mut out.theta_ps[lowest],
        );
        heating.add(
            ops,
            &state.ps,
            above,
            lowest,
            u_up,
            &du_up,
            &mut out.theta_ps[lowest],
        );
        self.dissipation = -(heating.rate(ops, &state.ps, lowest, gravity, u_low, &du_low)
            + heating.rate(ops, &state.ps, above, gravity, u_up, &du_up));
        for (k, du) in [(lowest, &du_low), (above, &du_up)] {
            out.u[k]
                .par_iter_mut()
                .with_min_len(MIN_CHUNK)
                .zip(du.par_iter())
                .for_each(|(o, &d)| *o += d);
        }
    }

    fn advance(&mut self, state: &LayeredState, dt: f64) {
        // The slab and the budgets use the same stored fluxes and
        // precipitation the atmosphere was given.
        let areas = &self.ops.geometry().cell_areas;
        let r2 = self.ops.radius() * self.ops.radius();
        let capacity = self.params.heat_capacity();
        let (radiation, surface, moist) = (&self.radiation, &self.surface, &self.moist);
        let net = |i: usize| {
            let s = surface[i];
            radiation[i].surface()
                - s.sensible
                - LATENT_HEAT * s.evaporation
                - FUSION_HEAT * moist[i].snow
        };
        if self.prescribed {
            let uptake: f64 = (0..areas.len()).map(|i| areas[i] * net(i)).sum();
            self.uptake += dt * r2 * uptake;
        } else {
            self.surface_temperature
                .par_iter_mut()
                .with_min_len(MIN_CHUNK)
                .enumerate()
                .for_each(|(i, t)| *t += dt * net(i) / capacity);
        }
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
        self.dissipated += dt * self.dissipation;
        self.time += dt;
        self.update_surface_fluxes(state);
        if self.time >= self.next_radiation - 1e-6 {
            self.radiate(state);
        }
    }

    fn budget(&self) -> Option<ForcingBudget> {
        let areas = &self.ops.geometry().cell_areas;
        let r2 = self.ops.radius() * self.ops.radius();
        let stored: f64 = if self.prescribed {
            self.uptake
        } else {
            self.surface_temperature
                .iter()
                .zip(areas)
                .map(|(t, a)| a * r2 * self.params.heat_capacity() * t)
                .sum()
        };
        Some(ForcingBudget {
            absorbed: self.absorbed,
            outgoing: self.outgoing,
            stored,
            evaporated: self.evaporated,
            precipitated: self.precipitated,
            dissipated: self.dissipated,
        })
    }

    fn cell_fields(&self) -> Vec<CellField> {
        let rad = |f: fn(&Fluxes) -> f64| self.radiation.iter().map(f).collect();
        let lowest = self.levels.count() - 1;
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
                self.radiation.iter().map(Fluxes::total_water).collect(),
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
                "atm.boundary_layer_moist_heating",
                "Moist-physics heating of the boundary layer",
                "K day-1",
                self.moist
                    .iter()
                    .map(|m| m.temperature[lowest] * DAY)
                    .collect(),
            ),
            field(
                "atm.surface_wind",
                "10 m wind speed",
                "m s-1",
                self.surface.iter().map(|s| s.wind).collect(),
            ),
            field(
                "atm.surface_stress",
                "Surface wind stress",
                "N m-2",
                self.surface.iter().map(|s| s.stress).collect(),
            ),
            field(
                "atm.shallow_cumulus",
                "Shallow cumulus mass flux out of the boundary layer",
                "g m-2 s-1",
                self.surface.iter().map(|s| 1e3 * s.cumulus).collect(),
            ),
            field(
                "moist.condensation_boundary",
                "Large-scale condensation in the boundary layer",
                "mm day-1",
                self.moist
                    .iter()
                    .map(|m| m.condensed[lowest] * DAY)
                    .collect(),
            ),
            field(
                "moist.condensation_above",
                "Large-scale condensation in the layer above the boundary layer",
                "mm day-1",
                self.moist
                    .iter()
                    .map(|m| m.condensed[lowest - 1] * DAY)
                    .collect(),
            ),
            field(
                "moist.condensation_upper",
                "Large-scale condensation in the higher layers",
                "mm day-1",
                self.moist
                    .iter()
                    .map(|m| m.condensed[..lowest - 1].iter().sum::<f64>() * DAY)
                    .collect(),
            ),
            field(
                "moist.reevaporation_boundary",
                "Rain re-evaporated in the boundary layer",
                "mm day-1",
                self.moist
                    .iter()
                    .map(|m| m.reevaporated[lowest] * DAY)
                    .collect(),
            ),
            field(
                "moist.reevaporation_above",
                "Rain re-evaporated above the boundary layer",
                "mm day-1",
                self.moist
                    .iter()
                    .map(|m| m.reevaporated[..lowest].iter().sum::<f64>() * DAY)
                    .collect(),
            ),
            field(
                "moist.convective_boundary",
                "Water convection adds to the boundary layer",
                "mm day-1",
                self.moist
                    .iter()
                    .map(|m| m.convective_boundary * DAY)
                    .collect(),
            ),
            field(
                "atm.boundary_layer_export",
                "Water carried out of the boundary layer by entrainment and cumulus",
                "mm day-1",
                self.surface.iter().map(|s| s.export * DAY).collect(),
            ),
            field(
                "atm.surface_stress_east",
                "Eastward surface stress on the air",
                "N m-2",
                self.surface.iter().map(|s| s.stress_east).collect(),
            ),
            field(
                "atm.surface_stress_north",
                "Northward surface stress on the air",
                "N m-2",
                self.surface.iter().map(|s| s.stress_north).collect(),
            ),
            field(
                "rad.free_troposphere_heating",
                "Radiative heating of the layers above the boundary layer",
                "W m-2",
                self.radiation
                    .iter()
                    .map(|f| f.layers[..lowest].iter().sum())
                    .collect(),
            ),
            field(
                "rad.clear_free_troposphere_heating",
                "Clear-sky radiative heating of the layers above the boundary layer",
                "W m-2",
                self.clear
                    .iter()
                    .map(|f| f.layers[..lowest].iter().sum())
                    .collect(),
            ),
            field(
                "atm.convective_velocity",
                "Boundary-layer convective velocity w*",
                "m s-1",
                self.surface.iter().map(|s| s.convective_velocity).collect(),
            ),
            field(
                "atm.entrainment",
                "Entrainment velocity at the boundary-layer top",
                "mm s-1",
                self.surface
                    .iter()
                    .map(|s| 1e3 * s.entrainment_velocity)
                    .collect(),
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
                "Top-layer relative humidity",
                "1",
                moist(|m| m.relative_humidity[0]),
            ),
            field(
                "atm.relative_humidity_lower",
                "Lowest-layer relative humidity",
                "1",
                self.moist
                    .iter()
                    .map(|m| m.relative_humidity[lowest])
                    .collect(),
            ),
        ]
    }
}

/// A moist start for the 3-D model on `planet`: everywhere the single
/// column's equilibrium (with `params`' CO₂) at rest, with its fixed-profile
/// humidity, plus 0.1 K of seeded noise in θ to break the symmetry. Returns
/// the state and the surface temperature (K).
///
/// The single column has two equal-mass layers; other layers take the mean
/// temperature of its profile over their range (each layer's lapse rate
/// through its mean, isothermal above the upper layer's mean level).
///
/// Over terrain (the model's surface geopotential) the same temperatures
/// hold on pressure surfaces, so the column is hydrostatic and at rest: the
/// planet's surface pressure is the sea-level value, and each cell's follows
/// from its height.
pub fn initial_state(
    model: &LayeredModel,
    planet: &Planet,
    params: &PhysicsParams,
) -> (LayeredState, f64) {
    let mut column = SingleColumn::earth(planet.mean_flux());
    column.params.co2 = params.co2;
    column.ps = planet.params().surface_pressure;
    column.gravity = model.gravity();
    let eq = column.equilibrium();
    let (gas, levels, ps) = (model.gas(), model.levels(), column.ps);
    let exner = levels.layer_exner(gas, ps);
    let two = Levels::two(column.sigma);
    let temperature: Vec<f64> = if *levels == two {
        eq.temperature.to_vec()
    } else {
        let profile = Profile::new(
            &params.radiation,
            gas,
            &two,
            ps,
            column.gravity,
            &eq.temperature,
        );
        let cap = two.mean_pressures(gas, ps)[0];
        (0..levels.count())
            .map(|k| {
                let (a, b) = levels.bounds(k);
                GAUSS
                    .iter()
                    .map(|&(x, w)| {
                        let s = 0.5 * (a + b) + 0.5 * (b - a) * x;
                        let t = if s * ps < cap {
                            eq.temperature[0]
                        } else {
                            profile.at(s)
                        };
                        0.5 * w * t
                    })
                    .sum()
            })
            .collect()
    };
    let theta0: Vec<f64> = temperature.iter().zip(&exner).map(|(t, e)| t / e).collect();
    let q = fixed_humidity(params, gas, levels, model.gravity(), ps, &theta0);
    let mesh = model.operators().mesh();
    let (nc, ne) = (mesh.n_cells(), mesh.n_edges());
    let n = levels.count();
    let reference = RestProfile::new(
        gas,
        ps,
        theta0[n - 1] * gas.exner(ps),
        &levels.mean_pressures(gas, ps),
        &temperature,
    );
    // Surface pressure, θ and humidity of each column at rest.
    let columns: Vec<(f64, Vec<f64>, Vec<f64>)> = model
        .surface_geopotential()
        .par_iter()
        .with_min_len(MIN_CHUNK / 16)
        .map(|&phi| {
            if phi == 0.0 {
                return (ps, theta0.clone(), q.clone());
            }
            let p = reference.pressure_at(phi);
            let theta: Vec<f64> = levels
                .mean_pressures(gas, p)
                .iter()
                .zip(levels.layer_exner(gas, p))
                .map(|(&pk, e)| reference.temperature(pk) / e)
                .collect();
            let q = fixed_humidity(params, gas, levels, model.gravity(), p, &theta);
            (p, theta, q)
        })
        .collect();
    let mut x = 0x2545_F491_4F6C_DD1Du64;
    let theta: Vec<Vec<f64>> = (0..n)
        .map(|k| {
            columns
                .iter()
                .map(|c| {
                    x ^= x << 13;
                    x ^= x >> 7;
                    x ^= x << 17;
                    c.1[k] + 0.1 * ((x >> 11) as f64 / (1u64 << 52) as f64 - 1.0)
                })
                .collect()
        })
        .collect();
    let humidity = (0..n)
        .map(|k| columns.iter().map(|c| c.2[k]).collect())
        .collect();
    let surface_pressure = columns.iter().map(|c| c.0).collect();
    debug_assert_eq!(columns.len(), nc);
    let state =
        LayeredState::new(surface_pressure, vec![vec![0.0; ne]; n], theta).with_humidity(humidity);
    (state, eq.surface_temperature)
}

/// Starting surface temperature per cell (K): `sea_level` (K) at sea level,
/// falling with height at the radiation's lapse rate.
pub fn initial_surface_temperature(
    model: &LayeredModel,
    params: &PhysicsParams,
    sea_level: f64,
) -> Vec<f64> {
    let lapse = params.radiation.lapse_rate / model.gravity();
    model
        .surface_geopotential()
        .iter()
        .map(|phi| sea_level - lapse * phi)
        .collect()
}

/// A horizontally uniform temperature profile: linear in ln p between the
/// surface at sea level and each layer's mean level, isothermal above the
/// top one, and the hydrostatic pressure at a given height in it.
struct RestProfile {
    gas_constant: f64,
    /// `(ln p, T)` from sea level upward.
    nodes: Vec<(f64, f64)>,
}

impl RestProfile {
    /// The profile with temperature `surface` (K) at sea-level pressure `ps`
    /// (Pa), and `temperature[k]` at pressures `pressure[k]`, layers ordered
    /// top down.
    fn new(gas: Gas, ps: f64, surface: f64, pressure: &[f64], temperature: &[f64]) -> Self {
        let mut nodes = vec![(fms_math::ln(ps), surface)];
        nodes.extend(
            pressure
                .iter()
                .zip(temperature)
                .rev()
                .map(|(&p, &t)| (fms_math::ln(p), t)),
        );
        Self {
            gas_constant: gas.gas_constant,
            nodes,
        }
    }

    /// Temperature (K) at pressure `p` (Pa).
    fn temperature(&self, p: f64) -> f64 {
        let x = fms_math::ln(p);
        for pair in self.nodes.windows(2) {
            let ((xa, ta), (xb, tb)) = (pair[0], pair[1]);
            if x >= xb {
                return ta + (tb - ta) * (xa - x) / (xa - xb);
            }
        }
        self.nodes[self.nodes.len() - 1].1
    }

    /// Pressure (Pa) at geopotential `phi` (m² s⁻², at or above sea level).
    ///
    /// Within a segment `T = T_a + s d` with `d = ln p_a − ln p`, so the
    /// hypsometric equation `c = (Φ − Φ_a)/R = T_a d + s d²/2` has the
    /// root `d = 2c / (T_a + √(T_a² + 2 s c))`, regular at `s = 0`.
    fn pressure_at(&self, phi: f64) -> f64 {
        let mut c = phi / self.gas_constant;
        for pair in self.nodes.windows(2) {
            let ((xa, ta), (xb, tb)) = (pair[0], pair[1]);
            let span = xa - xb;
            let thickness = 0.5 * (ta + tb) * span;
            if c <= thickness {
                let s = (tb - ta) / span;
                let d = 2.0 * c / (ta + (ta * ta + 2.0 * s * c).sqrt());
                return fms_math::exp(xa - d);
            }
            c -= thickness;
        }
        let (x, t) = self.nodes[self.nodes.len() - 1];
        fms_math::exp(x - c / t)
    }
}

/// Gauss–Legendre nodes on [−1, 1].
const GAUSS: [(f64, f64); 4] = [
    (-0.861_136_311_594_052_6, 0.347_854_845_137_453_8),
    (-0.339_981_043_584_856_3, 0.652_145_154_862_546_1),
    (0.339_981_043_584_856_3, 0.652_145_154_862_546_1),
    (0.861_136_311_594_052_6, 0.347_854_845_137_453_8),
];

/// Initial humidity (kg kg⁻¹, per layer) for a column at fixed relative
/// humidity: the Manabe–Wetherald profile of the single column.
pub fn fixed_humidity(
    params: &PhysicsParams,
    gas: Gas,
    levels: &Levels,
    gravity: f64,
    ps: f64,
    theta: &[f64],
) -> Vec<f64> {
    let f = column_fluxes(
        &params.radiation,
        gas,
        &Column {
            ps,
            surface_temperature: 288.0,
            theta,
            levels,
            gravity,
            insolation: 0.0,
            co2: params.co2,
            humidity: Humidity::Fixed,
            clouds: &[Cloud::NONE; 2],
        },
    );
    let ds = levels.dsigma();
    (0..levels.count())
        .map(|k| f.water[k] * gravity / (ds[k] * ps))
        .collect()
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
#[derive(Clone, Debug, PartialEq)]
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

    /// The column's two layers.
    pub fn levels(&self) -> Levels {
        Levels::two(self.sigma)
    }

    /// Fluxes of the current state.
    pub fn fluxes(&self) -> Fluxes {
        column_fluxes(
            &self.params.radiation,
            self.gas,
            &Column {
                ps: self.ps,
                surface_temperature: self.surface_temperature,
                theta: &self.theta,
                levels: &self.levels(),
                gravity: self.gravity,
                insolation: self.insolation,
                co2: self.params.co2,
                humidity: Humidity::Fixed,
                clouds: &[
                    Cloud::deck(self.sigma, self.params.radiation.cloud_albedo),
                    Cloud::NONE,
                ],
            },
        )
    }

    fn layer_exner(&self) -> [f64; 2] {
        let e = self.levels().layer_exner(self.gas, self.ps);
        [e[0], e[1]]
    }

    /// Sensible and latent heat flux (W m⁻²).
    pub fn surface_fluxes(&self) -> (f64, f64) {
        let air = surface_air_temperature(
            &self.params,
            self.gas,
            self.gravity,
            self.ps,
            self.layer_exner()[1],
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
        let exner = self.layer_exner();
        self.theta = [t[0] / exner[0], t[1] / exner[1]];
    }

    /// One step of `dt` seconds (forward Euler, radiation every step).
    pub fn step(&mut self, dt: f64) {
        let f = self.fluxes();
        let (sensible, latent) = self.surface_fluxes();
        let exner = self.layer_exner();
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
        let exner = self.layer_exner();
        Equilibrium {
            fluxes: self.fluxes(),
            surface_temperature: self.surface_temperature,
            temperature: [self.theta[0] * exner[0], self.theta[1] * exner[1]],
            surface_fluxes: self.surface_fluxes(),
            steps,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const GAS: Gas = Gas::EARTH_AIR;

    #[test]
    fn an_isothermal_rest_profile_has_the_barometric_pressure() {
        let t = 260.0;
        let profile = RestProfile::new(GAS, 1.0e5, t, &[3.0e4, 7.0e4, 9.0e4], &[t; 3]);
        for z in [100.0, 1500.0, 6000.0, 12_000.0] {
            let phi = 9.81 * z;
            let expected = 1.0e5 * fms_math::exp(-phi / (GAS.gas_constant * t));
            let p = profile.pressure_at(phi);
            assert!(
                (p / expected - 1.0).abs() < 1e-12,
                "z {z}: {p} vs {expected}"
            );
        }
    }

    #[test]
    fn pressure_at_inverts_the_hypsometric_integral() {
        let (p, t) = ([3.0e4, 7.0e4, 9.0e4], [225.0, 265.0, 282.0]);
        let profile = RestProfile::new(GAS, 1.0e5, 290.0, &p, &t);
        // Layer thicknesses by the trapezoid in ln p, in segments of 1 Pa.
        for phi in [500.0, 8_000.0, 30_000.0, 60_000.0] {
            let target = profile.pressure_at(phi);
            let mut integral = 0.0;
            let mut a = 1.0e5;
            while a > target {
                let b = (a - 1.0_f64).max(target);
                let mid = 0.5 * (profile.temperature(a) + profile.temperature(b));
                integral += GAS.gas_constant * mid * fms_math::ln(a / b);
                a = b;
            }
            assert!((integral / phi - 1.0).abs() < 1e-4, "{integral} vs {phi}");
        }
    }

    #[test]
    fn the_profile_passes_through_its_nodes() {
        let (p, t) = ([3.0e4, 7.0e4, 9.0e4], [225.0, 265.0, 282.0]);
        let profile = RestProfile::new(GAS, 1.0e5, 290.0, &p, &t);
        for (p, t) in p.iter().zip(t) {
            assert!((profile.temperature(*p) - t).abs() < 1e-9);
        }
        assert_eq!(profile.temperature(1.0e3), 225.0);
    }
}
