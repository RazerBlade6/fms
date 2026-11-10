//! Sea ice after Codron (2012) and the code of Bhatnagar et al. (2026): each
//! cell has one fraction, a mean thickness and a surface temperature, with a
//! linear profile down to the freezing point at the base.
//!
//! The temperature and the heat for growth and melting change every step
//! ([`SeaIce::step`]); thickness and area change only in
//! [`SeaIce::rearrange`], which the atmosphere calls with its radiation, so
//! that the areas its stored sunlight was divided over stay the areas that
//! receive it.
//!
//! Heat is measured from water at the freezing point, so ice holds
//! `−e(T) h` per unit area with `e` from [`IceParams::melting_energy`], and
//! every path below moves heat by that same `e`. New ice takes the
//! temperature of the ice already there, as in the authors' code.

/// How sea ice forms, melts and looks; [`Default`] is the authors'.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IceParams {
    /// Freezing point of seawater (K).
    pub freezing: f64,
    /// Melting point of the ice's surface (K).
    pub melting: f64,
    /// Density (kg m⁻³).
    pub density: f64,
    /// Specific heat (J kg⁻¹ K⁻¹).
    pub specific_heat: f64,
    /// Thermal conductivity (W m⁻¹ K⁻¹).
    pub conductivity: f64,
    /// Latent heat of fusion (J kg⁻¹).
    pub fusion: f64,
    /// Ice thinner than this (m) is removed.
    pub min_thickness: f64,
    /// Thickness new ice starts from, at least (m).
    pub thin: f64,
    /// Thickness growth reaches before it spreads, and melting before it
    /// shrinks area (m).
    pub thick: f64,
    /// Thickness new ice starts from, at most (m).
    pub new: f64,
    /// The thickest ice can be (m).
    pub max_thickness: f64,
    /// Less cover than this is removed.
    pub min_fraction: f64,
    /// The most cover a cell can have.
    pub max_fraction: f64,
    /// Albedo of thick bare ice in each band of sunlight. The paper's
    /// visible 0.65 is for a narrower band; 0.57 in ours keeps its broadband
    /// 0.48.
    pub albedo_limit: [f64; 2],
    /// Thickness over which albedo rises from the water's to its limit (m).
    pub albedo_depth: f64,
    /// Momentum roughness of the surface (m).
    pub roughness: f64,
    /// Bulk coefficient of the ocean-to-ice heat flux.
    pub ocean_exchange: f64,
    /// Smallest friction velocity under the ice (m s⁻¹).
    pub min_friction_velocity: f64,
}

impl Default for IceParams {
    fn default() -> Self {
        Self {
            freezing: 271.35,
            melting: 273.15,
            density: 917.0,
            specific_heat: 2067.0,
            conductivity: 2.17,
            fusion: 3.34e5,
            min_thickness: 0.01,
            thin: 0.15,
            thick: 1.5,
            new: 1.0,
            max_thickness: 10.0,
            min_fraction: 0.001,
            max_fraction: 0.999_999,
            albedo_limit: [0.57, 0.31],
            albedo_depth: 0.3,
            roughness: 1e-3,
            ocean_exchange: 0.0057,
            min_friction_velocity: 5e-4,
        }
    }
}

impl IceParams {
    /// Energy to turn a cubic metre of ice with surface temperature `t` (K)
    /// into water at the freezing point (J m⁻³).
    pub fn melting_energy(&self, t: f64) -> f64 {
        self.density * (self.fusion + 0.5 * self.specific_heat * (self.freezing - t))
    }

    /// Heat capacity of the ice's temperature, with its linear profile, per
    /// unit area of ice `h` metres thick (J m⁻² K⁻¹).
    pub fn capacity(&self, h: f64) -> f64 {
        0.5 * self.density * self.specific_heat * h
    }
}

/// What the atmosphere does to one cell's ice over a step (per unit area of
/// ice).
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct IceForcing {
    /// Net heat into the surface (W m⁻²).
    pub net: f64,
    /// The surface temperature `net` was computed at (K).
    pub skin: f64,
    /// Surface wind stress (N m⁻²).
    pub stress: f64,
    /// Part of `net`'s loss that pays to thaw what sublimates (W m⁻²); it
    /// comes back as ice removed from the top.
    pub sublimation: f64,
}

/// The ice of every cell.
#[derive(Clone, Debug)]
pub(crate) struct SeaIce {
    pub(crate) params: IceParams,
    /// Share of each cell's ocean covered.
    pub(crate) fraction: Vec<f64>,
    /// Mean thickness where covered (m).
    pub(crate) thickness: Vec<f64>,
    /// Surface temperature (K).
    pub(crate) temperature: Vec<f64>,
    /// Heat gained at the base since the last rearrangement (J per m² of
    /// ice); negative grows ice.
    pub(crate) bottom: Vec<f64>,
    /// Heat that has reached the top for melting since then (J per m² of
    /// ice).
    top: Vec<f64>,
    /// Heat the ice limit has removed from the record of ice cut off (J per
    /// m² of ocean).
    pub(crate) limited: Vec<f64>,
}

/// What [`SeaIce::step`] needs of the ocean under the ice.
pub(crate) struct Water<'a> {
    pub(crate) mixed: &'a mut [f64],
    pub(crate) area: &'a [f64],
    /// Heat capacity of the mixed layer (J m⁻² K⁻¹).
    pub(crate) capacity: f64,
    /// Density (kg m⁻³) and specific heat (J kg⁻¹ K⁻¹) of seawater.
    pub(crate) density: f64,
    pub(crate) specific_heat: f64,
}

impl SeaIce {
    pub(crate) fn new(params: IceParams, n: usize) -> Self {
        Self {
            params,
            fraction: vec![0.0; n],
            thickness: vec![0.0; n],
            temperature: vec![params.freezing; n],
            bottom: vec![0.0; n],
            top: vec![0.0; n],
            limited: vec![0.0; n],
        }
    }

    /// Heat held by cell `i`'s ice (J per m² of ocean), from water at the
    /// freezing point.
    pub(crate) fn heat(&self, i: usize) -> f64 {
        let (f, h) = (self.fraction[i], self.thickness[i]);
        f * (self.bottom[i] + self.top[i] - self.params.melting_energy(self.temperature[i]) * h)
            - self.limited[i]
    }

    /// Advance the temperatures `dt` seconds, take the ocean's heat, and
    /// record what is left for melting and growth.
    pub(crate) fn step(
        &mut self,
        dt: f64,
        water: &mut Water,
        forcing: impl Fn(usize) -> IceForcing,
    ) {
        let p = self.params;
        for i in 0..self.fraction.len() {
            let (f, h) = (self.fraction[i], self.thickness[i]);
            if f <= 0.0 || water.area[i] <= 0.0 {
                continue;
            }
            let IceForcing {
                net,
                skin,
                stress,
                sublimation,
            } = forcing(i);
            // The ocean gives up heat only from above freezing.
            let excess = water.mixed[i] - p.freezing;
            let friction = ((1.0 - f) * stress / water.density)
                .sqrt()
                .max(p.min_friction_velocity);
            let turbulent =
                (p.ocean_exchange * water.density * water.specific_heat * friction * excess)
                    .min(water.capacity * excess / dt)
                    .max(0.0);
            water.mixed[i] -= dt * f * turbulent / water.capacity;
            let conducted = p.conductivity * (p.freezing - skin) / h;
            let capacity = p.capacity(h);
            let mut t = self.temperature[i] + dt * (net + conducted) / capacity;
            if t > p.melting {
                self.top[i] += capacity * (t - p.melting);
                t = p.melting;
            }
            self.temperature[i] = t;
            self.top[i] += dt * sublimation;
            self.bottom[i] += dt * (turbulent - conducted);
        }
    }

    /// Spend the heat recorded since the last call on growth and melting,
    /// freeze open water that has fallen below freezing, and remove what is
    /// too little to keep.
    pub(crate) fn rearrange(&mut self, mixed: &mut [f64], area: &[f64], capacity: f64) {
        let p = self.params;
        for i in 0..self.fraction.len() {
            if area[i] <= 0.0 {
                continue;
            }
            let (start, mut f, mut h) = (self.fraction[i], self.fraction[i], self.thickness[i]);
            let mut t = self.temperature[i];
            let e = p.melting_energy(t);
            let (bottom, top) = (self.bottom[i], self.top[i]);
            self.bottom[i] = 0.0;
            self.top[i] = 0.0;
            // Volumes are per unit area of ocean.
            let mut to_ocean = 0.0;
            if f > 0.0 {
                if bottom < 0.0 {
                    let mut v = -bottom * f / e;
                    let spread = ((h - p.thick) / (p.max_thickness - p.thick)).clamp(0.0, 1.0);
                    let dh = (p.thick - h).min(v / f).max(0.0);
                    h += dh;
                    v -= f * dh;
                    if v > 0.0 {
                        let df = (p.max_fraction - f).min(v * spread / h).max(0.0);
                        f += df;
                        v -= df * h;
                        let added = v.min(((p.max_thickness - h) * f).max(0.0));
                        h += added / f;
                        self.limited[i] += e * (v - added);
                    }
                } else {
                    let mut v = bottom * f / e;
                    let lateral = ((h - p.thin) / (p.thick - p.thin)).clamp(0.0, 1.0);
                    let dh = (h - p.thick).min(v / f).max(0.0);
                    h -= dh;
                    v -= f * dh;
                    if v > 0.0 {
                        let dh = (h - p.min_thickness).min(v * lateral / f).max(0.0);
                        h -= dh;
                        v -= f * dh;
                        let df = f.min(v / h);
                        f -= df;
                        v -= df * h;
                        to_ocean += e * v.max(0.0);
                    }
                }
                if f > 0.0 {
                    let mut v = top * start / e;
                    let dh = (h - p.min_thickness).min(v / f).max(0.0);
                    h -= dh;
                    v -= f * dh;
                    if v > 0.0 {
                        let df = f.min(v / h);
                        f -= df;
                        v -= df * h;
                        to_ocean += e * v.max(0.0);
                    }
                } else {
                    to_ocean += top * start;
                }
            }
            mixed[i] += to_ocean / capacity;
            if mixed[i] < p.freezing {
                if f <= 0.0 {
                    t = p.freezing;
                }
                let formed = p.melting_energy(t);
                let v = capacity * (p.freezing - mixed[i]) / formed;
                let new = p.thin.max(p.new.min(h));
                let df = (p.max_fraction - f)
                    .min(p.min_fraction.max(v / new))
                    .max(0.0);
                let volume = f * h + v;
                f += df;
                h = volume / f;
                if h > p.max_thickness {
                    self.limited[i] += formed * (h - p.max_thickness) * f;
                    h = p.max_thickness;
                }
                mixed[i] = p.freezing;
            }
            if f > 0.0 && (f < p.min_fraction || h < p.min_thickness) {
                mixed[i] -= p.melting_energy(t) * h * f / capacity;
                f = 0.0;
                h = 0.0;
            }
            if f <= 0.0 {
                f = 0.0;
                h = 0.0;
                t = p.freezing;
            }
            self.fraction[i] = f;
            self.thickness[i] = h;
            self.temperature[i] = t;
        }
    }
}
