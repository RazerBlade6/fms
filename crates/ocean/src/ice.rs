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
//!
//! Snow lies on the ice as a mass per unit area of ice. Above
//! [`IceParams::snow_critical`] it is its own layer, with a surface
//! temperature, over an interface temperature that is the top of the ice;
//! below that the two temperatures are one. Each layer is linear between
//! its two boundary temperatures.

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
    /// Density of snow (kg m⁻³).
    pub snow_density: f64,
    /// Thermal conductivity of snow (W m⁻¹ K⁻¹).
    pub snow_conductivity: f64,
    /// Snow above this (kg m⁻²) is a separate layer, and covers the ice's
    /// albedo completely.
    pub snow_critical: f64,
    /// Albedo of snow, in both bands.
    pub snow_albedo: f64,
    /// Share of the snow landing on a cell's ice that the wind blows into
    /// the water, at no ice cover; it falls with the cover.
    pub snow_blown: f64,
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
    /// Density of the seawater ice floats in (kg m⁻³).
    pub sea_density: f64,
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
            snow_density: 300.0,
            snow_conductivity: 0.31,
            snow_critical: 0.05 * 300.0,
            snow_albedo: 0.65,
            snow_blown: 0.4,
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
            sea_density: 1026.0,
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

    /// Energy to turn a kilogram of snow with mean temperature `mean` (K)
    /// into water at the freezing point (J kg⁻¹).
    pub fn snow_melting_energy(&self, mean: f64) -> f64 {
        self.fusion + self.specific_heat * (self.freezing - mean)
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
    /// comes back as snow and ice removed from the top.
    pub sublimation: f64,
    /// Snow landing (kg m⁻² s⁻¹).
    pub snowfall: f64,
}

/// The ice of every cell.
#[derive(Clone, Debug)]
pub(crate) struct SeaIce {
    pub(crate) params: IceParams,
    /// Share of each cell's ocean covered.
    pub(crate) fraction: Vec<f64>,
    /// Mean thickness where covered (m).
    pub(crate) thickness: Vec<f64>,
    /// Temperature of the surface, snow or ice (K).
    pub(crate) temperature: Vec<f64>,
    /// Temperature of the top of the ice (K).
    pub(crate) interface: Vec<f64>,
    /// Snow on the ice (kg per m² of ice).
    pub(crate) snow: Vec<f64>,
    /// Snow that has landed since the last rearrangement (kg per m² of ice).
    pub(crate) fallen: Vec<f64>,
    /// Heat gained at the base since the last rearrangement (J per m² of
    /// ice); negative grows ice.
    pub(crate) bottom: Vec<f64>,
    /// Heat that has reached the top for melting since then (J per m² of
    /// ice).
    pub(crate) top: Vec<f64>,
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

/// The conductances and capacities of one cell's column, per unit area of
/// ice.
struct Column {
    /// Conductance of the ice (W m⁻² K⁻¹).
    ice: f64,
    /// Conductance of the snow layer, if there is one (W m⁻² K⁻¹).
    snow: Option<f64>,
    /// Heat capacity of the interface temperature (J m⁻² K⁻¹).
    interface: f64,
    /// Heat capacity of the surface temperature of a snow layer (J m⁻² K⁻¹).
    surface: f64,
}

impl SeaIce {
    pub(crate) fn new(params: IceParams, n: usize) -> Self {
        Self {
            params,
            fraction: vec![0.0; n],
            thickness: vec![0.0; n],
            temperature: vec![params.freezing; n],
            interface: vec![params.freezing; n],
            snow: vec![0.0; n],
            fallen: vec![0.0; n],
            bottom: vec![0.0; n],
            top: vec![0.0; n],
            limited: vec![0.0; n],
        }
    }

    fn column(&self, h: f64, snow: f64) -> Column {
        let p = &self.params;
        let ice = p.capacity(h);
        let half_snow = 0.5 * p.specific_heat * snow;
        if snow > p.snow_critical {
            Column {
                ice: p.conductivity / h,
                snow: Some(p.snow_conductivity * p.snow_density / snow),
                interface: ice + half_snow,
                surface: half_snow,
            }
        } else {
            Column {
                ice: p.conductivity / h,
                snow: None,
                interface: ice + 2.0 * half_snow,
                surface: 0.0,
            }
        }
    }

    /// Heat held by cell `i`'s ice and snow (J per m² of ocean), from water
    /// at the freezing point.
    pub(crate) fn heat(&self, i: usize) -> f64 {
        let p = &self.params;
        let (f, h) = (self.fraction[i], self.thickness[i]);
        let mean = 0.5 * (self.temperature[i] + self.interface[i]);
        f * (self.bottom[i] + self.top[i]
            - p.fusion * self.fallen[i]
            - p.melting_energy(self.interface[i]) * h
            - p.snow_melting_energy(mean) * self.snow[i])
            - self.limited[i]
    }

    /// What the surface exchanges with the column below over `dt` seconds:
    /// a conductance (W m⁻² K⁻¹) and the temperature it pulls towards (K),
    /// with the interface solved along with it.
    pub(crate) fn balance(&self, i: usize, dt: f64) -> (f64, f64) {
        let p = &self.params;
        let column = self.column(self.thickness[i], self.snow[i]);
        if let Some(snow) = column.snow {
            let (a, b) = self.interface_response(i, &column, snow, dt);
            let storage = column.surface / dt;
            let conductance = storage + snow * (1.0 - b);
            let pull = (storage * self.temperature[i] + snow * a) / conductance;
            (conductance, pull)
        } else {
            let storage = column.interface / dt;
            let conductance = storage + column.ice;
            let pull = (storage * self.temperature[i] + column.ice * p.freezing) / conductance;
            (conductance, pull)
        }
    }

    /// The interface temperature after `dt` seconds is `a + b × surface`.
    fn interface_response(&self, i: usize, column: &Column, snow: f64, dt: f64) -> (f64, f64) {
        let storage = column.interface / dt;
        let denominator = storage + column.ice + snow;
        (
            (storage * self.interface[i] + column.ice * self.params.freezing) / denominator,
            snow / denominator,
        )
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
                snowfall,
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
            // What blows off melts in the water.
            let blown = p.snow_blown * (1.0 - f);
            water.mixed[i] -= dt * f * blown * snowfall * p.fusion / water.capacity;
            self.fallen[i] += dt * (1.0 - blown) * snowfall;

            let column = self.column(h, self.snow[i]);
            let conducted = if let Some(snow) = column.snow {
                let (a, b) = self.interface_response(i, &column, snow, dt);
                let interface = a + b * skin;
                let mut surface =
                    self.temperature[i] + dt * (net + snow * (interface - skin)) / column.surface;
                let mut interface_now = interface;
                if surface > p.melting {
                    self.top[i] += column.surface * (surface - p.melting);
                    surface = p.melting;
                }
                if interface_now > p.melting {
                    self.top[i] += column.interface * (interface_now - p.melting);
                    interface_now = p.melting;
                }
                self.temperature[i] = surface;
                self.interface[i] = interface_now;
                column.ice * (p.freezing - interface)
            } else {
                let conducted = column.ice * (p.freezing - skin);
                let mut t = self.interface[i] + dt * (net + conducted) / column.interface;
                if t > p.melting {
                    self.top[i] += column.interface * (t - p.melting);
                    t = p.melting;
                }
                self.temperature[i] = t;
                self.interface[i] = t;
                conducted
            };
            self.top[i] += dt * sublimation;
            self.bottom[i] += dt * (turbulent - conducted);
        }
    }

    /// Spend the heat recorded since the last call on growth and melting,
    /// lay the snow that has landed, freeze open water that has fallen below
    /// freezing, and remove what is too little to keep.
    pub(crate) fn rearrange(&mut self, mixed: &mut [f64], area: &[f64], capacity: f64) {
        let p = self.params;
        for i in 0..self.fraction.len() {
            if area[i] <= 0.0 {
                continue;
            }
            let (start, mut f, mut h) = (self.fraction[i], self.fraction[i], self.thickness[i]);
            let (mut surface, mut interface) = (self.temperature[i], self.interface[i]);
            let (bottom, top) = (self.bottom[i], self.top[i]);
            let fallen = self.fallen[i];
            self.bottom[i] = 0.0;
            self.top[i] = 0.0;
            self.fallen[i] = 0.0;
            let mut s = self.snow[i];

            // New snow counts as water at the freezing point less its
            // fusion; cooling it to the snow's temperature warms the snow.
            if f > 0.0 && fallen > 0.0 {
                let warmth = fallen * p.specific_heat * (p.freezing - 0.5 * (surface + interface));
                s += fallen;
                let column = self.column(h, s);
                if column.snow.is_some() {
                    surface += warmth / column.surface;
                } else {
                    interface += warmth / column.interface;
                    surface = interface;
                }
            }
            let snow_energy = p.snow_melting_energy(0.5 * (surface + interface));
            let e = p.melting_energy(interface);

            // Masses and volumes are per unit area of ocean from here on.
            let mut mass = start * s;
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
                // Melting takes the snow first, then the ice.
                let mut energy = top * start;
                let melted = mass.min(energy / snow_energy).max(0.0);
                mass -= melted;
                energy -= melted * snow_energy;
                if f > 0.0 {
                    let mut v = energy / e;
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
                    to_ocean += energy;
                }
            }
            mixed[i] += to_ocean / capacity;

            // Snow weighing the ice below the waterline floods and freezes;
            // the sea water's freezing warms the top of the ice.
            if f > 0.0 && h > p.min_thickness {
                let load = mass / f + p.density * h - p.sea_density * h;
                if load > 0.0 {
                    let dh = load / (p.sea_density - p.density + p.snow_density);
                    mass -= f * p.snow_density * dh;
                    h += dh;
                    let warmth = dh * (e - p.snow_density * snow_energy);
                    interface += warmth / (p.capacity(h) + 0.5 * p.specific_heat * mass / f);
                }
            }

            if mixed[i] < p.freezing {
                if f <= 0.0 {
                    surface = p.freezing;
                    interface = p.freezing;
                }
                let formed = p.melting_energy(interface);
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
                let mean = 0.5 * (surface + interface);
                mixed[i] -= (p.melting_energy(interface) * h * f
                    + p.snow_melting_energy(mean) * mass)
                    / capacity;
                f = 0.0;
                mass = 0.0;
            }
            if f <= 0.0 {
                // Snow with no ice left to lie on melts in the water.
                let mean = 0.5 * (surface + interface);
                mixed[i] -= p.snow_melting_energy(mean) * mass / capacity;
                f = 0.0;
                h = 0.0;
                s = 0.0;
                surface = p.freezing;
                interface = p.freezing;
            } else {
                s = mass / f;
                // Too little snow for a layer of its own: one temperature,
                // with the heat the two held.
                if s <= p.snow_critical && surface != interface {
                    let ice = p.capacity(h);
                    let snow = p.specific_heat * s;
                    let mean =
                        (ice * interface + 0.5 * snow * (surface + interface)) / (ice + snow);
                    surface = mean;
                    interface = mean;
                }
            }
            self.fraction[i] = f;
            self.thickness[i] = h;
            self.snow[i] = s;
            self.temperature[i] = surface;
            self.interface[i] = interface;
        }
    }
}
