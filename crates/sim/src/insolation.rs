//! Systems that compute top-of-atmosphere sunlight on the climate grid.

use fms_fields::FieldId;
use fms_geo::coords::to_lat_lon;
use fms_planet::Planet;
use fms_sched::{Access, Rate, StepContext, System, SystemDesc, Tiers};
use glam::DVec3;
use rayon::prelude::*;

/// Instantaneous TOA insolation at every climate cell, every step.
pub struct Insolation {
    planet: Planet,
    centers: Vec<DVec3>,
    out: FieldId<f64>,
}

impl Insolation {
    /// A system writing to `out` (one value per climate cell centre).
    pub fn new(planet: Planet, centers: Vec<DVec3>, out: FieldId<f64>) -> Self {
        Self {
            planet,
            centers,
            out,
        }
    }
}

impl System for Insolation {
    fn describe(&self) -> SystemDesc {
        SystemDesc {
            name: "planet.insolation",
            access: Access::new().write(self.out),
            rate: Rate::Every(1),
            tiers: Tiers::ALL,
        }
    }

    fn run(&mut self, ctx: &StepContext<'_>) {
        let t = ctx.date().step as f64 * ctx.dt();
        self.planet
            .insolation_into(&self.centers, t, &mut ctx.write(self.out));
    }
}

/// Daily-mean TOA insolation at every climate cell, once per day.
pub struct DailyInsolation {
    planet: Planet,
    latitudes: Vec<f64>,
    out: FieldId<f64>,
}

impl DailyInsolation {
    /// A system writing to `out` (one value per climate cell centre).
    pub fn new(planet: Planet, centers: Vec<DVec3>, out: FieldId<f64>) -> Self {
        let latitudes = centers.iter().map(|&p| to_lat_lon(p).0).collect();
        Self {
            planet,
            latitudes,
            out,
        }
    }
}

impl System for DailyInsolation {
    fn describe(&self) -> SystemDesc {
        SystemDesc {
            name: "planet.daily_insolation",
            access: Access::new().write(self.out),
            rate: Rate::Daily,
            tiers: Tiers::ALL,
        }
    }

    fn run(&mut self, ctx: &StepContext<'_>) {
        // Evaluate at mid-day, the best single representative of the day.
        let t =
            (ctx.date().step as f64 + 0.5 * f64::from(ctx.calendar().steps_per_day())) * ctx.dt();
        let planet = &self.planet;
        ctx.write(self.out)
            .par_iter_mut()
            .zip(&self.latitudes)
            .for_each(|(o, &lat)| *o = planet.daily_mean_insolation(lat, t));
    }
}
