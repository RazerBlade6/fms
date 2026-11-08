//! The two-layer model driven by radiation and moist physics over a slab
//! aquaplanet (Earth's orbit and CO₂): global means, cloud effects, the
//! energy and water budgets, and zonal means. For calibration.
//!
//! ```sh
//! cargo run --release -p fms-atmos --example climate [frequency] [years]
//! ```
//!
//! Environment overrides (see `PhysicsParams` and `MoistureParams`):
//! `FMS_CO2`, `FMS_MAX_WIND`, `FMS_ALBEDO_HIGH`, `FMS_ALBEDO_LOW`,
//! `FMS_EMISSIVITY_HIGH`, `FMS_EMISSIVITY_LOW`, `FMS_RHC_HIGH`, `FMS_RHC_LOW`,
//! `FMS_ANVIL`, `FMS_CONDENSE_HIGH`, `FMS_CONDENSE_LOW`, `FMS_CONV_RH`,
//! `FMS_RHO_HIGH`, `FMS_RHO_LOW`, `FMS_REEVAP`, `FMS_CONV_TIME`,
//! `FMS_TRANSFER`, `FMS_SLAB`, `FMS_RAYLEIGH`, `FMS_REPORTS` (progress lines
//! per year).

use std::time::Instant;

use fms_atmos::column_physics::{LATENT_HEAT, initial_state};
use fms_atmos::williamson::{GRAVITY, RADIUS};
use fms_atmos::{
    ColumnPhysics, Gas, HyperdiffusionOptions, PhysicsParams, TwoLayer, TwoLayerSemiImplicit,
    TwoLayerSemiImplicitOptions,
};
use fms_geo::{LatitudeBands, LloydOptions, MeshGeometry, Operators, SphereMesh};
use fms_planet::Planet;

const DAY: f64 = 86_400.0;

fn env(name: &str, default: f64) -> f64 {
    std::env::var(name)
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(default)
}

fn main() {
    let mut args = std::env::args().skip(1);
    let k: u32 = args.next().map_or(16, |s| s.parse().unwrap());
    let years: f64 = args.next().map_or(2.0, |s| s.parse().unwrap());
    let planet = Planet::earth();
    let (mesh, _) = SphereMesh::centroidal(k, LloydOptions::default());
    let geometry = MeshGeometry::new(&mesh);
    let mut params = PhysicsParams::default();
    params.co2 = env("FMS_CO2", params.co2);
    params.transfer = env("FMS_TRANSFER", params.transfer);
    params.slab_depth = env("FMS_SLAB", params.slab_depth);
    params.radiation.rayleigh = env("FMS_RAYLEIGH", params.radiation.rayleigh);
    {
        let m = &mut params.moisture;
        let c = &mut m.clouds;
        c.albedo[0] = env("FMS_ALBEDO_HIGH", c.albedo[0]);
        c.albedo[1] = env("FMS_ALBEDO_LOW", c.albedo[1]);
        c.emissivity[0] = env("FMS_EMISSIVITY_HIGH", c.emissivity[0]);
        c.emissivity[1] = env("FMS_EMISSIVITY_LOW", c.emissivity[1]);
        c.critical_humidity[0] = env("FMS_RHC_HIGH", c.critical_humidity[0]);
        c.critical_humidity[1] = env("FMS_RHC_LOW", c.critical_humidity[1]);
        c.overcast_humidity[0] = env("FMS_RHO_HIGH", c.overcast_humidity[0]);
        c.overcast_humidity[1] = env("FMS_RHO_LOW", c.overcast_humidity[1]);
        c.convective = env("FMS_ANVIL", c.convective);
        m.condensation_humidity[0] = env("FMS_CONDENSE_HIGH", m.condensation_humidity[0]);
        m.condensation_humidity[1] = env("FMS_CONDENSE_LOW", m.condensation_humidity[1]);
        m.convection_humidity = env("FMS_CONV_RH", m.convection_humidity);
        m.reevaporation = env("FMS_REEVAP", m.reevaporation);
        m.convection_time = env("FMS_CONV_TIME", m.convection_time);
    }
    println!("{:?}", params.moisture);

    // Start from the single-column equilibrium (fixed relative humidity),
    // at rest.
    let model = TwoLayer::new(
        Operators::new(&mesh, &geometry, RADIUS),
        GRAVITY,
        Gas::EARTH_AIR,
        0.5,
        mesh.vertices()
            .iter()
            .map(|p| 2.0 * planet.rotation_rate() * p.z)
            .collect(),
        vec![0.0; mesh.n_cells()],
    )
    .with_humidity_profile(params.radiation);
    let nc = mesh.n_cells();
    let (state, surface_temperature) = initial_state(&model, &planet, &params);
    println!("single column: T_s {surface_temperature:.1} K");
    let dt = TwoLayerSemiImplicit::stable_dt(&model, env("FMS_MAX_WIND", 120.0), 0.6);
    let dt = DAY / (DAY / dt).ceil();
    let diffusion = HyperdiffusionOptions::strongest_for_step(model.operators(), 4, dt);
    let physics = ColumnPhysics::new(
        &model,
        &planet,
        params,
        0.0,
        &state,
        vec![surface_temperature; nc],
    );
    let model = model
        .with_hyperdiffusion(diffusion)
        .with_forcing(Box::new(physics));
    let mut si = TwoLayerSemiImplicit::new(
        model,
        state,
        dt,
        TwoLayerSemiImplicitOptions::default(),
        0.5,
    );
    let start = si.model().diagnostics(si.state());
    let water0 = si.water();
    let budget0 = si.model().forcing().unwrap().budget().unwrap();
    let (energy0, water_budget0) = (si.energy_budget(), si.water_budget());
    let year = planet.year_seconds();
    let per_year = (year / dt).round() as u64;
    let total = (years * per_year as f64) as u64;
    let nb = 18;
    let bands = LatitudeBands::new(&mesh, &geometry, nb);
    let area = 4.0 * std::f64::consts::PI * RADIUS * RADIUS;
    let mean = |v: &[f64]| -> f64 {
        v.iter()
            .zip(&geometry.cell_areas)
            .map(|(t, a)| t * a)
            .sum::<f64>()
            / (4.0 * std::f64::consts::PI)
    };
    // Final-year sums of global means (by field) and zonal means.
    let names: Vec<&str> = si
        .model()
        .forcing()
        .unwrap()
        .cell_fields()
        .iter()
        .map(|f| f.name)
        .collect();
    let mut sums = vec![0.0; names.len()];
    let mut zonal = vec![vec![0.0; nb]; names.len()];
    let mut samples = 0.0;
    let mut last_year = None;
    // Zonal-mean winds (east, north) per layer.
    let mut winds = vec![[0.0f64; 4]; nb];
    let t = Instant::now();
    while si.steps() < total {
        si.step();
        if let Err(b) = si.check() {
            println!("BLOW-UP day {:.1}: {b}", si.steps() as f64 * dt / DAY);
            return;
        }
        let sampling = si.steps() + per_year >= total && si.steps() % 8 == 0;
        let reporting = si.steps() % (per_year / env("FMS_REPORTS", 8.0) as u64) == 0;
        if si.steps() + per_year == total {
            last_year = Some(si.model().forcing().unwrap().budget().unwrap());
        }
        if !(sampling || reporting) {
            continue;
        }
        let fields = si.model().forcing().unwrap().cell_fields();
        if sampling {
            for (j, f) in fields.iter().enumerate() {
                sums[j] += mean(&f.values);
                for (z, m) in zonal[j].iter_mut().zip(bands.mean(&geometry, &f.values)) {
                    *z += m;
                }
            }
            samples += 1.0;
            for k in 0..2 {
                let v = geometry.reconstruct_vectors(&mesh, &si.state().u[k]);
                let (east, north): (Vec<f64>, Vec<f64>) = mesh
                    .cell_centers()
                    .iter()
                    .zip(&v)
                    .map(|(&p, &v)| fms_geo::zonal::east_north(p, v))
                    .unzip();
                for (b, (e, n)) in bands
                    .mean(&geometry, &east)
                    .into_iter()
                    .zip(bands.mean(&geometry, &north))
                    .enumerate()
                {
                    winds[b][2 * k] += e;
                    winds[b][2 * k + 1] += n;
                }
            }
        }
        if reporting {
            let b = si.model().forcing().unwrap().budget().unwrap();
            let d = si.model().diagnostics(si.state());
            let e = si.energy_budget();
            let w = si.water_budget();
            let water = si.water();
            let elapsed = si.steps() as f64 * dt;
            let closure = (d.energy - start.energy)
                + (b.stored - budget0.stored)
                + LATENT_HEAT * (water - water0)
                - (b.absorbed - budget0.absorbed - b.outgoing + budget0.outgoing)
                - (e.forcing.kinetic - energy0.forcing.kinetic)
                - (e.diffusion.total() - energy0.diffusion.total());
            let water_closure = (water - water0)
                - (b.evaporated - budget0.evaporated - b.precipitated + budget0.precipitated)
                - (w.filter - water_budget0.filter)
                - (w.fixer - water_budget0.fixer);
            let h = si.check().unwrap();
            let get =
                |name: &str| mean(&fields[names.iter().position(|n| *n == name).unwrap()].values);
            println!(
                "day {:>6.1}: T_s {:.2} K  ASR {:.1} OLR {:.1} CRE {:+.0}/{:+.0}  RH {:.2}/{:.2}  cloud {:.2}/{:.2}  E {:.2} P {:.2} (conv {:.2}) mm/d  water {:.1} kg/m²  max |u| {:.0}  energy residual {:+.4} W/m²  water residual {:+.1e} (filter {:+.1e}, fixer {:+.1e}, substeps {})",
                elapsed / DAY,
                get("rad.surface_temperature"),
                get("rad.absorbed_shortwave"),
                get("rad.outgoing_longwave"),
                get("rad.shortwave_cloud_effect"),
                get("rad.longwave_cloud_effect"),
                get("atm.relative_humidity_upper"),
                get("atm.relative_humidity_lower"),
                get("atm.high_cloud"),
                get("atm.low_cloud"),
                get("atm.evaporation"),
                get("atm.precipitation"),
                get("atm.convective_precipitation"),
                water / area,
                h.max_wind,
                closure / elapsed / area,
                water_closure / water,
                (w.filter - water_budget0.filter) / water,
                (w.fixer - water_budget0.fixer) / water,
                w.max_substeps,
            );
        }
    }
    println!("{years} years in {:.0} s", t.elapsed().as_secs_f64());
    let Some(b0) = last_year else { return };
    let b = si.model().forcing().unwrap().budget().unwrap();
    let wm2 = |j: f64| j / year / area;
    let asr = wm2(b.absorbed - b0.absorbed);
    let olr = wm2(b.outgoing - b0.outgoing);
    println!(
        "last year: ASR {asr:.1}, OLR {olr:.1}, imbalance {:+.2} W/m², slab uptake {:+.2} W/m², planetary albedo {:.3}; E {:.3} P {:.3} mm/day",
        asr - olr,
        wm2(b.stored - b0.stored),
        1.0 - asr / (planet.mean_flux() / 4.0),
        (b.evaporated - b0.evaporated) / year / area * DAY,
        (b.precipitated - b0.precipitated) / year / area * DAY,
    );
    for (j, name) in names.iter().enumerate() {
        println!("  {name:<34} {:>9.3}", sums[j] / samples);
    }
    let lat = bands.latitudes();
    let pick = [
        "rad.surface_temperature",
        "atm.precipitation",
        "atm.convective_precipitation",
        "atm.evaporation",
        "atm.high_cloud",
        "atm.low_cloud",
        "atm.relative_humidity_lower",
        "atm.relative_humidity_upper",
        "rad.outgoing_longwave",
        "rad.absorbed_shortwave",
    ];
    println!(
        "  lat      T_s      P   conv      E   high    low   RH_1   RH_0    OLR    ASR    u_0    v_0    u_1    v_1"
    );
    for b in 0..nb {
        let row: Vec<String> = pick
            .iter()
            .map(|p| {
                let j = names.iter().position(|n| n == p).unwrap();
                format!("{:>6.2}", zonal[j][b] / samples)
            })
            .collect();
        let w: Vec<String> = winds[b]
            .iter()
            .map(|w| format!("{:>6.2}", w / samples))
            .collect();
        println!(
            "  {:>5.1}°  {} {}",
            lat[b].to_degrees(),
            row.join(" "),
            w.join(" ")
        );
    }
}
