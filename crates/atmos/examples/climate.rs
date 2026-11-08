//! The layered model (the climate's three layers unless `FMS_LEVELS` gives
//! other interior interfaces, comma-separated) driven by radiation and moist
//! physics over a slab aquaplanet (Earth's orbit and CO₂): global means,
//! cloud effects, the energy and water budgets, and zonal means. For
//! calibration.
//!
//! ```sh
//! cargo run --release -p fms-atmos --example climate [frequency] [years]
//! ```
//!
//! `FMS_SST=control` (or `qobs`) holds the sea-surface temperature at that
//! Aqua-Planet Experiment profile on a perpetual equinox (it settles in about
//! a month; run e.g. `0.3` years with `FMS_AVERAGE=60`; compare with
//! Blackburn et al. 2013); `FMS_AVERAGE` sets the
//! final averaging window in days (default a year).
//!
//! Environment overrides (see `PhysicsParams` and `MoistureParams`):
//! `FMS_CO2`, `FMS_MAX_WIND`, `FMS_ALBEDO_HIGH`, `FMS_ALBEDO_LOW`,
//! `FMS_EMISSIVITY_HIGH`, `FMS_EMISSIVITY_LOW`, `FMS_RHC_HIGH`, `FMS_RHC_LOW`,
//! `FMS_ANVIL`, `FMS_CONDENSE_HIGH`, `FMS_CONDENSE_LOW`, `FMS_CONDENSE_BL`, `FMS_CONV_RH`,
//! `FMS_RHO_HIGH`, `FMS_RHO_LOW`, `FMS_REEVAP`, `FMS_CONV_TIME`,
//! `FMS_TRANSFER`, `FMS_SLAB`, `FMS_RAYLEIGH`, `FMS_REPORTS` (progress lines
//! per year).

use std::time::Instant;

use fms_atmos::column_physics::{ApeSst, CLIMATE_LEVELS, LATENT_HEAT, initial_state};
use fms_atmos::williamson::{GRAVITY, RADIUS};
use fms_atmos::{
    ColumnPhysics, Gas, HyperdiffusionOptions, LayeredModel, LayeredSemiImplicit,
    LayeredSemiImplicitOptions, Levels, PhysicsParams,
};
use fms_geo::{LatitudeBands, LloydOptions, MeshGeometry, Operators, SphereMesh};
use fms_planet::{Planet, PlanetParams};

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
    // Fixed sea-surface temperature (the Aqua-Planet Experiment's QOBS
    // profile on a perpetual equinox), or Earth's orbit over the slab.
    let fixed_sst = match std::env::var("FMS_SST").as_deref() {
        Ok("control") => Some(ApeSst::Control),
        Ok("qobs") => Some(ApeSst::Qobs),
        Ok(other) => panic!("FMS_SST: control or qobs, not {other}"),
        Err(_) => None,
    };
    let planet = if fixed_sst.is_some() {
        Planet::new(PlanetParams {
            obliquity_deg: 0.0,
            eccentricity: 0.0,
            ..PlanetParams::earth()
        })
        .0
    } else {
        Planet::earth()
    };
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
        m.condensation_humidity[2] = env("FMS_CONDENSE_BL", m.condensation_humidity[2]);
        m.convection_humidity = env("FMS_CONV_RH", m.convection_humidity);
        m.reevaporation = env("FMS_REEVAP", m.reevaporation);
        m.convection_time = env("FMS_CONV_TIME", m.convection_time);
    }
    println!("{:?}", params.moisture);

    let levels = std::env::var("FMS_LEVELS").map_or(Levels::new(&CLIMATE_LEVELS), |v| {
        let interior: Vec<f64> = v.split(',').map(|s| s.parse().unwrap()).collect();
        Levels::new(&interior)
    });
    let n_layers = levels.count();
    println!("levels {:?}", levels.interfaces());
    // Start from the single-column equilibrium (fixed relative humidity),
    // at rest.
    let model = LayeredModel::new(
        Operators::new(&mesh, &geometry, RADIUS),
        GRAVITY,
        Gas::EARTH_AIR,
        levels,
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
    let dt = LayeredSemiImplicit::stable_dt(&model, env("FMS_MAX_WIND", 120.0), 0.6);
    let dt = DAY / (DAY / dt).ceil();
    let diffusion = HyperdiffusionOptions::strongest_for_step(model.operators(), 4, dt);
    let start_surface = if let Some(sst) = fixed_sst {
        mesh.cell_centers()
            .iter()
            .map(|&p| sst.at(fms_geo::coords::to_lat_lon(p).0))
            .collect()
    } else {
        vec![surface_temperature; nc]
    };
    let mut physics = ColumnPhysics::new(&model, &planet, params, 0.0, &state, start_surface);
    if let Some(sst) = fixed_sst {
        physics = physics.with_prescribed_surface();
        println!("fixed SST: APE {sst:?}, perpetual equinox");
    }
    let model = model
        .with_hyperdiffusion(diffusion)
        .with_dissipative_heating()
        .with_forcing(Box::new(physics));
    let mut si =
        LayeredSemiImplicit::new(model, state, dt, LayeredSemiImplicitOptions::default(), 0.5);
    let start = si.model().diagnostics(si.state());
    let water0 = si.water();
    let budget0 = si.model().forcing().unwrap().budget().unwrap();
    let (energy0, water_budget0) = (si.energy_budget(), si.water_budget());
    let year = planet.year_seconds();
    let per_year = (year / dt).round() as u64;
    let total = (years * per_year as f64) as u64;
    // Averaging window at the end of the run (days; default a year).
    let average = env("FMS_AVERAGE", year / DAY);
    let window = ((average * DAY / dt).round() as u64).min(total - 1);
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
    let mut winds = vec![vec![0.0f64; 2 * n_layers]; nb];
    let t = Instant::now();
    while si.steps() < total {
        si.step();
        if let Err(b) = si.check() {
            println!("BLOW-UP day {:.1}: {b}", si.steps() as f64 * dt / DAY);
            return;
        }
        let sampling = si.steps() + window >= total && si.steps() % 8 == 0;
        let reporting = si.steps() % (per_year / env("FMS_REPORTS", 8.0) as u64) == 0;
        if si.steps() + window == total {
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
            for k in 0..n_layers {
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
                - (b.dissipated - budget0.dissipated)
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
    let span = window as f64 * dt;
    let wm2 = |j: f64| j / span / area;
    let asr = wm2(b.absorbed - b0.absorbed);
    let olr = wm2(b.outgoing - b0.outgoing);
    println!(
        "final {average:.0} days: ASR {asr:.1}, OLR {olr:.1}, imbalance {:+.2} W/m², surface uptake {:+.2} W/m², planetary albedo {:.3}; E {:.3} P {:.3} mm/day",
        asr - olr,
        wm2(b.stored - b0.stored),
        1.0 - asr / (planet.mean_flux() / 4.0),
        (b.evaporated - b0.evaporated) / span / area * DAY,
        (b.precipitated - b0.precipitated) / span / area * DAY,
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
        "atm.surface_wind",
        "atm.entrainment",
        "atm.sensible_heat_flux",
        "atm.boundary_layer_moist_heating",
        "atm.shallow_cumulus",
    ];
    let layer_heads: String = (0..n_layers)
        .map(|k| format!("    u_{k}    v_{k}"))
        .collect();
    println!(
        "  lat      T_s      P   conv      E   high    low   RH_b   RH_0    OLR    ASR   V_10    w_e     SH  Qmoist    M_c{layer_heads}    psi"
    );
    // Mass streamfunction at the boundary layer's top (10⁹ kg/s, positive
    // northward flow above it): 2πa cos φ Σ v̄ₖ δσₖ p_s / g over the layers
    // above.
    let ds = si.model().dsigma().to_vec();
    let ps_mean = mean(&si.state().ps);
    let mut psi_max: f64 = 0.0;
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
        let flux: f64 = (0..n_layers - 1)
            .map(|k| winds[b][2 * k + 1] / samples * ds[k] * ps_mean / GRAVITY)
            .sum();
        let psi = 2.0 * std::f64::consts::PI * RADIUS * fms_math::cos(lat[b]) * flux / 1e9;
        psi_max = psi_max.max(psi.abs());
        println!(
            "  {:>5.1}°  {} {} {:>6.0}",
            lat[b].to_degrees(),
            row.join(" "),
            w.join(" "),
            psi
        );
    }
    println!("  largest |psi| at the boundary layer's top: {psi_max:.0}e9 kg/s");
}
