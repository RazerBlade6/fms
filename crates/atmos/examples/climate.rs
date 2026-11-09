//! The layered model (the climate's four layers unless `FMS_LEVELS` gives
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
//! `FMS_TOP_LAPSE` scales the lapse rate that extrapolates the layer above
//! to the boundary layer's top.
//! `FMS_DRAG` scales the surface drag coefficient and `FMS_DIFFUSION_TIME`
//! the hyperdiffusion's e-folding time, for sensitivity experiments. The run
//! ends with Hadley-cell diagnostics per latitude band.
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
use fms_grid::{LatitudeBands, LloydOptions, MeshGeometry, Operators, SphereMesh};
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
    params.surface_layer.drag_factor = env("FMS_DRAG", 1.0);
    params.surface_layer.top_lapse_factor = env("FMS_TOP_LAPSE", 1.0);
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
    let mut diffusion = HyperdiffusionOptions::strongest_for_step(model.operators(), 4, dt);
    diffusion.e_folding *= env("FMS_DIFFUSION_TIME", 1.0);
    let start_surface = if let Some(sst) = fixed_sst {
        mesh.cell_centers()
            .iter()
            .map(|&p| sst.at(fms_grid::coords::to_lat_lon(p).0))
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
    // Per band and layer: zonal means of u·v, θ, q, Φ, T, v·q, and the
    // product of the zonal means of v and q.
    let mut layer_means = vec![vec![[0.0f64; 7]; n_layers]; nb];
    // Time-mean northward wind and humidity per layer and cell.
    let mut cell_v = vec![vec![0.0f64; nc]; n_layers];
    let mut cell_q = vec![vec![0.0f64; nc]; n_layers];
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
                    .map(|(&p, &v)| fms_grid::zonal::east_north(p, v))
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
                let uv: Vec<f64> = east.iter().zip(&north).map(|(e, n)| e * n).collect();
                for (b, x) in bands.mean(&geometry, &uv).into_iter().enumerate() {
                    layer_means[b][k][0] += x;
                }
                let st = si.state();
                let q: Vec<f64> = (0..nc).map(|i| st.q_ps[k][i] / st.ps[i]).collect();
                let vq: Vec<f64> = north.iter().zip(&q).map(|(n, q)| n * q).collect();
                for (b, x) in bands.mean(&geometry, &vq).into_iter().enumerate() {
                    layer_means[b][k][5] += x;
                }
                let (north_band, q_band) =
                    (bands.mean(&geometry, &north), bands.mean(&geometry, &q));
                for b in 0..nb {
                    layer_means[b][k][6] += north_band[b] * q_band[b];
                }
                for i in 0..nc {
                    cell_v[k][i] += north[i];
                    cell_q[k][i] += q[i];
                }
            }
            {
                let st = si.state();
                let model = si.model();
                let (mut th, mut q, mut phi, mut t) = (
                    vec![vec![0.0; nc]; n_layers],
                    vec![vec![0.0; nc]; n_layers],
                    vec![vec![0.0; nc]; n_layers],
                    vec![vec![0.0; nc]; n_layers],
                );
                for i in 0..nc {
                    let ps = st.ps[i];
                    let theta: Vec<f64> = (0..n_layers).map(|k| st.theta_ps[k][i] / ps).collect();
                    let g = model.geopotentials_at(ps, &theta, 0.0);
                    let ex = model.layer_exner(ps);
                    for k in 0..n_layers {
                        th[k][i] = theta[k];
                        q[k][i] = st.q_ps[k][i] / ps;
                        phi[k][i] = g[k];
                        t[k][i] = theta[k] * ex[k];
                    }
                }
                for k in 0..n_layers {
                    for (j, f) in [&th[k], &q[k], &phi[k], &t[k]].into_iter().enumerate() {
                        for (b, x) in bands.mean(&geometry, f).into_iter().enumerate() {
                            layer_means[b][k][j + 1] += x;
                        }
                    }
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

    // Hadley-cell diagnostics: the balances that set the overturning.
    let field =
        |name: &str, b: usize| zonal[names.iter().position(|n| *n == name).unwrap()][b] / samples;
    let cp = Gas::EARTH_AIR.heat_capacity;
    let two_pi_a = 2.0 * std::f64::consts::PI * RADIUS;
    let omega = planet.rotation_rate();
    let m: Vec<f64> = ds.iter().map(|d| d * ps_mean / GRAVITY).collect();
    let bottom = n_layers - 1;
    let mean_of = |b: usize, k: usize, j: usize| layer_means[b][k][j] / samples;
    let wind = |b: usize, k: usize, c: usize| winds[b][2 * k + c] / samples;
    let mse = |b: usize, k: usize| {
        cp * mean_of(b, k, 4) + mean_of(b, k, 3) + LATENT_HEAT * mean_of(b, k, 2)
    };
    let dse = |b: usize, k: usize| cp * mean_of(b, k, 4) + mean_of(b, k, 3);
    // Column momentum flux (kg s⁻²): eddy and total.
    let flux = |b: usize, eddy: bool| -> f64 {
        (0..n_layers)
            .map(|k| {
                let total = mean_of(b, k, 0);
                m[k] * if eddy {
                    total - wind(b, k, 0) * wind(b, k, 1)
                } else {
                    total
                }
            })
            .sum()
    };
    // Atmospheric column energy input (W m⁻²).
    let input = |b: usize| {
        let toa = field("rad.absorbed_shortwave", b) - field("rad.outgoing_longwave", b);
        let surface = field("rad.surface_shortwave", b)
            - field("rad.surface_longwave", b)
            - field("atm.sensible_heat_flux", b)
            - field("atm.latent_heat_flux", b);
        toa - surface
    };
    let edges: Vec<f64> = (0..=nb)
        .map(|j| -std::f64::consts::FRAC_PI_2 + j as f64 * std::f64::consts::PI / nb as f64)
        .collect();
    let band_area =
        |b: usize| two_pi_a * RADIUS * (fms_math::sin(edges[b + 1]) - fms_math::sin(edges[b]));
    let mean_input: f64 = (0..nb).map(|b| input(b) * band_area(b)).sum::<f64>() / area;
    // Implied northward energy transport (PW) at each band's northern edge.
    let mut implied = vec![0.0; nb];
    let mut acc = 0.0;
    for (b, x) in implied.iter_mut().enumerate() {
        acc += (input(b) - mean_input) * band_area(b);
        *x = acc / 1e15;
    }
    println!(
        "\n  Hadley diagnostics (psi in 1e9 kg/s, northward-positive flow above the level; stress mN/m²; [u'v'] m²/s²; heating K/day of the free troposphere; descent g/m²/s; energy W/m², PW, kJ/kg)"
    );
    println!(
        "  lat   psi_0.5 psi_bl psi_ek  tau_x tau_tot tau_edd  uv'_0  uv'_1  Q_clr  Q_all descent  F_atm  T_impl  T_mmc   dh dh_eff"
    );
    let dp_free = levels_bottom_top(&ds) * ps_mean;
    // The layers above σ = 0.5, for a streamfunction that means the same
    // thing however the layers are divided.
    let upper = {
        let mut sigma = 0.0;
        ds.iter()
            .take_while(|d| {
                sigma += **d;
                sigma <= 0.5 + 1e-9
            })
            .count()
            .max(1)
    };
    for b in 0..nb {
        let phi = lat[b];
        let cosp = fms_math::cos(phi);
        let f = 2.0 * omega * fms_math::sin(phi);
        let psi_top: f64 =
            two_pi_a * cosp * (0..upper).map(|k| m[k] * wind(b, k, 1)).sum::<f64>() / 1e9;
        let psi_bl: f64 =
            two_pi_a * cosp * (0..bottom).map(|k| m[k] * wind(b, k, 1)).sum::<f64>() / 1e9;
        let tau = field("atm.surface_stress_east", b);
        let psi_ek = if phi.abs() > 8f64.to_radians() {
            two_pi_a * cosp * tau / f / 1e9
        } else {
            f64::NAN
        };
        // τ = (1/(a cos²φ)) ∂(cos²φ F)/∂φ, centred over neighbouring bands.
        let implied_tau = |eddy: bool| {
            if b == 0 || b + 1 == nb {
                return f64::NAN;
            }
            let g = |c: usize| fms_math::cos(lat[c]).powi(2) * flux(c, eddy);
            (g(b + 1) - g(b - 1)) / (lat[b + 1] - lat[b - 1]) / (RADIUS * cosp * cosp)
        };
        let to_k_day = |w: f64| w * GRAVITY / (cp * dp_free) * DAY;
        let q_clr = field("rad.clear_free_troposphere_heating", b);
        let q_all = field("rad.free_troposphere_heating", b);
        let descent = -q_all / (dse(b, 0) - dse(b, bottom - 1)) * 1e3;
        let t_mmc: f64 = two_pi_a
            * cosp
            * (0..n_layers)
                .map(|k| m[k] * wind(b, k, 1) * mse(b, k))
                .sum::<f64>()
            / 1e15;
        let t_impl = if b == 0 {
            0.5 * implied[0]
        } else {
            0.5 * (implied[b - 1] + implied[b])
        };
        // Gross moist stability, whatever the layering (kJ/kg): the h of the
        // flow above the boundary layer, weighted by its mass flux, minus the
        // boundary layer's.
        let dh_eff = (0..bottom)
            .map(|k| m[k] * wind(b, k, 1) * (mse(b, k) - mse(b, bottom)))
            .sum::<f64>()
            / (0..bottom).map(|k| m[k] * wind(b, k, 1)).sum::<f64>()
            / 1e3;
        println!(
            "  {:>5.1} {:>7.0} {:>6.0} {:>6.0} {:>6.1} {:>7.1} {:>7.1} {:>6.1} {:>6.1} {:>6.2} {:>6.2} {:>7.2} {:>6.1} {:>7.2} {:>6.2} {:>5.1} {:>6.1}",
            phi.to_degrees(),
            psi_top,
            psi_bl,
            psi_ek,
            1e3 * tau,
            1e3 * implied_tau(false),
            1e3 * implied_tau(true),
            mean_of(b, 0, 0) - wind(b, 0, 0) * wind(b, 0, 1),
            mean_of(b, 1, 0) - wind(b, 1, 0) * wind(b, 1, 1),
            to_k_day(q_clr),
            to_k_day(q_all),
            descent,
            input(b),
            t_impl,
            t_mmc,
            (mse(b, 0) - mse(b, bottom)) / 1e3,
            dh_eff,
        );
    }
    let rows: Vec<(f64, Vec<[f64; 4]>)> = (nb / 2..nb)
        .map(|b| {
            let layers = (0..n_layers)
                .map(|k| {
                    [
                        mse(b, k) / 1e3,
                        cp * mean_of(b, k, 4) / 1e3,
                        mean_of(b, k, 3) / 1e3,
                        LATENT_HEAT * mean_of(b, k, 2) / 1e3,
                    ]
                })
                .collect();
            (lat[b].to_degrees(), layers)
        })
        .collect();
    print_mse(&rows);

    println!(
        "\n  Water budget of the rain (mm/day). Boundary layer: E − export + conv − cond + reevap + adv = 0 (adv is the residual)"
    );
    println!(
        "  lat        E export   conv   cond reevap    adv | cond_abv cond_up reev_abv  deep      P"
    );
    for (b, latitude) in lat.iter().enumerate() {
        let f = |name: &str| field(name, b);
        let (e, export, conv, cond, re) = (
            f("atm.evaporation"),
            f("atm.boundary_layer_export"),
            f("moist.convective_boundary"),
            f("moist.condensation_boundary"),
            f("moist.reevaporation_boundary"),
        );
        println!(
            "  {:>5.1} {:>6.2} {:>6.2} {:>6.2} {:>6.2} {:>6.2} {:>6.2} | {:>8.2} {:>7.2} {:>8.2} {:>5.2} {:>6.2}",
            latitude.to_degrees(),
            e,
            export,
            conv,
            cond,
            re,
            -(e - export + conv - cond + re),
            f("moist.condensation_above"),
            f("moist.condensation_upper"),
            f("moist.reevaporation_above"),
            f("atm.convective_precipitation"),
            f("atm.precipitation"),
        );
    }

    // Northward water transport (10⁸ kg/s): the mean flow in the boundary
    // layer and above it, the eddies, and what E − P implies.
    let mut implied_water = vec![0.0; nb];
    let net: f64 = (0..nb)
        .map(|b| (field("atm.evaporation", b) - field("atm.precipitation", b)) * band_area(b))
        .sum::<f64>()
        / area;
    let mut acc = 0.0;
    for (b, x) in implied_water.iter_mut().enumerate() {
        acc += (field("atm.evaporation", b) - field("atm.precipitation", b) - net) * band_area(b)
            / DAY;
        *x = acc / 1e8;
    }
    // [v̄* q̄*] per band and layer: the stationary eddies' flux.
    let stationary: Vec<Vec<f64>> = (0..n_layers)
        .map(|k| {
            let product: Vec<f64> = (0..nc)
                .map(|i| cell_v[k][i] / samples * cell_q[k][i] / samples)
                .collect();
            let means = bands.mean(&geometry, &product);
            (0..nb)
                .map(|b| means[b] - wind(b, k, 1) * mean_of(b, k, 2))
                .collect()
        })
        .collect();
    println!(
        "\n  Northward water transport (1e8 kg/s) through each band's centre: time-mean zonal flow (boundary layer, above), its variations in time, stationary and transient eddies (and their boundary-layer part), total, E − P implied"
    );
    println!("  lat   mean_bl mean_abv  var_zm station transnt tran_bl   total implied");
    for b in 0..nb {
        let factor = two_pi_a * fms_math::cos(lat[b]) / 1e8;
        let mean_flow = |k: usize| m[k] * wind(b, k, 1) * mean_of(b, k, 2);
        let sum =
            |f: &dyn Fn(usize) -> f64| factor * (0..n_layers).map(|k| m[k] * f(k)).sum::<f64>();
        let mean_bl = factor * mean_flow(bottom);
        let mean_above = factor * (0..bottom).map(mean_flow).sum::<f64>();
        let total = sum(&|k| mean_of(b, k, 5));
        let variation = sum(&|k| mean_of(b, k, 6) - wind(b, k, 1) * mean_of(b, k, 2));
        let station = sum(&|k| stationary[k][b]);
        let transient = sum(&|k| mean_of(b, k, 5) - mean_of(b, k, 6) - stationary[k][b]);
        let transient_bl = factor
            * m[bottom]
            * (mean_of(b, bottom, 5) - mean_of(b, bottom, 6) - stationary[bottom][b]);
        let t_impl = if b == 0 {
            0.5 * implied_water[0]
        } else {
            0.5 * (implied_water[b - 1] + implied_water[b])
        };
        println!(
            "  {:>5.1} {:>7.2} {:>8.2} {:>7.2} {:>7.2} {:>7.2} {:>7.2} {:>7.2} {:>7.2}",
            lat[b].to_degrees(),
            mean_bl,
            mean_above,
            variation,
            station,
            transient,
            transient_bl,
            total,
            t_impl,
        );
    }
}

/// Moist static energy per layer (kJ/kg) and its parts, at each band.
fn print_mse(rows: &[(f64, Vec<[f64; 4]>)]) {
    println!("\n  Moist static energy by layer (kJ/kg): h = c_p T + Φ + L q");
    for (lat, layers) in rows {
        let parts: Vec<String> = layers
            .iter()
            .enumerate()
            .map(|(k, [h, t, phi, lq])| format!("h_{k} {h:.1} (c_pT {t:.1} Φ {phi:.1} Lq {lq:.1})"))
            .collect();
        println!("  {lat:>5.1}  {}", parts.join("  "));
    }
}

/// σ thickness of the layers above the boundary layer.
fn levels_bottom_top(ds: &[f64]) -> f64 {
    ds[..ds.len() - 1].iter().sum()
}
