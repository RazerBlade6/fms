//! A baroclinic storm life cycle with and without hyperdiffusion: survival,
//! peak winds, the energy budget, and the damping of large-scale waves.
//!
//! ```sh
//! cargo run --release -p fms-atmos --example storm [frequency] [days]
//! ```
//!
//! Environment: `FMS_ORDER` (default 4), `FMS_TAU` (grid-scale e-folding
//! time in hours; default: the strongest the step allows with a 2× margin), `FMS_SHEAR` (default 40 m/s), `FMS_MAX_WIND` (design wind for the
//! step, default 80 m/s), `FMS_INVISCID`
//! (also run without diffusion).

use fms_atmos::williamson::{GRAVITY, OMEGA, RADIUS};
use fms_atmos::{
    Gas, HyperdiffusionOptions, ShearedFlow, TwoLayer, TwoLayerSemiImplicit,
    TwoLayerSemiImplicitOptions,
};
use fms_geo::coords::to_lat_lon;
use fms_geo::{LloydOptions, MeshGeometry, Operators, SphereMesh};

const DAY: f64 = 86_400.0;

fn env<T: std::str::FromStr>(name: &str, default: T) -> T {
    std::env::var(name)
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(default)
}

fn main() {
    let mut args = std::env::args().skip(1);
    let k: u32 = args.next().map_or(32, |s| s.parse().unwrap());
    let days: f64 = args.next().map_or(30.0, |s| s.parse().unwrap());
    let order: u32 = env("FMS_ORDER", 4);
    let tau: Option<f64> = std::env::var("FMS_TAU").ok().map(|s| s.parse().unwrap());
    let shear: f64 = env("FMS_SHEAR", 40.0);
    let (mesh, _) = SphereMesh::centroidal(k, LloydOptions::default());
    let geometry = MeshGeometry::new(&mesh);
    println!("frequency {k}: {} cells", mesh.n_cells());

    let inviscid = std::env::var("FMS_INVISCID").is_ok();
    for diffuse in [false, true] {
        if !diffuse && !inviscid {
            continue;
        }
        let mut model = TwoLayer::new(
            Operators::new(&mesh, &geometry, RADIUS),
            GRAVITY,
            Gas::EARTH_AIR,
            0.5,
            mesh.vertices().iter().map(|p| 2.0 * OMEGA * p.z).collect(),
            vec![0.0; mesh.n_cells()],
        );
        let dt = TwoLayerSemiImplicit::stable_dt(&model, env("FMS_MAX_WIND", 80.0), 0.6);
        let dt = DAY / (DAY / dt).ceil();
        if diffuse {
            let options = match tau {
                Some(h) => HyperdiffusionOptions {
                    order,
                    e_folding: h * 3600.0,
                },
                None => HyperdiffusionOptions::strongest_for_step(model.operators(), order, dt),
            };
            model = model.with_hyperdiffusion(options);
            let d = model.hyperdiffusion().unwrap();
            let [lc, lv] = d.largest_eigenvalues();
            let spacing = mesh.mean_spacing() * RADIUS;
            println!(
                "∇{order}, τ = {:.1} h: λ_max·d² = {:.2} (cells), {:.2} (vertices); ν = {:.2e}; fastest e-folding {:.0} s",
                options.e_folding / 3600.0,
                lc * spacing * spacing,
                lv * spacing * spacing,
                d.coefficient(),
                1.0 / d.max_rate(),
            );
            for l in [6, 12, 24, 48] {
                println!(
                    "  degree {l:>2}: e-folding {:.1} days",
                    1.0 / (d.rate_at_degree(l) * DAY)
                );
            }
        } else {
            println!("inviscid");
        }
        let flow = ShearedFlow {
            speed: [10.0 + shear, 10.0],
            stability: 40.0,
            theta_equator: 315.0,
            ps_equator: 1.0e5,
            omega: OMEGA,
        };
        let mut state = flow.state(&model);
        for k in 0..2 {
            for (i, &p) in mesh.cell_centers().iter().enumerate() {
                let (lat, lon) = to_lat_lon(p);
                let band =
                    |c: f64| fms_math::exp(-((lat - c.to_radians()) / 15f64.to_radians()).powi(2));
                let waves: f64 = (4..=8)
                    .map(|m| fms_math::cos(f64::from(m) * lon + 0.7 * f64::from(m)))
                    .sum();
                state.theta_ps[k][i] += state.ps[i] * 0.5 * (band(45.0) + band(-45.0)) * waves;
            }
        }
        let start = model.diagnostics(&state);
        let mut si = TwoLayerSemiImplicit::new(
            model,
            state,
            dt,
            TwoLayerSemiImplicitOptions::default(),
            0.5,
        );
        println!("dt {dt:.0} s, Courant limit {:.0} m/s", si.courant_limit());
        let per_day = (DAY / dt).round() as u64;
        let total = (days * per_day as f64).round() as u64;
        let mut peak: f64 = 0.0;
        while si.steps() < total {
            si.step();
            match si.check() {
                Ok(h) => peak = peak.max(h.max_wind),
                Err(b) => {
                    println!("  BLOW-UP day {:.2}: {b}", si.steps() as f64 * dt / DAY);
                    break;
                }
            }
            if si.steps() % per_day == 0 {
                let d = si.model().diagnostics(si.state());
                let h = si.check().unwrap();
                let de = d.energy - start.energy;
                let (e, layer) = h.max_wind_at;
                let (lat, lon) = to_lat_lon(geometry.edge_points[e]);
                // Neighbouring-edge contrast: how grid-scale is the peak?
                let u = &si.state().u[layer];
                let mut zeta = vec![0.0; mesh.n_vertices()];
                si.model().operators().curl(u, &mut zeta);
                let zmax = zeta.iter().fold(0.0f64, |a, z| a.max(z.abs()));
                println!(
                    "      peak at {:.0}°, {:.0}° (layer {layer}); max |ζ| {:.1e} s⁻¹",
                    lat.to_degrees(),
                    lon.to_degrees(),
                    zmax
                );
                println!(
                    "  day {:>4.1}: max |u| {:>5.1}  KE {:.3e} J  ΔE {:+.2e}  dissipated {:+.2e}  residual {:+.2e} ({:+.1e} of E)  p_s {:.0}–{:.0} hPa",
                    si.steps() as f64 * dt / DAY,
                    h.max_wind,
                    d.kinetic_energy,
                    de,
                    si.dissipated_energy(),
                    de - si.dissipated_energy(),
                    (de - si.dissipated_energy()) / start.energy,
                    h.surface_pressure.0 / 100.0,
                    h.surface_pressure.1 / 100.0,
                );
            }
        }
        let end = si.model().diagnostics(si.state());
        println!(
            "  peak |u| {peak:.1} m/s; mass {:+.1e}, Θ {:+.1e}",
            end.mass / start.mass - 1.0,
            end.theta / start.theta - 1.0
        );
    }
}
