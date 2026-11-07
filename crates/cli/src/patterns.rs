//! Analytic test fields for checking the rendering of edge, vertex and vector
//! fields before any real dynamics exist.

use fms_fields::{ColorScale, FieldKey, FieldSpec, Location, Quantity};
use fms_geo::coords::to_lat_lon;
use fms_sim::World;
use glam::DVec3;

/// Rossby–Haurwitz wave parameters (Williamson et al. 1992, test case 6).
const RH_OMEGA: f64 = 7.848e-6;
const RH_K: f64 = 7.848e-6;
const RH_WAVENUMBER: i32 = 4;

/// Stream function of the Rossby–Haurwitz wave (m² s⁻¹) at a unit vector.
fn rossby_haurwitz_stream(p: DVec3, radius: f64) -> f64 {
    let (lat, lon) = to_lat_lon(p);
    let (sin_lat, cos_lat) = fms_math::sin_cos(lat);
    let a2 = radius * radius;
    -a2 * RH_OMEGA * sin_lat
        + a2 * RH_K
            * cos_lat.powi(RH_WAVENUMBER)
            * sin_lat
            * fms_math::cos(f64::from(RH_WAVENUMBER) * lon)
}

/// Register the test fields on the climate grid and return their keys, plus
/// the key of the wind field to overlay on each.
pub fn register(world: &mut World) -> Result<Vec<(FieldKey, Option<FieldKey>)>, String> {
    let radius = world.planet.params().radius;
    let grid = &world.climate;
    let (mesh, geom) = (&grid.mesh, &grid.geometry);
    let err = |e: fms_fields::FieldError| e.to_string();

    // Solid-body rotation about an axis tilted 45°, 40 m/s at its equator.
    let axis = DVec3::new(1.0, 0.0, 1.0).normalize();
    let solid: Vec<f64> = (0..mesh.n_edges())
        .map(|e| axis.cross(geom.edge_points[e]).dot(geom.edge_normals[e]) * 40.0)
        .collect();

    // Rossby–Haurwitz wave from a vertex stream function. Using
    // u = −Δψ / l makes the discrete flow exactly non-divergent.
    let psi: Vec<f64> = mesh
        .vertices()
        .iter()
        .map(|&v| rossby_haurwitz_stream(v, radius))
        .collect();
    let rh: Vec<f64> = mesh
        .edges()
        .iter()
        .enumerate()
        .map(|(e, edge)| {
            -(psi[edge.vertices[1] as usize] - psi[edge.vertices[0] as usize])
                / (geom.edge_lengths[e] * radius)
        })
        .collect();
    // Relative vorticity by discrete circulation around each dual triangle.
    let vorticity: Vec<f64> = (0..mesh.n_vertices() as u32)
        .map(|v| {
            let circulation: f64 = mesh
                .vertex_edges(v)
                .iter()
                .zip(geom.vertex_edge_signs(v))
                .map(|(&e, s)| s * geom.dual_edge_lengths[e as usize] * radius * rh[e as usize])
                .sum();
            circulation / (geom.vertex_areas[v as usize] * radius * radius)
        })
        .collect();
    let lengths: Vec<f64> = geom
        .edge_lengths
        .iter()
        .map(|l| l * radius / 1000.0)
        .collect();

    let wind = |name, description| FieldSpec {
        description,
        persistent: false,
        quantity: Quantity::NormalComponent,
        ..FieldSpec::new(name, "m s-1", Location::CLIMATE_EDGES)
    };
    let f = &mut world.fields;
    let solid_key = f
        .register_with(
            wind("test.solid_body_wind", "Solid-body rotation, tilted axis"),
            solid,
        )
        .map_err(err)?
        .key();
    let rh_key = f
        .register_with(
            wind("test.rossby_haurwitz_wind", "Rossby-Haurwitz wave 4"),
            rh,
        )
        .map_err(err)?
        .key();
    let vort_key = f
        .register_with(
            FieldSpec {
                description: "Rossby-Haurwitz relative vorticity",
                color_scale: ColorScale::Diverging,
                persistent: false,
                ..FieldSpec::new(
                    "test.rossby_haurwitz_vorticity",
                    "s-1",
                    Location::CLIMATE_VERTICES,
                )
            },
            vorticity,
        )
        .map_err(err)?
        .key();
    let psi_key = f
        .register_with(
            FieldSpec {
                description: "Rossby-Haurwitz stream function",
                color_scale: ColorScale::Diverging,
                persistent: false,
                ..FieldSpec::new(
                    "test.rossby_haurwitz_stream",
                    "m2 s-1",
                    Location::CLIMATE_VERTICES,
                )
            },
            psi,
        )
        .map_err(err)?
        .key();
    let len_key = f
        .register_with(
            FieldSpec {
                description: "Voronoi edge length",
                persistent: false,
                ..FieldSpec::new("test.edge_length", "km", Location::CLIMATE_EDGES)
            },
            lengths,
        )
        .map_err(err)?
        .key();

    Ok(vec![
        (solid_key, Some(solid_key)),
        (rh_key, Some(rh_key)),
        (vort_key, Some(rh_key)),
        (psi_key, Some(rh_key)),
        (len_key, None),
    ])
}

#[cfg(test)]
mod tests {
    use super::*;
    use fms_geo::coords::to_lat_lon;
    use fms_sim::WorldConfig;

    fn world_at(climate_frequency: u32) -> World {
        let mut w = World::new(WorldConfig {
            terrain_frequency: 2 * climate_frequency,
            climate_frequency,
            ..WorldConfig::default()
        });
        register(&mut w).unwrap();
        w
    }

    fn world() -> World {
        world_at(16)
    }

    /// (max, rms, distance of the max-error vertex to the nearest pentagon in
    /// cell spacings) of the discrete vorticity against Williamson TC6:
    /// ζ = 2ω sinφ − K sinφ cos^R φ (R² + 3R + 2) cos Rλ.
    fn vorticity_error(climate_frequency: u32) -> (f64, f64, f64) {
        let w = world_at(climate_frequency);
        let id = w
            .fields
            .id::<f64>("test.rossby_haurwitz_vorticity")
            .unwrap();
        let zeta = w.fields.read(id);
        let r = f64::from(RH_WAVENUMBER);
        let mesh = &w.climate.mesh;
        let (mut max_err, mut sq, mut worst) = (0.0f64, 0.0, glam::DVec3::Z);
        for (v, &p) in mesh.vertices().iter().enumerate() {
            let (lat, lon) = to_lat_lon(p);
            let (s, c) = fms_math::sin_cos(lat);
            let exact = 2.0 * RH_OMEGA * s
                - RH_K
                    * s
                    * c.powi(RH_WAVENUMBER)
                    * (r * r + 3.0 * r + 2.0)
                    * fms_math::cos(r * lon);
            let err = (zeta[v] - exact).abs();
            sq += err * err;
            if err > max_err {
                (max_err, worst) = (err, p);
            }
        }
        let near = (0..mesh.n_cells() as u32)
            .filter(|&c| mesh.cell_neighbors(c).len() == 5)
            .map(|c| fms_geo::spherical::arc_length(worst, mesh.cell_centers()[c as usize]))
            .fold(f64::INFINITY, f64::min);
        (
            max_err,
            (sq / mesh.n_vertices() as f64).sqrt(),
            near / mesh.mean_spacing(),
        )
    }

    #[test]
    fn discrete_vorticity_matches_the_analytic_wave() {
        // Peak |ζ| is 7.4e-5. A sign error would give errors ~1.5e-4; truncation
        // error must shrink with resolution.
        let (coarse, coarse_rms, coarse_at) = vorticity_error(16);
        let (fine, fine_rms, fine_at) = vorticity_error(32);
        println!(
            "vorticity error: k=16 max {coarse:.3e} rms {coarse_rms:.3e} (worst {coarse_at:.1} cells from a pentagon); \
             k=32 max {fine:.3e} rms {fine_rms:.3e} (worst {fine_at:.1} cells from a pentagon)"
        );
        // As with TRiSK and Perot, the max error sits at the 12 pentagons and
        // does not shrink (Peixoto 2016); the RMS error is first order.
        assert!(coarse < 1e-5 && fine < 1e-5, "{coarse}, {fine}");
        assert!(
            coarse_at < 2.0 && fine_at < 2.0,
            "max error should be at a pentagon"
        );
        assert!(
            fine_rms < coarse_rms / 1.8,
            "should converge: {coarse_rms} -> {fine_rms}"
        );
    }

    #[test]
    fn stream_function_wind_is_non_divergent() {
        let w = world();
        let id = w.fields.id::<f64>("test.rossby_haurwitz_wind").unwrap();
        let u = w.fields.read(id);
        let (mesh, geom) = (&w.climate.mesh, &w.climate.geometry);
        for c in 0..mesh.n_cells() as u32 {
            let flux: f64 = mesh
                .cell_edges(c)
                .iter()
                .zip(geom.cell_edge_signs(c))
                .map(|(&e, s)| s * geom.edge_lengths[e as usize] * u[e as usize])
                .sum();
            assert!(flux.abs() < 1e-9, "cell {c}: {flux}");
        }
    }

    #[test]
    fn views_resample_every_element_kind() {
        let w = world();
        for name in [
            "test.solid_body_wind",
            "test.rossby_haurwitz_vorticity",
            "test.edge_length",
        ] {
            let key = w.fields.key(name).unwrap();
            let v = crate::view::view(&w, key);
            assert_eq!(v.values.len(), w.climate.mesh.n_cells(), "{name}");
            assert!(v.values.iter().all(|x| x.is_finite()));
            assert_eq!(v.vectors.is_some(), crate::view::is_vector(&w, key));
        }
        // Solid-body speed peaks at 40 m/s on the rotation equator.
        let key = w.fields.key("test.solid_body_wind").unwrap();
        let max = crate::view::view(&w, key)
            .values
            .into_iter()
            .fold(0.0, f64::max);
        assert!((max - 40.0).abs() < 1.0, "{max}");
    }
}
