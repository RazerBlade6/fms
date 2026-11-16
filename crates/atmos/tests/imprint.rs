//! The grid-imprint diagnostic on a real mesh (M2.14.2), with synthetic
//! fields: the mesh repeats under a half-turn about the axis and a mirror in
//! the equator; removing the zonal mean leaves almost nothing of a function
//! of latitude; a shared pattern with the grid's symmetry is "locked" while
//! one without it, or independent noise, is not; bumps at the pentagons
//! show there; and the sampler splits a run into alternating blocks.

use std::sync::OnceLock;

use fms_atmos::imprint::FIELDS;
use fms_atmos::williamson::{GRAVITY, OMEGA, RADIUS};
use fms_atmos::{Gas, GridImprint, ImprintAnalysis, LayeredModel, LayeredState, Levels};
use fms_grid::coords::to_lat_lon;
use fms_grid::{LloydOptions, MeshGeometry, Operators, SphereMesh};
use glam::DVec3;

/// The CLI's default for 2,562 cells: `2 · max(8, round(√(n / 10)))`.
const BANDS: usize = 32;

fn mesh() -> &'static (SphereMesh, MeshGeometry) {
    static MESH: OnceLock<(SphereMesh, MeshGeometry)> = OnceLock::new();
    MESH.get_or_init(|| {
        let (mesh, _) = SphereMesh::centroidal(16, LloydOptions::default());
        let geometry = MeshGeometry::new(&mesh);
        (mesh, geometry)
    })
}

/// A field from a function of (latitude, longitude) at the cell centres.
fn field(mesh: &SphereMesh, f: impl Fn(f64, f64) -> f64) -> Vec<f64> {
    mesh.cell_centers()
        .iter()
        .map(|&p| {
            let (lat, lon) = to_lat_lon(p);
            f(lat, lon)
        })
        .collect()
}

/// A small deterministic generator (SplitMix64).
struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9e37_79b9_7f4a_7c15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        z ^ (z >> 31)
    }

    /// Uniform noise with standard deviation `sd`.
    fn noise(&mut self, sd: f64) -> f64 {
        let u = (self.next() >> 11) as f64 / (1u64 << 53) as f64;
        (2.0 * u - 1.0) * 3f64.sqrt() * sd
    }
}

fn with_noise(values: &[f64], rng: &mut Rng, sd: f64) -> Vec<f64> {
    values.iter().map(|v| v + rng.noise(sd)).collect()
}

fn zonal(lat: f64) -> f64 {
    300.0 - 40.0 * fms_math::sin(lat).powi(2)
}

/// A temperature-like field: `(name, units, mirror sign)`.
const T: (&str, &str, f64) = ("t", "K", 1.0);

/// The mean of cos⁴φ over the sphere is 8/15, so `a cos²φ cos mλ` has an
/// RMS of `a √(4/15)`.
fn cos2_rms(a: f64) -> f64 {
    a * (4.0f64 / 15.0).sqrt()
}

#[test]
fn the_mesh_repeats_under_a_half_turn_and_a_mirror() {
    let (mesh, geometry) = mesh();
    let analysis = ImprintAnalysis::new(mesh, geometry, BANDS);
    let spacing = (4.0 * std::f64::consts::PI / mesh.n_cells() as f64).sqrt();
    eprintln!(
        "symmetry error {:.3e} rad ({:.2e} of a spacing)",
        analysis.symmetry_error(),
        analysis.symmetry_error() / spacing
    );
    // Measured: 3.0e-8 rad, 4e-7 of a grid spacing.
    assert!(analysis.symmetry_error() < 1e-6 * spacing * 10.0);
    assert!(analysis.symmetries_are_permutations());
    // The parts of any field sum to it, and each has its symmetry.
    let mut rng = Rng(3);
    let e: Vec<f64> = (0..mesh.n_cells()).map(|_| rng.noise(1.0)).collect();
    let parts = analysis.symmetry_parts(&e, -1.0);
    for (i, v) in e.iter().enumerate() {
        let sum: f64 = parts.iter().map(|p| p[i]).sum();
        assert!((sum - v).abs() < 1e-12);
    }
}

#[test]
fn a_function_of_latitude_leaves_almost_no_residual() {
    let (mesh, geometry) = mesh();
    let analysis = ImprintAnalysis::new(mesh, geometry, BANDS);
    let t = field(mesh, |lat, _| zonal(lat));
    let zero = vec![0.0; mesh.n_cells()];
    let r = analysis.analyze(T, [&t, &t], [10, 10], &zero);
    let poles: Vec<usize> = (0..mesh.n_cells())
        .filter(|&c| mesh.cell_centers()[c].z.abs() > 1.0 - 1e-9)
        .collect();
    assert_eq!(poles.len(), 2);
    let at_poles = poles
        .iter()
        .fold(0.0f64, |m, &c| m.max(r.residual[c].abs()));
    let largest = r.residual.iter().fold(0.0f64, |m, v| m.max(v.abs()));
    let strongest = r.waves.iter().fold(0.0f64, |m, v| m.max(v.abs()));
    eprintln!(
        "latitude only: rms {:.3e} K, zonal part {:.3e} K, largest {largest:.3e} K, at the poles \
         {at_poles:.3e} K, pentagons {:.3e} K, locked {:.3e} K, strongest wave {strongest:.3e} K²",
        r.rms, r.zonal, r.pentagon, r.locked
    );
    // The range is 40 K; 0.5% of it is 0.2 K. Measured: RMS 1.0e-3 K, of
    // which band means 9.5e-4 K; largest 2.7e-3 K; 7.4e-5 K at the pole
    // cells; 6.5e-4 K over the pentagons; the strongest wave 2.3e-8 K².
    // What is left is a function of latitude, so all of it is "locked"
    // (1.0e-3 K).
    assert!(largest < 0.01, "{largest}");
    assert!(at_poles < 0.001, "{at_poles}");
    assert!(r.rms < 0.002, "{}", r.rms);
    assert!(strongest < 1e-6, "{:?}", r.waves);
}

#[test]
fn a_shared_pattern_with_the_grids_symmetry_is_locked() {
    let (mesh, geometry) = mesh();
    let analysis = ImprintAnalysis::new(mesh, geometry, BANDS);
    let (a, sd) = (2.0, 0.4);
    // cos 2λ repeats under the half-turn, cos²φ under the mirror.
    let clean = field(mesh, |lat, lon| {
        zonal(lat) + a * fms_math::cos(lat).powi(2) * fms_math::cos(2.0 * lon)
    });
    let mut rng = Rng(7);
    let (ha, hb) = (
        with_noise(&clean, &mut rng, sd),
        with_noise(&clean, &mut rng, sd),
    );
    let zero = vec![0.0; mesh.n_cells()];
    let r = analysis.analyze(T, [&ha, &hb], [10, 10], &zero);
    let analytic = cos2_rms(a);
    eprintln!(
        "symmetric pattern + noise: locked {:.4} vs analytic {analytic:.4}, persistent {:.4}, \
         parts {:?}, correlation {:.3}, noise {:.4} vs {:.4}, wave 2 share {:.3}",
        r.locked,
        r.persistent,
        r.symmetric,
        r.correlation,
        r.noise,
        sd / 2f64.sqrt(),
        r.waves[1] / r.covariance
    );
    // Measured: locked 1.031 vs 1.033; covariance 1.063 K² in the
    // symmetric part against 1e-4 to 2e-3 K² in the others; correlation
    // 0.868; noise 0.2836 vs 0.2828; wave 2 holds 99.6% of the covariance.
    assert!((r.locked / analytic - 1.0).abs() < 0.1, "{}", r.locked);
    let sum: f64 = r.symmetric.iter().sum();
    assert!((sum - r.covariance).abs() < 1e-9 * r.covariance.abs());
    // Half the difference of two halves with noise `sd` each has `sd / √2`.
    assert!(
        (r.noise / (sd / 2f64.sqrt()) - 1.0).abs() < 0.1,
        "{}",
        r.noise
    );
}

#[test]
fn a_shared_pattern_without_the_grids_symmetry_is_not_locked() {
    let (mesh, geometry) = mesh();
    let analysis = ImprintAnalysis::new(mesh, geometry, BANDS);
    let a = 2.0;
    let zero = vec![0.0; mesh.n_cells()];
    // Odd under the half-turn; odd under the mirror.
    for (name, pattern) in [
        (
            "cos²φ cos 3λ",
            Box::new(|lat: f64, lon: f64| fms_math::cos(lat).powi(2) * fms_math::cos(3.0 * lon))
                as Box<dyn Fn(f64, f64) -> f64>,
        ),
        (
            "sin φ cos²φ cos 2λ",
            Box::new(|lat: f64, lon: f64| {
                fms_math::sin(lat) * fms_math::cos(lat).powi(2) * fms_math::cos(2.0 * lon)
            }),
        ),
    ] {
        let both = field(mesh, |lat, lon| zonal(lat) + a * pattern(lat, lon));
        let r = analysis.analyze(T, [&both, &both], [10, 10], &zero);
        eprintln!(
            "{name}: persistent {:.4}, locked {:.3e}, parts {:?}",
            r.persistent, r.locked, r.symmetric
        );
        // Measured: persistent 1.033 (cos 3λ) and 0.389 (sin φ, whose RMS is
        // a √(8/105) = 0.390); both wholly in their own part, so locked is 0.
        assert!(r.persistent > 0.3, "{}", r.persistent);
        assert!(r.locked < 0.02 * r.persistent, "{}", r.locked);
    }
}

#[test]
fn independent_noise_shares_no_pattern() {
    let (mesh, geometry) = mesh();
    let analysis = ImprintAnalysis::new(mesh, geometry, BANDS);
    let sd = 1.0;
    let base = field(mesh, |lat, _| zonal(lat));
    let mut rng = Rng(42);
    let (ha, hb) = (
        with_noise(&base, &mut rng, sd),
        with_noise(&base, &mut rng, sd),
    );
    let zero = vec![0.0; mesh.n_cells()];
    let r = analysis.analyze(T, [&ha, &hb], [10, 10], &zero);
    eprintln!(
        "noise only: correlation {:.4}, persistent {:.4}, locked {:.4}, noise {:.4}, parts {:?}",
        r.correlation, r.persistent, r.locked, r.noise, r.symmetric
    );
    // Measured: correlation 0.012, persistent 0.107, locked 0 (the
    // symmetric part's covariance came out −0.007 K², the others +0.004 to
    // +0.009 K²), noise 0.696.
    assert!(r.correlation.abs() < 0.1, "{}", r.correlation);
    assert!(r.persistent < 0.35 * r.noise);
    assert!(r.locked < 0.35 * r.noise);
}

#[test]
fn bumps_on_the_pentagons_show_there() {
    let (mesh, geometry) = mesh();
    let analysis = ImprintAnalysis::new(mesh, geometry, BANDS);
    assert_eq!(analysis.pentagons().len(), 12);
    assert_eq!(analysis.ring().len(), 12 * 6);
    let mut bumps = field(mesh, |lat, _| zonal(lat));
    for &c in analysis.pentagons() {
        bumps[c] += 1.0;
    }
    let zero = vec![0.0; mesh.n_cells()];
    let r = analysis.analyze(T, [&bumps, &bumps], [10, 10], &zero);
    eprintln!(
        "pentagon bumps: pentagon {:.4}, ring {:.4}, persistent {:.4}, locked {:.4}, parts {:?}",
        r.pentagon, r.ring, r.persistent, r.locked, r.symmetric
    );
    // Measured: pentagon 0.980, ring 0.400 (√(1/6) = 0.41), persistent
    // 0.0596 (√(12/2562) = 0.068 less what the zonal-mean fit takes), and
    // locked 0.0596: the 12 pentagons have the grid's symmetry, so all of
    // the covariance is in the symmetric part.
    assert!(r.pentagon > 0.9, "{}", r.pentagon);
    assert!(r.pentagon > 5.0 * r.persistent);
    assert!(r.ring < r.pentagon && r.ring > r.persistent);
    assert!((r.locked / r.persistent - 1.0).abs() < 0.01);
}

#[test]
fn the_sampler_alternates_blocks_and_keeps_means_and_spread() {
    let (mesh, geometry) = mesh();
    let model = LayeredModel::new(
        Operators::new(mesh, geometry, RADIUS),
        GRAVITY,
        Gas::EARTH_AIR,
        Levels::two(0.5),
        mesh.vertices().iter().map(|p| 2.0 * OMEGA * p.z).collect(),
        vec![0.0; mesh.n_cells()],
    );
    let (upper, lower) = (0, 1);
    let n = mesh.n_cells();
    // Seven samples in blocks of two: A A B B A A B.
    let ps = [1.00e5, 1.01e5, 0.99e5, 1.02e5, 0.98e5, 1.00e5, 1.03e5];
    let temperature = [280.0, 282.0, 284.0, 286.0, 288.0, 290.0, 292.0];
    let wind = [10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 22.0];
    // Solid-body eastward wind `U cos φ` in the upper layer.
    let solid = |speed: f64| -> Vec<f64> {
        geometry
            .edge_points
            .iter()
            .zip(&geometry.edge_normals)
            .map(|(&p, &normal)| {
                let east = DVec3::Z.cross(p);
                speed * east.dot(normal)
            })
            .collect()
    };
    let mut imprint = GridImprint::new(&model, upper, lower, 2, BANDS);
    for t in 0..ps.len() {
        let exner = model.layer_exner(ps[t])[lower];
        let state = LayeredState::new(
            vec![ps[t]; n],
            vec![
                solid(wind[t]),
                vec![0.0; model.operators().mesh().n_edges()],
            ],
            vec![vec![300.0; n], vec![temperature[t] / exner; n]],
        );
        imprint.sample(&model, &state);
    }
    assert_eq!(imprint.samples(), 7);
    assert_eq!(imprint.complete_blocks(), 3);
    assert_eq!(imprint.counts(), [4, 3]);
    let halves = |v: &[f64]| {
        let a = [0, 1, 4, 5].iter().map(|&i| v[i]).sum::<f64>() / 4.0;
        let b = [2, 3, 6].iter().map(|&i| v[i]).sum::<f64>() / 3.0;
        [a, b]
    };
    let spread = |v: &[f64]| {
        let mean = v.iter().sum::<f64>() / v.len() as f64;
        (v.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / v.len() as f64).sqrt()
    };
    for (f, values, tolerance) in [(0, &ps, 1e-6), (1, &temperature, 1e-9)] {
        let [a, b] = imprint.means(f);
        let expected = halves(values);
        assert!(a.iter().all(|x| (x - expected[0]).abs() < tolerance));
        assert!(b.iter().all(|x| (x - expected[1]).abs() < tolerance));
        let sd = spread(values);
        assert!(
            imprint
                .transient(f)
                .iter()
                .all(|x| (x - sd).abs() < 1e-6 * sd),
            "{}",
            FIELDS[f].0
        );
    }
    // The upper layer's eastward wind (reconstructed: within a few percent
    // of U cos φ); the lower layer's northward wind is zero.
    let [a, _] = imprint.means(2);
    let expected = halves(&wind)[0];
    for (&p, &u) in mesh.cell_centers().iter().zip(&a) {
        let cos = DVec3::Z.cross(p).length();
        assert!((u - expected * cos).abs() < 0.05 * expected, "{u}");
    }
    assert!(imprint.means(3)[0].iter().all(|v| v.abs() < 1e-9));
    let results = imprint.results();
    assert_eq!(results.len(), FIELDS.len());
    assert_eq!(results[0].name, "surface_pressure");
    // Uniform fields leave no residual.
    assert!(results[0].rms < 1e-6 && results[1].rms < 1e-9);
}
