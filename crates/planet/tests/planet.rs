//! Validation, derived quantities, orbit and insolation against analytic results.

use std::f64::consts::{FRAC_PI_2, PI, TAU};

use fms_geo::coords::from_lat_lon;
use fms_planet::{
    Adjustment, Planet, PlanetParams, eccentric_anomaly, eccentric_from_true, true_anomaly,
};

/// A planet with a circular orbit (analytic cases are cleanest there).
fn circular(obliquity_deg: f64) -> Planet {
    Planet::new(PlanetParams {
        eccentricity: 0.0,
        obliquity_deg,
        ..PlanetParams::earth()
    })
    .0
}

fn close(a: f64, b: f64, tol: f64) -> bool {
    (a - b).abs() <= tol * b.abs().max(1.0)
}

#[test]
fn earth_needs_no_adjustments() {
    let (_, adjustments) = Planet::new(PlanetParams::earth());
    assert!(adjustments.is_empty(), "{adjustments:?}");
}

#[test]
fn out_of_range_values_are_clamped_and_reported() {
    let params = PlanetParams {
        eccentricity: 0.9,
        obliquity_deg: f64::NAN,
        perihelion_longitude_deg: -90.0,
        ocean_fraction: 1.5,
        ..PlanetParams::earth()
    };
    let (planet, adjustments) = Planet::new(params);
    let p = planet.params();
    assert_eq!(p.eccentricity, 0.3);
    assert_eq!(p.obliquity_deg, 23.44);
    assert_eq!(p.perihelion_longitude_deg, 270.0);
    assert_eq!(p.ocean_fraction, 1.0);
    assert!(adjustments.contains(&Adjustment::Clamped {
        parameter: "eccentricity",
        from: 0.9,
        to: 0.3
    }));
    assert!(adjustments.contains(&Adjustment::Replaced {
        parameter: "obliquity_deg",
        to: 23.44
    }));
    assert!(adjustments.contains(&Adjustment::Wrapped {
        from: -90.0,
        to: 270.0
    }));
    assert_eq!(adjustments.len(), 4);
    // Messages are human-readable.
    assert!(
        adjustments[0].to_string().contains("supported range")
            || adjustments[0].to_string().contains("finite")
    );
}

#[test]
fn extreme_flux_is_flagged_not_clamped() {
    let (planet, adjustments) = Planet::new(PlanetParams {
        semi_major_axis_au: 5.0,
        ..PlanetParams::earth()
    });
    assert!(close(planet.mean_flux(), 1361.0 / 25.0, 1e-12));
    assert!(matches!(
        adjustments.as_slice(),
        [Adjustment::ExtremeFlux { .. }]
    ));
}

#[test]
fn earth_derived_quantities() {
    let earth = Planet::earth();
    // Sidereal day 23 h 56 min 4.1 s.
    assert!(
        (earth.sidereal_day() - 86_164.09).abs() < 0.05,
        "{}",
        earth.sidereal_day()
    );
    assert!(close(earth.rotation_rate(), 7.292_115e-5, 1e-6));
    assert!(close(earth.coriolis(FRAC_PI_2 / 3.0), 7.292e-5, 1e-3)); // f at 30° = Ω
    assert_eq!(earth.coriolis(0.0), 0.0);
    assert!(close(
        earth.coriolis_at(from_lat_lon(0.7, 1.0)),
        earth.coriolis(0.7),
        1e-12
    ));
    assert!(close(earth.mean_flux(), 1361.0, 1e-12));
}

#[test]
fn kepler_solver_is_accurate() {
    for &e in &[0.0, 0.0167, 0.1, 0.3, 0.9] {
        for i in 0..100 {
            let m = TAU * f64::from(i) / 100.0;
            let ecc = eccentric_anomaly(m, e);
            assert!(
                (ecc - e * fms_math::sin(ecc) - m).abs() < 1e-12,
                "e={e} M={m}"
            );
            let nu = true_anomaly(ecc, e);
            // Round trip true → eccentric anomaly, compared as angles.
            let back = eccentric_from_true(nu, e);
            assert!(
                fms_math::sin(back - ecc).abs() < 1e-9 && fms_math::cos(back - ecc) > 0.0,
                "e={e} M={m}"
            );
        }
    }
}

#[test]
fn equinox_solstice_and_perihelion() {
    let earth = Planet::earth();
    let year = earth.year_seconds();
    let s0 = earth.orbit_state(0.0);
    assert!(s0.solar_longitude.abs() < 1e-12 || (s0.solar_longitude - TAU).abs() < 1e-12);
    assert!(s0.declination.abs() < 1e-12);

    // Find the northern summer solstice: maximum declination equals the obliquity.
    let max_dec = (0..10_000)
        .map(|i| {
            earth
                .orbit_state(year * f64::from(i) / 10_000.0)
                .declination
        })
        .fold(f64::NEG_INFINITY, f64::max);
    assert!((max_dec - 23.44f64.to_radians()).abs() < 1e-6);

    // Perihelion: closest approach, at solar longitude 282.9°, distance 1 − e.
    let (t_peri, r_min) = (0..100_000)
        .map(|i| {
            let t = year * f64::from(i) / 100_000.0;
            (t, earth.orbit_state(t).distance_ratio)
        })
        .fold((0.0, f64::INFINITY), |a, b| if b.1 < a.1 { b } else { a });
    assert!((r_min - (1.0 - 0.0167)).abs() < 1e-8);
    let lon = earth.orbit_state(t_peri).solar_longitude.to_degrees();
    assert!((lon - 282.9).abs() < 0.01, "{lon}");
    // Earth's perihelion is ~13 days after the December solstice (~286 days after
    // the March equinox).
    let days_after_equinox = t_peri / earth.solar_day();
    assert!(
        (days_after_equinox - 288.0).abs() < 4.0,
        "{days_after_equinox}"
    );
}

#[test]
fn eccentric_orbits_have_unequal_seasons() {
    // Kepler's second law: the half-year containing perihelion is shorter.
    let earth = Planet::earth();
    let year = earth.year_seconds();
    let n = 100_000;
    let northern_summer_half = (0..n)
        .filter(|&i| {
            earth
                .orbit_state(year * f64::from(i) / f64::from(n))
                .declination
                > 0.0
        })
        .count() as f64
        / f64::from(n);
    // On Earth, March–September equinox spans ~186 of 365.24 days.
    assert!(
        (northern_summer_half * 365.2422 - 186.4).abs() < 1.0,
        "{}",
        northern_summer_half * 365.2422
    );
}

#[test]
fn daily_mean_insolation_analytic_values() {
    let p = circular(23.44);
    let s0 = p.mean_flux();
    let quarter = p.year_seconds() / 4.0;
    // Equator at equinox: S₀/π.
    assert!(close(p.daily_mean_insolation(0.0, 0.0), s0 / PI, 1e-12));
    // North pole at the summer solstice: S₀ sin ε (Sun circles at elevation ε).
    assert!(close(
        p.daily_mean_insolation(FRAC_PI_2, quarter),
        s0 * fms_math::sin(23.44f64.to_radians()),
        1e-9
    ));
    // Polar night.
    assert_eq!(p.daily_mean_insolation(-FRAC_PI_2, quarter), 0.0);
    assert_eq!(
        p.daily_mean_insolation((-80.0f64).to_radians(), quarter),
        0.0
    );
    // No tilt: no seasons, so daily-mean insolation is the same all year.
    let flat = circular(0.0);
    for i in 0..8 {
        let t = flat.year_seconds() * f64::from(i) / 8.0;
        assert!(close(
            flat.daily_mean_insolation(0.5, t),
            flat.daily_mean_insolation(0.5, 0.0),
            1e-9
        ));
    }
}

/// Area-weighted global mean of a function of latitude (midpoint rule in sin φ).
fn global_mean(f: impl Fn(f64) -> f64) -> f64 {
    let n = 20_000;
    (0..n)
        .map(|i| {
            let mu = -1.0 + 2.0 * (f64::from(i) + 0.5) / f64::from(n);
            f(fms_math::asin(mu))
        })
        .sum::<f64>()
        / f64::from(n)
}

#[test]
fn global_mean_insolation_is_a_quarter_of_the_flux() {
    let earth = Planet::earth();
    for i in 0..6 {
        let t = earth.year_seconds() * f64::from(i) / 6.0;
        let r = earth.orbit_state(t).distance_ratio;
        let expected = earth.mean_flux() / 4.0 / (r * r);
        let mean = global_mean(|lat| earth.daily_mean_insolation(lat, t));
        assert!(close(mean, expected, 1e-6), "t={t}: {mean} vs {expected}");
    }
}

#[test]
fn instantaneous_insolation_integrates_to_daily_mean() {
    let earth = Planet::earth();
    let t0 = earth.year_seconds() * 0.3; // arbitrary day
    let day = earth.solar_day();
    for &lat_deg in &[-70.0, -30.0, 0.0, 45.0, 75.0] {
        let lat = f64::to_radians(lat_deg);
        let n = 4000;
        // Freeze the orbit over the day: daily_mean_insolation assumes that too.
        let mean: f64 = (0..n)
            .map(|i| {
                let t = t0 + day * (f64::from(i) + 0.5) / f64::from(n);
                let sun = earth.sun_direction(t);
                let state = earth.orbit_state(t0);
                let p = from_lat_lon(lat, 1.3);
                // Use the day-start declination for an apples-to-apples comparison.
                let sun = from_lat_lon(state.declination, fms_math::atan2(sun.y, sun.x));
                earth.mean_flux() / state.distance_ratio.powi(2) * p.dot(sun).max(0.0)
            })
            .sum::<f64>()
            / f64::from(n);
        let expected = earth.daily_mean_insolation(lat, t0);
        assert!(
            (mean - expected).abs() < 1e-3 * earth.mean_flux(),
            "lat {lat_deg}: {mean} vs {expected}"
        );
    }
}

#[test]
fn sun_is_overhead_at_noon_on_the_meridian() {
    let p = circular(0.0);
    let day = p.solar_day();
    // Equinox noon at longitude 0: Sun straight overhead at (0, 0).
    assert!((p.sun_direction(0.5 * day) - from_lat_lon(0.0, 0.0)).length() < 1e-9);
    // Six hours later the subsolar point has moved 90° west.
    assert!((p.sun_direction(0.75 * day) - from_lat_lon(0.0, -FRAC_PI_2)).length() < 1e-6);
    // Midnight: Sun at the antimeridian, so no light at (0, 0).
    assert_eq!(p.insolation(from_lat_lon(0.0, 0.0), 0.0), 0.0);
    assert!(close(
        p.insolation(from_lat_lon(0.0, PI), 0.0),
        p.mean_flux(),
        1e-9
    ));
}

#[test]
fn annual_mean_insolation_is_hemispherically_symmetric() {
    // True even for eccentric orbits: Kepler's second law balances the faster,
    // closer perihelion passage exactly.
    let p = Planet::new(PlanetParams {
        eccentricity: 0.2,
        ..PlanetParams::earth()
    })
    .0;
    let year = p.year_seconds();
    let n = 20_000;
    let annual = |lat: f64| {
        (0..n)
            .map(|i| p.daily_mean_insolation(lat, year * (f64::from(i) + 0.5) / f64::from(n)))
            .sum::<f64>()
            / f64::from(n)
    };
    for &lat_deg in &[10.0, 40.0, 70.0, 89.0] {
        let lat = f64::to_radians(lat_deg);
        let (north, south) = (annual(lat), annual(-lat));
        assert!(
            (north - south).abs() < 1e-3 * north,
            "{lat_deg}: {north} vs {south}"
        );
    }
}

#[test]
fn insolation_into_matches_pointwise() {
    let earth = Planet::earth();
    let points: Vec<_> = (0..500)
        .map(|i| from_lat_lon(fms_math::sin(f64::from(i) * 0.37), f64::from(i) * 1.1))
        .collect();
    let mut out = vec![0.0; points.len()];
    let t = 12_345_678.0;
    earth.insolation_into(&points, t, &mut out);
    for (p, o) in points.iter().zip(&out) {
        assert_eq!(*o, earth.insolation(*p, t));
    }
}

#[test]
fn earth_seasons_have_their_real_lengths() {
    // Present-day northern seasons: spring 92.76, summer 93.65, autumn 89.84,
    // winter 88.99 days.
    let lengths = Planet::earth().season_lengths_days();
    for (got, want) in lengths.iter().zip([92.76, 93.65, 89.84, 88.99]) {
        assert!((got - want).abs() < 0.3, "{lengths:?}");
    }
    assert!((lengths.iter().sum::<f64>() - 365.2422).abs() < 1e-9);
}

#[test]
fn circular_orbits_have_equal_seasons() {
    for len in circular(23.44).season_lengths_days() {
        assert!((len - 365.2422 / 4.0).abs() < 1e-9);
    }
}

#[test]
fn seasons_switch_at_equinoxes_and_solstices() {
    use fms_planet::{Hemisphere, Season};
    let earth = Planet::earth();
    let starts = earth.season_starts();
    assert_eq!(starts[0], 0.0);
    for (q, &t) in starts.iter().enumerate() {
        // Solar longitude at each start is exactly a multiple of 90°.
        let lon = earth.orbit_state(t).solar_longitude;
        let target = q as f64 * FRAC_PI_2;
        assert!(fms_math::sin(lon - target).abs() < 1e-9 && fms_math::cos(lon - target) > 0.0);
        let just_after = t + 60.0;
        assert_eq!(earth.season(just_after, Hemisphere::North), Season::ALL[q]);
        assert_eq!(
            earth.season(just_after, Hemisphere::South),
            Season::ALL[(q + 2) % 4]
        );
    }
    // Maximum declination falls in the summer solstice's instant.
    let solstice = earth.orbit_state(starts[1]).declination;
    assert!((solstice - 23.44f64.to_radians()).abs() < 1e-9);
}
