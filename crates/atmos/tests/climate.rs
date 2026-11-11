//! Radiation and moisture (M2 steps 8–9) in 3-D: the two-layer model over a
//! slab aquaplanet with Earth's orbit. The planet's energy and water budgets
//! close, and (slow, run with `--ignored`) the equilibrium climate is that of
//! an ice-free, Earth-like aquaplanet.

use fms_atmos::column_physics::{
    ApeSst, CLIMATE_LEVELS, LATENT_HEAT, LandParams, OceanSpinUp, initial_state,
    initial_surface_temperature_by_latitude,
};
use fms_atmos::williamson::{GRAVITY, RADIUS};
use fms_atmos::{
    ColumnPhysics, EARTH_ROTATION, EddyParams, EkmanParams, Forcing, Gas, HyperdiffusionOptions,
    LayeredModel, LayeredSemiImplicit, LayeredSemiImplicitOptions, Levels, OceanParams,
    PhysicsParams,
};
use fms_grid::{LatitudeBands, LloydOptions, MeshGeometry, Operators, SphereMesh};
use fms_planet::{Planet, PlanetParams};

const DAY: f64 = 86_400.0;
/// A callback run after each step.
type EachStep<'a> = &'a mut dyn FnMut(&LayeredSemiImplicit, &MeshGeometry, f64);
/// No terrain: a cell's sub-grid height spread and mean height (m).
const FLAT: fn(glam::DVec3) -> (f64, f64) = |_| (0.0, 0.0);

/// Run `days` on frequency `k` over a 10 m slab, which settles in a few
/// years, from the single-column equilibrium, calling `each` after every
/// step.
fn run(k: u32, days: f64, each: impl FnMut(&LayeredSemiImplicit, &MeshGeometry, f64)) {
    let params = PhysicsParams {
        ocean: OceanParams::slab(10.0),
        ..PhysicsParams::default()
    };
    let levels = Levels::new(&CLIMATE_LEVELS);
    run_over(params, k, days, levels, |_| 0.0, FLAT, each);
}

/// Land fraction of a continent within 45° of the prime meridian, with
/// partly land cells along its coasts.
fn continent_fraction(p: glam::DVec3) -> f64 {
    let lon = fms_math::atan2(p.y, p.x).abs().to_degrees();
    ((45.0 - lon) / 10.0 + 0.5).clamp(0.0, 1.0)
}

/// [`run`] on other layers, land and relief, over the default ocean.
fn run_with(
    k: u32,
    days: f64,
    levels: Levels,
    land: impl Fn(glam::DVec3) -> f64,
    relief: impl Fn(glam::DVec3) -> (f64, f64),
    each: impl FnMut(&LayeredSemiImplicit, &MeshGeometry, f64),
) {
    let params = PhysicsParams::default();
    run_over(params, k, days, levels, land, relief, each);
}

/// [`run_with`] with any physics.
fn run_over(
    params: PhysicsParams,
    k: u32,
    days: f64,
    levels: Levels,
    land: impl Fn(glam::DVec3) -> f64,
    relief: impl Fn(glam::DVec3) -> (f64, f64),
    each: impl FnMut(&LayeredSemiImplicit, &MeshGeometry, f64),
) {
    run_on(
        &Planet::earth(),
        params,
        k,
        days,
        levels,
        land,
        relief,
        each,
    );
}

/// [`run_over`] on any planet.
#[allow(clippy::too_many_arguments)]
fn run_on(
    planet: &Planet,
    params: PhysicsParams,
    k: u32,
    days: f64,
    levels: Levels,
    land: impl Fn(glam::DVec3) -> f64,
    relief: impl Fn(glam::DVec3) -> (f64, f64),
    mut each: impl FnMut(&LayeredSemiImplicit, &MeshGeometry, f64),
) {
    let (mesh, _) = SphereMesh::centroidal(k, LloydOptions::default());
    let geometry = MeshGeometry::new(&mesh);
    let model = LayeredModel::new(
        Operators::new(&mesh, &geometry, RADIUS),
        GRAVITY,
        Gas::EARTH_AIR,
        levels,
        mesh.vertices()
            .iter()
            .map(|p| 2.0 * planet.rotation_rate() * p.z)
            .collect(),
        mesh.cell_centers()
            .iter()
            .map(|&p| GRAVITY * relief(p).1)
            .collect(),
    );
    let (state, ts) = initial_state(&model, planet, &params);
    // Radiatively driven jets and storms reach ~100 m/s.
    let dt = LayeredSemiImplicit::stable_dt(&model, 120.0, 0.6);
    let dt = DAY / (DAY / dt).ceil();
    let diffusion = HyperdiffusionOptions::strongest_for_step(model.operators(), 4, dt);
    let physics = ColumnPhysics::with_land(
        &model,
        planet,
        params,
        0.0,
        &state,
        vec![ts; mesh.n_cells()],
        mesh.cell_centers().iter().map(|&p| land(p)).collect(),
    )
    .with_relief(
        mesh.cell_centers().iter().map(|&p| relief(p).0).collect(),
        &mesh
            .cell_centers()
            .iter()
            .map(|&p| relief(p).1)
            .collect::<Vec<_>>(),
    );
    let model = model
        .with_hyperdiffusion(diffusion)
        .with_dissipative_heating()
        .with_forcing(Box::new(physics));
    let mut si =
        LayeredSemiImplicit::new(model, state, dt, LayeredSemiImplicitOptions::default(), 0.5);
    while (si.steps() as f64) * dt < days * DAY {
        si.step();
        si.check().expect("the moist climate blew up");
        each(&si, &geometry, dt);
    }
}

/// FNV-1a hash of every prognostic field's bits.
fn state_hash(s: &fms_atmos::LayeredState) -> u64 {
    let mut h = 0xcbf2_9ce4_8422_2325u64;
    let fields = std::iter::once(&s.ps)
        .chain(&s.u)
        .chain(&s.theta_ps)
        .chain(&s.q_ps);
    for x in fields.flatten() {
        h = (h ^ x.to_bits()).wrapping_mul(0x0000_0100_0000_01b3);
    }
    h
}

#[test]
fn a_world_without_land_is_the_aquaplanet_bit_for_bit() {
    // Pinned from the aquaplanet before land tiles existed (M2 step 10.2.2).
    // If this fails, first ask whether the change was meant to alter the
    // ocean-only climate; only then paste in the new value.
    const GOLDEN: u64 = 0xe70e_3500_47bd_d6b0;
    let mut hash = 0;
    run(8, 5.0, |si, _, _| hash = state_hash(si.state()));
    assert_eq!(hash, GOLDEN, "new hash: {hash:#018x}");
}

#[test]
fn the_planets_energy_and_water_budgets_close() {
    budgets_close(16, 90.0, |_| 0.0, FLAT, |_, _, _| ());
}

#[test]
fn the_planets_energy_and_water_budgets_close_with_land_tiles() {
    budgets_close(16, 90.0, continent_fraction, FLAT, |_, _, _| ());
}

/// Run `days` at frequency `k`, calling `each` after every step, and check
/// the planet's budgets over the run.
fn budgets_close(
    k: u32,
    days: f64,
    land: impl Fn(glam::DVec3) -> f64,
    relief: impl Fn(glam::DVec3) -> (f64, f64),
    mut each: impl FnMut(&LayeredSemiImplicit, &MeshGeometry, f64),
) {
    // What the atmosphere (enthalpy and latent energy) and the ocean and
    // land gain equals sunlight absorbed minus longwave emitted (friction's
    // and diffusion's lost kinetic energy comes back as heat); only the
    // time-stepping error is left. A 50 m ocean starting at one uniform
    // temperature is still far from balance at the end of these runs, which
    // makes that error larger than the old 10 m slab's: −0.042 W/m² over 90
    // days (10 m: −0.022) and −0.064 over 425 days at 640 cells (10 m:
    // −0.024). It is not the ocean's transport: a 50 m slab with none gives
    // −0.064 too, and halving the step takes that to −0.022. The
    // atmosphere's water changes by evaporation minus precipitation (and
    // the time filter's trace), to round-off.
    let mut first = None;
    let mut last = None;
    let levels = Levels::new(&CLIMATE_LEVELS);
    run_with(k, days, levels, land, relief, |si, geometry, dt| {
        each(si, geometry, dt);
        let budget = si.model().forcing().unwrap().budget().unwrap();
        let row = (
            si.model().diagnostics(si.state()),
            si.energy_budget(),
            budget,
            si.steps() as f64 * dt,
            si.water(),
            si.water_budget(),
        );
        if first.is_none() {
            first = Some(row);
        }
        last = Some(row);
    });
    let ((d0, e0, b0, t0, w0, wb0), (d1, e1, b1, t1, w1, wb1)) = (first.unwrap(), last.unwrap());
    let area = 4.0 * std::f64::consts::PI * RADIUS * RADIUS;
    let wm2 = |j: f64| j / (t1 - t0) / area;
    let toa = (b1.absorbed - b0.absorbed) - (b1.outgoing - b0.outgoing);
    let uptake = (d1.energy - d0.energy) + (b1.stored - b0.stored) + LATENT_HEAT * (w1 - w0);
    // Friction's and diffusion's lost kinetic energy comes back as heat.
    let removed = (e1.forcing.kinetic - e0.forcing.kinetic)
        + (b1.dissipated - b0.dissipated)
        + (e1.diffusion.total() - e0.diffusion.total());
    let residual = toa - uptake + removed;
    eprintln!(
        "{days} days: net in at the top {:+.3} W/m², uptake {:+.3}, kinetic energy removed {:.3}, residual {:+.5}",
        wm2(toa),
        wm2(uptake),
        -wm2(removed),
        wm2(residual)
    );
    assert!(
        wm2(residual).abs() < 0.08,
        "residual {} W/m²",
        wm2(residual)
    );
    assert!((d1.mass / d0.mass - 1.0).abs() < 1e-13);

    let exchanged = (b1.evaporated - b0.evaporated) - (b1.precipitated - b0.precipitated);
    let (filter, fixer) = (wb1.filter - wb0.filter, wb1.fixer - wb0.fixer);
    let water_residual = (w1 - w0) - exchanged - filter - fixer;
    let mm_day = |kg: f64| kg / area / (t1 - t0) * DAY;
    eprintln!(
        "water: {:.1} kg/m², evaporation − precipitation {:+.4} mm/day, filter {:+.1e}, fixer {:+.1e}, residual {:+.1e} of the water",
        w1 / area,
        mm_day(exchanged),
        filter / w1,
        fixer / w1,
        water_residual / w1
    );
    assert!(
        water_residual.abs() < 1e-12 * w1,
        "water residual {water_residual:e} kg"
    );
    // The filter's trace wanders during spin-up, fully accounted for above.
    // Over the 50 m ocean: measured 6.1e-3 of the water at day 90 (the 10 m
    // slab's 1.6e-3), where the air has dried to 10 kg/m². Its sum grows with
    // the run: 2.1e-2 over 425 days at 640 cells (10 m: 7.3e-3).
    let filter_limit = 1e-2 * (days / 90.0).max(1.0);
    assert!(filter.abs() < filter_limit * w1 && fixer.abs() < 1e-9 * w1);

    // The land's share: atmosphere, soil, snow and what ran off together
    // change only by what the ocean tiles exchange with the air.
    let (soil, runoff) = (b1.soil_water - b0.soil_water, b1.runoff - b0.runoff);
    let snow = b1.snow - b0.snow;
    let ocean = b1.ocean_exchange - b0.ocean_exchange;
    let total_residual = (w1 - w0) + soil + snow + runoff - ocean - filter - fixer;
    eprintln!(
        "soil water {:+.4} mm/day, snow {:+.4} mm/day, runoff {:.4} mm/day, residual {:+.1e} of the water",
        mm_day(soil),
        mm_day(snow),
        mm_day(runoff),
        total_residual / w1
    );
    assert!(
        total_residual.abs() < 1e-11 * w1,
        "total residual {total_residual:e} kg"
    );
}

#[test]
fn an_ocean_jump_moves_the_sea_and_stays_out_of_the_energy_budget() {
    let (planet, _) = Planet::new(PlanetParams {
        year_days: 20.0,
        ..PlanetParams::earth()
    });
    let params = |jumps| PhysicsParams {
        spin_up: OceanSpinUp {
            jumps,
            ..OceanSpinUp::default()
        },
        ..PhysicsParams::default()
    };
    // A 20-day year, so the jump comes at the end of day 40.
    let with = budget_run(&planet, params(1), 45.0, |_| 0.0);
    let without = budget_run(&planet, params(0), 45.0, |_| 0.0);
    eprintln!(
        "residual {:+.5} W/m², without the jump {:+.5}",
        with.energy, without.energy
    );
    assert!(
        with.largest_change > 0.2,
        "largest step change {} K",
        with.largest_change
    );
    assert!(
        (40.0..40.1).contains(&with.largest_day),
        "at day {}",
        with.largest_day
    );
    // The jump adds tens of W/m² of heat over this run; none of it may show.
    assert!(
        (with.energy - without.energy).abs() < 0.03,
        "{} against {}",
        with.energy,
        without.energy
    );
}

#[test]
fn energy_and_water_close_through_freezing_and_melting() {
    // A dim planet with fast seasons and a shallow sea, so that it freezes
    // and thaws within the run.
    let (planet, _) = Planet::new(PlanetParams {
        year_days: 40.0,
        luminosity: 0.9,
        ..PlanetParams::earth()
    });
    // One shallow layer, which can't carry currents: still water.
    let params = PhysicsParams {
        ocean: OceanParams {
            mixed_depth: 5.0,
            deep_depth: 0.0,
            ..OceanParams::still()
        },
        ..PhysicsParams::default()
    };
    let run = budget_run(&planet, params, 90.0, |_| 0.0);
    eprintln!(
        "energy residual {:+.5} W/m², water residual {:+.1e}, ice cover {:.3} to {:.3} of the planet",
        run.energy, run.water, run.least_ice, run.most_ice
    );
    assert!(
        run.energy.abs() < 0.1,
        "energy residual {} W/m²",
        run.energy
    );
    assert!(run.water.abs() < 1e-11, "water residual {}", run.water);
    assert!(run.most_ice > 0.05, "ice cover {}", run.most_ice);
    assert!(
        run.most_ice > run.least_ice + 0.005,
        "the ice never changed"
    );
}

#[test]
fn energy_and_water_close_as_land_becomes_an_ice_sheet() {
    // Polar land on a dim planet whose short year lets cells turn into ice
    // sheets within the run.
    let (planet, _) = Planet::new(PlanetParams {
        year_days: 40.0,
        luminosity: 0.85,
        ..PlanetParams::earth()
    });
    let params = PhysicsParams {
        land: LandParams {
            snow_covering: 2.0,
            ice_sheet_years: 0.5,
            ice_sheet_start_temperature: 0.0,
            ..LandParams::default()
        },
        ..PhysicsParams::default()
    };
    let polar = |p: glam::DVec3| if p.z.abs() > 0.75 { 1.0 } else { 0.0 };
    let run = budget_run(&planet, params, 90.0, polar);
    eprintln!(
        "energy residual {:+.5} W/m², water residual {:+.1e}, ice sheets on {:.3} of the planet",
        run.energy, run.water, run.most_sheet
    );
    assert!(
        run.energy.abs() < 0.1,
        "energy residual {} W/m²",
        run.energy
    );
    assert!(run.water.abs() < 1e-11, "water residual {}", run.water);
    assert!(run.most_sheet > 0.02, "ice sheets {}", run.most_sheet);
}

#[test]
fn energy_and_water_close_with_a_seeded_ice_sheet() {
    // Everything starts colder than the threshold, so all the land starts
    // as ice sheet, and a pack this thin can melt away within the run.
    let (planet, _) = Planet::new(PlanetParams {
        year_days: 40.0,
        ..PlanetParams::earth()
    });
    let params = PhysicsParams {
        land: LandParams {
            ice_sheet_start_temperature: 400.0,
            ice_sheet_start_pack: 30.0,
            ..LandParams::default()
        },
        ..PhysicsParams::default()
    };
    let polar = |p: glam::DVec3| if p.z.abs() > 0.75 { 1.0 } else { 0.0 };
    let run = budget_run(&planet, params, 45.0, polar);
    eprintln!(
        "energy residual {:+.5} W/m², water residual {:+.1e}, ice sheets on {:.3} of the planet",
        run.energy, run.water, run.most_sheet
    );
    assert!(
        run.energy.abs() < 0.1,
        "energy residual {} W/m²",
        run.energy
    );
    assert!(run.water.abs() < 1e-11, "water residual {}", run.water);
    assert!(run.most_sheet > 0.1, "ice sheets {}", run.most_sheet);
}

#[test]
fn energy_and_water_close_with_ekman_transport() {
    // The wind alone moves the ocean, against a continent's coasts.
    let params = PhysicsParams {
        ocean: OceanParams {
            ekman: Some(EkmanParams::default()),
            eddies: None,
            ..OceanParams::default()
        },
        ..PhysicsParams::default()
    };
    let run = budget_run(&Planet::earth(), params, 45.0, continent_fraction);
    eprintln!(
        "energy residual {:+.5} W/m², water residual {:+.1e}, upwelling up to {:.1e} kg/m²/s",
        run.energy, run.water, run.most_upwelling
    );
    assert!(
        run.energy.abs() < 0.1,
        "energy residual {} W/m²",
        run.energy
    );
    assert!(run.water.abs() < 1e-11, "water residual {}", run.water);
    assert!(run.most_upwelling > 0.0, "the water never moved");
}

#[test]
fn energy_and_water_close_with_eddy_transport() {
    // Eddies and the wind move the ocean together, as by default.
    let params = PhysicsParams::default();
    assert_eq!(
        params.ocean.eddies,
        Some(EddyParams::default()),
        "the default ocean has eddies"
    );
    let run = budget_run(&Planet::earth(), params, 45.0, continent_fraction);
    eprintln!(
        "energy residual {:+.5} W/m², water residual {:+.1e}, upwelling up to {:.1e} kg/m²/s",
        run.energy, run.water, run.most_upwelling
    );
    assert!(
        run.energy.abs() < 0.1,
        "energy residual {} W/m²",
        run.energy
    );
    assert!(run.water.abs() < 1e-11, "water residual {}", run.water);
}

#[test]
fn the_oceans_earth_turns_as_the_planets_does() {
    // So that the eddy coefficients are exactly their Earth values on Earth.
    assert_eq!(Planet::earth().rotation_rate(), EARTH_ROTATION);
}

/// The budgets of a run on [`CLIMATE_LEVELS`] at 642 cells.
struct BudgetRun {
    /// Energy residual (W/m²).
    energy: f64,
    /// Water residual, as a share of the atmosphere's water.
    water: f64,
    /// The sea's largest change in one step (K), and its day.
    largest_change: f64,
    largest_day: f64,
    /// The least and most of the planet under ice (after the first day).
    least_ice: f64,
    most_ice: f64,
    /// The most of the planet that was ever ice sheet.
    most_sheet: f64,
    /// The fastest upwelling anywhere (kg m⁻² s⁻¹).
    most_upwelling: f64,
}

fn budget_run(
    planet: &Planet,
    params: PhysicsParams,
    days: f64,
    land: impl Fn(glam::DVec3) -> f64,
) -> BudgetRun {
    let mut first = None;
    let mut last = None;
    let (mut previous, mut largest, mut when): (Option<Vec<f64>>, f64, f64) = (None, 0.0, 0.0);
    let (mut least_ice, mut most_ice) = (f64::INFINITY, 0.0f64);
    let mut most_sheet = 0.0f64;
    let mut most_upwelling = 0.0f64;
    let levels = Levels::new(&CLIMATE_LEVELS);
    run_on(
        planet,
        params,
        8,
        days,
        levels,
        land,
        FLAT,
        |si, geometry, dt| {
            let forcing = si.model().forcing().unwrap();
            let row = (
                si.model().diagnostics(si.state()),
                si.energy_budget(),
                forcing.budget().unwrap(),
                si.steps() as f64 * dt,
                si.water(),
                si.water_budget(),
            );
            first.get_or_insert(row);
            last = Some(row);
            let fields = forcing.cell_fields();
            let field = |name: &str| &fields.iter().find(|f| f.name == name).unwrap().values;
            let sea = field("ocean.mixed_temperature");
            if let Some(before) = &previous {
                let change = sea
                    .iter()
                    .zip(before)
                    .map(|(a, b): (&f64, &f64)| (a - b).abs())
                    .fold(0.0, f64::max);
                if change > largest {
                    (largest, when) = (change, row.3 / DAY);
                }
            }
            previous = Some(sea.clone());
            let sheet = field("land.ice_sheet")
                .iter()
                .zip(&geometry.cell_areas)
                .map(|(f, a)| f * a)
                .sum::<f64>()
                / (4.0 * std::f64::consts::PI);
            most_sheet = most_sheet.max(sheet);
            most_upwelling = field("ocean.upwelling")
                .iter()
                .fold(most_upwelling, |a, &w| a.max(w));
            if row.3 > DAY {
                let ice = field("ocean.ice_fraction")
                    .iter()
                    .zip(&geometry.cell_areas)
                    .map(|(f, a)| f * a)
                    .sum::<f64>()
                    / (4.0 * std::f64::consts::PI);
                least_ice = least_ice.min(ice);
                most_ice = most_ice.max(ice);
            }
        },
    );
    let ((d0, e0, b0, t0, w0, wb0), (d1, e1, b1, t1, w1, wb1)) = (first.unwrap(), last.unwrap());
    let area = 4.0 * std::f64::consts::PI * RADIUS * RADIUS;
    let toa = (b1.absorbed - b0.absorbed) - (b1.outgoing - b0.outgoing);
    let uptake = (d1.energy - d0.energy) + (b1.stored - b0.stored) + LATENT_HEAT * (w1 - w0);
    let removed = (e1.forcing.kinetic - e0.forcing.kinetic)
        + (b1.dissipated - b0.dissipated)
        + (e1.diffusion.total() - e0.diffusion.total());
    let water = (w1 - w0) - (b1.ocean_exchange - b0.ocean_exchange)
        + (b1.soil_water - b0.soil_water)
        + (b1.snow - b0.snow)
        + (b1.runoff - b0.runoff)
        - (wb1.filter - wb0.filter)
        - (wb1.fixer - wb0.fixer);
    BudgetRun {
        energy: (toa - uptake + removed) / (t1 - t0) / area,
        water: water / w1,
        largest_change: largest,
        largest_day: when,
        least_ice,
        most_ice,
        most_sheet,
        most_upwelling,
    }
}

#[test]
fn the_latitude_start_keeps_the_mean_and_warms_the_equator() {
    let planet = Planet::earth();
    let (mesh, _) = SphereMesh::centroidal(8, LloydOptions::default());
    let geometry = MeshGeometry::new(&mesh);
    let model = LayeredModel::new(
        Operators::new(&mesh, &geometry, RADIUS),
        GRAVITY,
        Gas::EARTH_AIR,
        Levels::new(&CLIMATE_LEVELS),
        vec![0.0; mesh.n_vertices()],
        vec![0.0; mesh.n_cells()],
    );
    let params = PhysicsParams::default();
    let start = initial_surface_temperature_by_latitude(&model, &planet, &params, 288.0);
    let total: f64 = geometry.cell_areas.iter().sum();
    let mean = start
        .iter()
        .zip(&geometry.cell_areas)
        .map(|(t, a)| t * a)
        .sum::<f64>()
        / total;
    assert!((mean - 288.0).abs() < 1e-9, "mean {mean}");
    let at = |z: f64| {
        let centers = mesh.cell_centers();
        let i = (0..centers.len())
            .min_by(|&a, &b| {
                (centers[a].z - z)
                    .abs()
                    .total_cmp(&(centers[b].z - z).abs())
            })
            .unwrap();
        start[i]
    };
    assert!((295.0..301.0).contains(&at(0.0)), "equator {}", at(0.0));
    assert!((262.0..275.0).contains(&at(1.0)), "pole {}", at(1.0));
}

#[test]
fn a_year_on_land_swings_more_than_ocean_and_snow_comes_and_goes() {
    // A year after two months of spin-up, in the mid-latitudes of the
    // continent and of the ocean around it. Melting happens here, so the
    // budgets are checked across it.
    let (mesh, _) = SphereMesh::centroidal(8, LloydOptions::default());
    let band = |p: &glam::DVec3| (30.0..60.0).contains(&fms_math::asin(p.z).abs().to_degrees());
    let cells = |tile: f64| -> Vec<usize> {
        (0..mesh.n_cells())
            .filter(|&i| {
                let p = mesh.cell_centers()[i];
                continent_fraction(p) == tile && band(&p)
            })
            .collect()
    };
    let (land, ocean) = (cells(1.0), cells(0.0));
    assert!(land.len() >= 8 && ocean.len() >= 8);
    let spin_up = 60.0;
    let days = spin_up + 365.0;
    let mut record: Vec<Vec<f64>> = vec![Vec::new(); mesh.n_cells()];
    let mut per_day = 0;
    let (mut driest, mut wettest) = (f64::MAX, f64::MIN);
    let (mut most_snow, mut least_snow) = (vec![f64::MIN; land.len()], vec![f64::MAX; land.len()]);
    budgets_close(8, days, continent_fraction, FLAT, |si, _, dt| {
        per_day = (DAY / dt).round() as usize;
        let fields = si.model().forcing().unwrap().cell_fields();
        let water = fields.iter().find(|f| f.name == "land.soil_water").unwrap();
        for &i in &land {
            driest = driest.min(water.values[i]);
            wettest = wettest.max(water.values[i]);
        }
        if si.steps() as f64 * dt <= spin_up * DAY {
            return;
        }
        let snow = fields.iter().find(|f| f.name == "land.snow").unwrap();
        for (j, &i) in land.iter().enumerate() {
            most_snow[j] = most_snow[j].max(snow.values[i]);
            least_snow[j] = least_snow[j].min(snow.values[i]);
        }
        let ts = fields
            .iter()
            .find(|f| f.name == "rad.surface_temperature")
            .unwrap();
        for &i in land.iter().chain(&ocean) {
            record[i].push(ts.values[i]);
        }
    });
    // The bucket stays within its limits and really does limit evaporation
    // somewhere: it starts at the critical level, 112.5 kg/m².
    let capacity = PhysicsParams::default().land.soil_water_capacity;
    eprintln!("soil water over mid-latitude land: {driest:.1} to {wettest:.1} kg/m²");
    assert!(driest >= 0.0 && wettest <= capacity + 1e-9);
    assert!(driest < 0.5 * capacity, "the soil never dried: {driest}");
    let average = |x: &[f64]| x.iter().sum::<f64>() / x.len() as f64;
    let (most, least) = (average(&most_snow), average(&least_snow));
    eprintln!(
        "snow over mid-latitude land: {least:.2} at the year's least, {most:.1} kg/m² at its most"
    );
    assert!(most > 5.0, "no winter snow: {most}");
    assert!(least < 0.5, "the snow never melted: {least}");
    let range = |x: &[f64]| {
        let max = x.iter().copied().fold(f64::MIN, f64::max);
        let min = x.iter().copied().fold(f64::MAX, f64::min);
        max - min
    };
    let mean = |x: &[f64]| x.iter().sum::<f64>() / x.len() as f64;
    let daily = |i: usize| {
        mean(
            &record[i]
                .chunks_exact(per_day)
                .map(range)
                .collect::<Vec<_>>(),
        )
    };
    let seasonal = |i: usize| {
        range(
            &record[i]
                .chunks_exact(10 * per_day)
                .map(mean)
                .collect::<Vec<_>>(),
        )
    };
    let over = |cells: &[usize], f: &dyn Fn(usize) -> f64| {
        cells.iter().map(|&i| f(i)).sum::<f64>() / cells.len() as f64
    };
    let (land_day, ocean_day) = (over(&land, &daily), over(&ocean, &daily));
    let (land_year, ocean_year) = (over(&land, &seasonal), over(&ocean, &seasonal));
    eprintln!(
        "daily range: land {land_day:.2} K, ocean {ocean_day:.2} K; \
         seasonal range: land {land_year:.1} K, ocean {ocean_year:.1} K"
    );
    assert!(
        land_day > 3.0 * ocean_day,
        "daily: land {land_day}, ocean {ocean_day}"
    );
    // Measured over the 50 m ocean: seasonal land 44.4 K, ocean 8.6 K, and
    // daily 10.3 K and 0.06 K.
    assert!(
        land_year > 3.0 * ocean_year,
        "seasonal: land {land_year}, ocean {ocean_year}"
    );
}

/// The sub-grid height spread and mean height (m) of a north–south ridge
/// along the prime meridian.
fn ridge(p: glam::DVec3) -> (f64, f64) {
    let lon = fms_math::atan2(p.y, p.x).to_degrees() / 8.0;
    let shape = fms_math::exp(-0.5 * lon * lon);
    (800.0 * shape, 2500.0 * shape)
}

#[test]
fn wind_up_a_slope_rains_on_the_windward_side_and_not_the_lee() {
    // A steady 10 m/s westerly in the boundary layer, on an all-land world
    // with the ridge's slope, and one evaluation of the physics.
    const WIND: f64 = 10.0;
    let planet = Planet::earth();
    let (mesh, _) = SphereMesh::centroidal(16, LloydOptions::default());
    let geometry = MeshGeometry::new(&mesh);
    let params = PhysicsParams::default();
    let centers = mesh.cell_centers();
    let model = LayeredModel::new(
        Operators::new(&mesh, &geometry, RADIUS),
        GRAVITY,
        Gas::EARTH_AIR,
        Levels::new(&CLIMATE_LEVELS),
        mesh.vertices()
            .iter()
            .map(|p| 2.0 * planet.rotation_rate() * p.z)
            .collect(),
        vec![0.0; mesh.n_cells()],
    );
    let (mut state, ts) = initial_state(&model, &planet, &params);
    for (u, edge) in state.u[CLIMATE_LEVELS.len()].iter_mut().zip(mesh.edges()) {
        let [a, b] = edge.cells.map(|c| centers[c as usize]);
        let east = glam::DVec3::Z.cross(a + b).normalize();
        *u = WIND * east.dot((b - a).normalize());
    }
    let mut physics = ColumnPhysics::with_land(
        &model,
        &planet,
        params,
        0.0,
        &state,
        vec![ts; mesh.n_cells()],
        vec![1.0; mesh.n_cells()],
    )
    .with_relief(
        centers.iter().map(|&p| ridge(p).0).collect(),
        &centers.iter().map(|&p| ridge(p).1).collect::<Vec<_>>(),
    );
    // The fluxes hold the ascent, so they need the relief first.
    physics.advance(&state, 1.0);
    let mut tendency = state.clone();
    physics.add_tendency(&state, &mut tendency);
    let fields = physics.cell_fields();
    let rain = &fields
        .iter()
        .find(|f| f.name == "atm.orographic_precipitation")
        .unwrap()
        .values;
    let flank = |west: f64, east: f64| {
        let cells: Vec<usize> = (0..mesh.n_cells())
            .filter(|&i| {
                let lon = fms_math::atan2(centers[i].y, centers[i].x).to_degrees();
                let lat = fms_math::asin(centers[i].z).abs().to_degrees();
                (west..east).contains(&lon) && lat < 60.0
            })
            .collect();
        cells.iter().map(|&i| rain[i]).sum::<f64>() / cells.len() as f64
    };
    let (windward, lee) = (flank(-16.0, -4.0), flank(4.0, 16.0));
    eprintln!("windward {windward:.2} mm/day, lee {lee:.4}");
    assert!(windward > 0.5, "windward rain {windward} mm/day");
    assert!(lee < 0.01 * windward, "lee rain {lee} mm/day");
    assert!(rain.iter().all(|r| *r >= 0.0));
}

#[test]
fn hidden_hills_slow_the_wind_and_return_its_energy_as_heat() {
    // Hills the atmosphere can't see, over land: nothing but their drag acts
    // on the wind. The budgets are checked across it, which closes only if
    // the wind's lost energy comes back as heat.
    let (mesh, _) = SphereMesh::centroidal(16, LloydOptions::default());
    let crest: Vec<usize> = (0..mesh.n_cells())
        .filter(|&i| ridge(mesh.cell_centers()[i]).0 > 400.0)
        .collect();
    let days = 60.0;
    let mesh = &mesh;
    let bottom = CLIMATE_LEVELS.len();
    let hills = |p: glam::DVec3| (ridge(p).0, 0.0);
    let measure = |run: &mut dyn FnMut(EachStep)| {
        let (mut speed, mut samples, mut heated) = (0.0, 0.0, 0.0);
        run(&mut |si, geometry, dt| {
            heated = si.model().forcing().unwrap().budget().unwrap().dissipated;
            if si.steps() as f64 * dt > 15.0 * DAY {
                let wind = geometry.reconstruct_vectors(mesh, &si.state().u[bottom]);
                speed += crest.iter().map(|&i| wind[i].length()).sum::<f64>() / crest.len() as f64;
                samples += 1.0;
            }
        });
        (speed / samples, heated)
    };
    let (speed, heated) =
        measure(&mut |each| budgets_close(16, days, continent_fraction, hills, each));
    let (flat_speed, flat_heated) = measure(&mut |each| {
        let levels = Levels::new(&CLIMATE_LEVELS);
        run_with(16, days, levels, continent_fraction, FLAT, each);
    });
    eprintln!(
        "wind on the crest {speed:.2} m/s, without the hills {flat_speed:.2}; \
         heat from friction {heated:.3e} J, without {flat_heated:.3e}"
    );
    assert!(
        speed < 0.9 * flat_speed,
        "wind {speed}, without {flat_speed}"
    );
    assert!(heated > flat_heated, "heat {heated}, without {flat_heated}");
}

#[test]
fn more_layers_run_the_moist_aquaplanet() {
    // Three layers (a thin bottom one) and five: 60 days from the single
    // column's equilibrium without blowing up, with the water budget exact,
    // no water conjured by the negative-humidity fixer, and the energy budget
    // closed over days 30–60 (measured: −0.007 and −0.003 W/m²; two layers
    // +0.023). The climate is calibrated for two layers only (O-25); step 10
    // recalibrates on three.
    for levels in [Levels::new(&[0.4, 0.85]), Levels::equal_mass(5)] {
        let n = levels.count();
        let mut first = None;
        let mut last = None;
        let mut fields: Vec<(String, f64)> = Vec::new();
        run_with(
            16,
            60.0,
            levels,
            |_| 0.0,
            FLAT,
            |si, geometry, dt| {
                let forcing = si.model().forcing().unwrap();
                let row = (
                    si.model().diagnostics(si.state()),
                    si.energy_budget(),
                    forcing.budget().unwrap(),
                    si.steps() as f64 * dt,
                    si.water(),
                    si.water_budget(),
                );
                if first.is_none() || (si.steps() as f64 * dt - 30.0 * 86_400.0).abs() < 0.5 * dt {
                    first = Some(row);
                }
                last = Some(row);
                let total: f64 = geometry.cell_areas.iter().sum();
                fields = forcing
                    .cell_fields()
                    .into_iter()
                    .map(|f| {
                        let mean = f
                            .values
                            .iter()
                            .zip(&geometry.cell_areas)
                            .map(|(v, a)| v * a)
                            .sum::<f64>()
                            / total;
                        (f.name.to_string(), mean)
                    })
                    .collect();
            },
        );
        let ((d0, e0, b0, t0, w0, wb0), (d1, e1, b1, t1, w1, wb1)) =
            (first.unwrap(), last.unwrap());
        let area = 4.0 * std::f64::consts::PI * RADIUS * RADIUS;
        let wm2 = |j: f64| j / (t1 - t0) / area;
        let toa = (b1.absorbed - b0.absorbed) - (b1.outgoing - b0.outgoing);
        let uptake = (d1.energy - d0.energy) + (b1.stored - b0.stored) + LATENT_HEAT * (w1 - w0);
        let removed = (e1.forcing.kinetic - e0.forcing.kinetic)
            + (b1.dissipated - b0.dissipated)
            + (e1.diffusion.total() - e0.diffusion.total());
        let residual = wm2(toa - uptake + removed);
        let exchanged = (b1.evaporated - b0.evaporated) - (b1.precipitated - b0.precipitated);
        let water_residual =
            (w1 - w0) - exchanged - (wb1.filter - wb0.filter) - (wb1.fixer - wb0.fixer);
        let get = |name: &str| fields.iter().find(|f| f.0 == name).unwrap().1;
        eprintln!(
            "{n} layers, day 60: T_s {:.1} K, OLR {:.1}, ASR {:.1} W/m², cloud {:.2}, water {:.1} kg/m², \
             precipitation {:.2} mm/day; days 30–60: energy residual {residual:+.4} W/m², water residual {:+.1e}",
            get("rad.surface_temperature"),
            get("rad.outgoing_longwave"),
            get("rad.absorbed_shortwave"),
            get("atm.cloud_cover"),
            get("atm.water_vapour"),
            get("atm.precipitation"),
            water_residual / w1
        );
        assert!(
            water_residual.abs() < 1e-12 * w1,
            "{n} layers: water residual {water_residual:e} kg"
        );
        assert!(
            (wb1.fixer - wb0.fixer).abs() < 1e-9 * w1,
            "{n} layers: fixer added water"
        );
        assert!(
            residual.abs() < 0.05,
            "{n} layers: energy residual {residual} W/m²"
        );
        assert!((d1.mass / d0.mass - 1.0).abs() < 1e-13);
    }
}

#[test]
#[ignore = "slow: five model years (~6 minutes); run with --ignored"]
fn the_equilibrium_climate_is_an_earthlike_aquaplanet() {
    let planet = Planet::earth();
    let year = planet.year_seconds();
    let years = 5.0;
    let mut start = None;
    let mut sums: Vec<f64> = Vec::new();
    let mut names: Vec<&str> = Vec::new();
    let mut n = 0.0;
    let mut end = None;
    run(16, years * year / DAY, |si, geometry, dt| {
        let t = si.steps() as f64 * dt;
        if t < (years - 1.0) * year {
            return;
        }
        let forcing = si.model().forcing().unwrap();
        let row = (
            forcing.budget().unwrap(),
            si.model().diagnostics(si.state()),
            si.energy_budget(),
            t,
            si.water(),
        );
        if start.is_none() {
            start = Some(row);
        }
        if si.steps() % 8 == 0 {
            let fields = forcing.cell_fields();
            let total: f64 = geometry.cell_areas.iter().sum();
            if sums.is_empty() {
                sums = vec![0.0; fields.len()];
                names = fields.iter().map(|f| f.name).collect();
            }
            for (s, f) in sums.iter_mut().zip(&fields) {
                *s += f
                    .values
                    .iter()
                    .zip(&geometry.cell_areas)
                    .map(|(v, a)| v * a)
                    .sum::<f64>()
                    / total;
            }
            n += 1.0;
        }
        end = Some(row);
    });
    let ((b0, d0, e0, t0, w0), (b1, d1, e1, t1, w1)) = (start.unwrap(), end.unwrap());
    let area = 4.0 * std::f64::consts::PI * RADIUS * RADIUS;
    let wm2 = |j: f64| j / (t1 - t0) / area;
    let (asr, olr) = (
        wm2(b1.absorbed - b0.absorbed),
        wm2(b1.outgoing - b0.outgoing),
    );
    let albedo = 1.0 - asr / (planet.mean_flux() / 4.0);
    let uptake = wm2((d1.energy - d0.energy) + (b1.stored - b0.stored) + LATENT_HEAT * (w1 - w0));
    let removed = -wm2((e1.forcing.kinetic - e0.forcing.kinetic)
        + (b1.dissipated - b0.dissipated)
        + (e1.diffusion.total() - e0.diffusion.total()));
    let mean = |name: &str| sums[names.iter().position(|n| *n == name).unwrap()] / n;
    let ts = mean("rad.surface_temperature");
    let (precipitation, cloud) = (mean("atm.precipitation"), mean("atm.cloud_cover"));
    let rh = [
        mean("atm.relative_humidity_upper"),
        mean("atm.relative_humidity_lower"),
    ];
    eprintln!(
        "final year: T_s {ts:.2} K, ASR {asr:.1}, OLR {olr:.1}, albedo {albedo:.3}, cloud {cloud:.2}, \
         RH {:.2}/{:.2}, P {precipitation:.2} mm/day, water {:.1} kg/m²; \
         ASR − OLR {:+.2} W/m² = uptake {uptake:+.2} + kinetic energy removed {removed:.2}",
        rh[0],
        rh[1],
        mean("atm.water_vapour"),
        asr - olr
    );
    // An ice-free aquaplanet with Earth's clouds runs warmer than Earth
    // (288 K), whose ice, snow and land brighten the clear sky; measured
    // 290.7 K (O-23). Earth: 240 W/m² each way, albedo 0.29–0.30, cloud
    // 0.6–0.7. On four layers "upper" is the nearly dry layer above 250 hPa
    // and "lower" the boundary layer.
    // Measured: OLR 230.7, albedo 0.322, cloud 0.71, RH 0.12/0.83, P 3.8 mm/day.
    assert!((288.0..294.0).contains(&ts), "T_s {ts}");
    assert!((225.0..250.0).contains(&olr), "OLR {olr}");
    assert!((0.27..0.33).contains(&albedo), "albedo {albedo}");
    assert!((0.55..0.75).contains(&cloud), "cloud {cloud}");
    assert!(
        (rh[0] - 0.12).abs() < 0.04 && (rh[1] - 0.83).abs() < 0.04,
        "RH {rh:?}"
    );
    // Too wet (O-24), but bounded.
    assert!((3.0..6.0).contains(&precipitation), "P {precipitation}");
    // Nearly settled: little heat still going into the atmosphere and slab
    // (measured 0.5 W/m²; the moist aquaplanet settles more slowly than the
    // fixed-humidity one did). The top of the atmosphere doesn't balance
    // exactly: friction and diffusion remove kinetic energy without
    // returning it as heat (O-21). The budget closes to the leapfrog's
    // time-stepping error, ~0.1 W/m² with moist physics (measured 0.12).
    assert!(uptake.abs() < 1.0, "still drifting: {uptake} W/m²");
    assert!((asr - olr - uptake - removed).abs() < 0.2);
}

#[test]
#[ignore = "slow: 180 days at fixed sea temperatures (~1.5 minutes); run with --ignored"]
fn the_fixed_sst_aquaplanet_has_trades_and_a_rain_belt() {
    // The Aqua-Planet Experiment's CONTROL sea temperatures on a perpetual
    // equinox, days 60–180 (the atmosphere settles in a month). Measured:
    // rain 2.5–2.9 mm/day at 25° and 5.9–6.1 at 5° (Earth's subtropics rain
    // ~2: O-22, a limit of the circulation rather than the boundary layer),
    // bottom-layer easterlies of about −6 m/s at 15° and westerlies of about
    // +6 at 35°.
    let planet = Planet::new(PlanetParams {
        obliquity_deg: 0.0,
        eccentricity: 0.0,
        ..PlanetParams::earth()
    })
    .0;
    let (mesh, _) = SphereMesh::centroidal(16, LloydOptions::default());
    let geometry = MeshGeometry::new(&mesh);
    let params = PhysicsParams::default();
    let model = LayeredModel::new(
        Operators::new(&mesh, &geometry, RADIUS),
        GRAVITY,
        Gas::EARTH_AIR,
        Levels::new(&CLIMATE_LEVELS),
        mesh.vertices()
            .iter()
            .map(|p| 2.0 * planet.rotation_rate() * p.z)
            .collect(),
        vec![0.0; mesh.n_cells()],
    )
    .with_humidity_profile(params.radiation);
    let (state, _) = initial_state(&model, &planet, &params);
    let dt = LayeredSemiImplicit::stable_dt(&model, 120.0, 0.6);
    let dt = DAY / (DAY / dt).ceil();
    let diffusion = HyperdiffusionOptions::strongest_for_step(model.operators(), 4, dt);
    let sst = mesh
        .cell_centers()
        .iter()
        .map(|&p| ApeSst::Control.at(fms_grid::coords::to_lat_lon(p).0))
        .collect();
    let physics =
        ColumnPhysics::new(&model, &planet, params, 0.0, &state, sst).with_prescribed_surface();
    let model = model
        .with_hyperdiffusion(diffusion)
        .with_dissipative_heating()
        .with_forcing(Box::new(physics));
    let mut si =
        LayeredSemiImplicit::new(model, state, dt, LayeredSemiImplicitOptions::default(), 0.5);
    const BANDS: usize = 18;
    let band = |lat: f64| ((lat + 90.0) / 10.0) as usize;
    let bands = LatitudeBands::new(&mesh, &geometry, BANDS);
    let (mut rain, mut wind, mut samples) = (vec![0.0; BANDS], vec![0.0; BANDS], 0.0);
    let bottom = si.model().levels().count() - 1;
    while (si.steps() as f64) * dt < 180.0 * DAY {
        si.step();
        si.check().expect("the fixed-SST aquaplanet blew up");
        if (si.steps() as f64) * dt > 60.0 * DAY && si.steps() % 8 == 0 {
            let fields = si.model().forcing().unwrap().cell_fields();
            let p = fields
                .iter()
                .find(|f| f.name == "atm.precipitation")
                .unwrap();
            for (s, x) in rain.iter_mut().zip(bands.mean(&geometry, &p.values)) {
                *s += x;
            }
            let v = geometry.reconstruct_vectors(&mesh, &si.state().u[bottom]);
            let east: Vec<f64> = mesh
                .cell_centers()
                .iter()
                .zip(&v)
                .map(|(&p, &v)| fms_grid::zonal::east_north(p, v).0)
                .collect();
            for (s, x) in wind.iter_mut().zip(bands.mean(&geometry, &east)) {
                *s += x;
            }
            samples += 1.0;
        }
    }
    let rain_at = |lat: f64| rain[band(lat)] / samples;
    let wind_at = |lat: f64| wind[band(lat)] / samples;
    eprintln!(
        "rain 5°/25°: {:.2}/{:.2} N, {:.2}/{:.2} S mm/day; bottom-layer wind 15°/35°: {:+.1}/{:+.1} N, {:+.1}/{:+.1} S m/s",
        rain_at(5.0),
        rain_at(25.0),
        rain_at(-5.0),
        rain_at(-25.0),
        wind_at(15.0),
        wind_at(35.0),
        wind_at(-15.0),
        wind_at(-35.0),
    );
    for hemisphere in [1.0, -1.0] {
        let (equator, subtropics) = (rain_at(hemisphere * 5.0), rain_at(hemisphere * 25.0));
        assert!(subtropics < 3.3, "subtropics rain {subtropics} mm/day");
        assert!(equator > 1.5 * subtropics, "no rain belt: {equator}");
        let (trades, westerlies) = (wind_at(hemisphere * 15.0), wind_at(hemisphere * 35.0));
        assert!(trades < -2.0, "trades {trades} m/s");
        assert!(westerlies > 2.0, "westerlies {westerlies} m/s");
    }
}
