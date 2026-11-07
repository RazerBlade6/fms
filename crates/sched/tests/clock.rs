//! Calendar arithmetic: days, years, months, seasons.

use fms_sched::{Calendar, CalendarError, Hemisphere, Season};

#[test]
fn rejects_invalid_parameters() {
    assert_eq!(
        Calendar::new(0.0, 48, 365.0),
        Err(CalendarError::DayLength(0.0))
    );
    assert_eq!(
        Calendar::new(86_400.0, 0, 365.0),
        Err(CalendarError::StepsPerDay)
    );
    assert_eq!(
        Calendar::new(86_400.0, 48, 0.5),
        Err(CalendarError::YearLength(0.5))
    );
    assert!(Calendar::new(f64::NAN, 48, 365.0).is_err());
}

#[test]
fn dt_divides_the_day() {
    let c = Calendar::new(36_000.0, 24, 200.0).unwrap();
    assert_eq!(c.dt(), 1500.0);
    assert_eq!(Calendar::earth().dt(), 1800.0);
}

#[test]
fn days_and_steps() {
    let c = Calendar::new(86_400.0, 48, 365.0).unwrap();
    let d = c.date(48 * 10 + 7);
    assert_eq!(
        (d.day, d.step_of_day, d.year, d.day_of_year),
        (10, 7, 0, 10)
    );
    assert!(!d.is_day_start());
    assert!(c.date(48 * 10).is_day_start());
}

#[test]
fn fractional_years_alternate_like_leap_years() {
    let c = Calendar::new(1.0, 1, 365.25).unwrap();
    let lengths: Vec<u32> = (0..8)
        .map(|y| c.date(c.year_start_day(y)).days_in_year)
        .collect();
    assert_eq!(lengths, [365, 365, 365, 366, 365, 365, 365, 366]);
    // Four years add up to exactly 1461 days.
    assert_eq!(c.year_start_day(4), 1461);
}

#[test]
fn every_day_maps_to_a_consistent_date() {
    let c = Calendar::new(1.0, 1, 365.2422).unwrap();
    let mut previous = c.date(0);
    assert!(previous.is_year_start());
    for day in 1..(400 * 366) {
        let d = c.date(day);
        assert!(d.day_of_year < d.days_in_year);
        assert!(d.month < 12);
        if d.year == previous.year {
            assert_eq!(d.day_of_year, previous.day_of_year + 1);
            assert!(d.month >= previous.month);
        } else {
            assert_eq!(d.year, previous.year + 1);
            assert_eq!(previous.day_of_year + 1, previous.days_in_year);
            assert!(d.is_year_start());
        }
        previous = d;
    }
    // 400 years of 365.2422 days.
    assert_eq!(c.year_start_day(400), 146_096);
}

#[test]
fn twelve_months_per_year() {
    let c = Calendar::new(1.0, 1, 365.0).unwrap();
    let starts: Vec<u32> = (0..365)
        .filter(|&d| c.date(d).is_month_start(&c))
        .map(|d| d as u32)
        .collect();
    assert_eq!(starts.len(), 12);
    assert_eq!(starts[0], 0);
    // Short years still get 12 months.
    let tiny = Calendar::new(1.0, 1, 12.0).unwrap();
    assert_eq!(
        (0..12)
            .filter(|&d| tiny.date(d).is_month_start(&tiny))
            .count(),
        12
    );
}

#[test]
fn seasons_follow_the_year_and_flip_by_hemisphere() {
    let c = Calendar::new(1.0, 4, 100.0).unwrap();
    let at = |day: u64| c.date(day * 4);
    assert_eq!(at(0).season(&c, Hemisphere::North), Season::Spring);
    assert_eq!(at(30).season(&c, Hemisphere::North), Season::Summer);
    assert_eq!(at(60).season(&c, Hemisphere::North), Season::Autumn);
    assert_eq!(at(90).season(&c, Hemisphere::North), Season::Winter);
    assert_eq!(at(30).season(&c, Hemisphere::South), Season::Winter);
    assert_eq!(at(0).season(&c, Hemisphere::South), Season::Autumn);
    let d = c.date(4 * 50 + 2);
    assert!((d.year_fraction(&c) - 50.5 / 100.0).abs() < 1e-12);
    assert!((d.day_fraction(&c) - 0.5).abs() < 1e-12);
}
