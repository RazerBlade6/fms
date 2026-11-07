//! Simulation clock and calendar.
//!
//! Time is an integer count of base steps, so it never drifts and replays are
//! exact. Everything else (day, year, month, season) is derived from that count.
//!
//! A **day** is the planet's solar day (one day/night cycle), and the base step
//! divides it evenly. A **year** may be a non-integer number of days; year `y`
//! starts on day `⌊y · days_per_year⌋`, so years alternate between `⌊d⌋` and
//! `⌈d⌉` days the way leap years do. Year fraction 0 is the northern vernal
//! equinox by convention. Astronomical seasons come from the orbit, in
//! `fms-planet`.

use std::fmt;

/// Fixed time structure of a world: base step, day length and year length.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Calendar {
    day_seconds: f64,
    steps_per_day: u32,
    days_per_year: f64,
}

/// Invalid calendar parameters.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum CalendarError {
    /// Day length must be positive and finite.
    DayLength(f64),
    /// There must be at least one step per day.
    StepsPerDay,
    /// A year must be at least one day long and finite.
    YearLength(f64),
}

impl fmt::Display for CalendarError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            CalendarError::DayLength(d) => write!(f, "invalid solar day length {d} s"),
            CalendarError::StepsPerDay => write!(f, "steps per day must be at least 1"),
            CalendarError::YearLength(y) => write!(f, "invalid year length {y} days"),
        }
    }
}

impl std::error::Error for CalendarError {}

impl Calendar {
    /// Build a calendar.
    ///
    /// - `day_seconds`: length of the planet's solar day, in seconds.
    /// - `steps_per_day`: base steps per solar day, so `dt = day_seconds / steps_per_day`.
    /// - `days_per_year`: year length in solar days (may be fractional).
    pub fn new(
        day_seconds: f64,
        steps_per_day: u32,
        days_per_year: f64,
    ) -> Result<Self, CalendarError> {
        if !(day_seconds.is_finite() && day_seconds > 0.0) {
            return Err(CalendarError::DayLength(day_seconds));
        }
        if steps_per_day == 0 {
            return Err(CalendarError::StepsPerDay);
        }
        if !(days_per_year.is_finite() && days_per_year >= 1.0) {
            return Err(CalendarError::YearLength(days_per_year));
        }
        Ok(Self {
            day_seconds,
            steps_per_day,
            days_per_year,
        })
    }

    /// Earth: 86 400 s day, 365.2422-day year, 48 steps per day (30 min).
    pub fn earth() -> Self {
        Self {
            day_seconds: 86_400.0,
            steps_per_day: 48,
            days_per_year: 365.2422,
        }
    }

    /// Base step `dt` in seconds.
    pub fn dt(&self) -> f64 {
        self.day_seconds / f64::from(self.steps_per_day)
    }

    /// Solar day length in seconds.
    pub fn day_seconds(&self) -> f64 {
        self.day_seconds
    }

    /// Base steps per solar day.
    pub fn steps_per_day(&self) -> u32 {
        self.steps_per_day
    }

    /// Year length in solar days.
    pub fn days_per_year(&self) -> f64 {
        self.days_per_year
    }

    /// First day of year `year`.
    pub fn year_start_day(&self, year: u64) -> u64 {
        (year as f64 * self.days_per_year).floor() as u64
    }

    /// Date of a given step.
    pub fn date(&self, step: u64) -> Date {
        let steps_per_day = u64::from(self.steps_per_day);
        let day = step / steps_per_day;
        let step_of_day = (step % steps_per_day) as u32;

        // Estimate, then correct for floating-point edge cases at boundaries.
        let mut year = (day as f64 / self.days_per_year).floor() as u64;
        while year > 0 && self.year_start_day(year) > day {
            year -= 1;
        }
        while self.year_start_day(year + 1) <= day {
            year += 1;
        }
        let year_start = self.year_start_day(year);
        let days_in_year = (self.year_start_day(year + 1) - year_start) as u32;
        let day_of_year = (day - year_start) as u32;
        let month = (u64::from(day_of_year) * 12 / u64::from(days_in_year)) as u8;

        Date {
            step,
            day,
            step_of_day,
            year,
            day_of_year,
            days_in_year,
            month,
        }
    }
}

/// Position in time, derived from a step count by [`Calendar::date`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Date {
    /// Base steps since the start of the simulation.
    pub step: u64,
    /// Solar days since the start.
    pub day: u64,
    /// Step within the current day, `0..steps_per_day`.
    pub step_of_day: u32,
    /// Years since the start (year 0 is the first).
    pub year: u64,
    /// Day within the current year, `0..days_in_year`.
    pub day_of_year: u32,
    /// Length of the current year in days (`⌊d⌋` or `⌈d⌉`).
    pub days_in_year: u32,
    /// Month `0..12`: the year split into 12 equal parts by day.
    pub month: u8,
}

impl Date {
    /// Fraction of the current year elapsed, `[0, 1)`, including the time of day.
    pub fn year_fraction(&self, calendar: &Calendar) -> f64 {
        let day = f64::from(self.day_of_year)
            + f64::from(self.step_of_day) / f64::from(calendar.steps_per_day());
        day / f64::from(self.days_in_year)
    }

    /// Fraction of the current day elapsed, `[0, 1)`.
    pub fn day_fraction(&self, calendar: &Calendar) -> f64 {
        f64::from(self.step_of_day) / f64::from(calendar.steps_per_day())
    }

    /// Calendar quarter `0..4`, by year fraction. Quarter 0 starts at the
    /// northern vernal equinox by convention. This is a calendar division only:
    /// for real seasons, which differ in length on eccentric orbits, use
    /// `fms_planet::Planet::season`.
    pub fn quarter(&self, calendar: &Calendar) -> u8 {
        (self.year_fraction(calendar) * 4.0).floor() as u8 % 4
    }

    /// Whether this is the first step of a day.
    pub fn is_day_start(&self) -> bool {
        self.step_of_day == 0
    }

    /// Whether this is the first step of a month.
    pub fn is_month_start(&self, calendar: &Calendar) -> bool {
        if !self.is_day_start() {
            return false;
        }
        if self.day_of_year == 0 {
            return true;
        }
        let previous = calendar.date(self.step - u64::from(calendar.steps_per_day()));
        previous.month != self.month
    }

    /// Whether this is the first step of a year.
    pub fn is_year_start(&self) -> bool {
        self.is_day_start() && self.day_of_year == 0
    }
}
