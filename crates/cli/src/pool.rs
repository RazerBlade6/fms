//! The climate run's thread pool: how many threads it has, and whether the
//! idle ones spin or sleep between a step's parallel regions.
//!
//! A step is a sequence of short parallel regions. A worker with nothing to
//! do yields a few times and then sleeps, and waking it for the next region
//! costs tens of microseconds and a context switch. Where a step has work
//! for every thread, at grids whose loops are split between threads, that
//! waiting is the largest single cost: keeping the workers spinning, so they
//! take the next region's work at once, is 18% faster at 20,252 cells. On a
//! small grid the dynamics has only its layers to share out; spinning gains
//! nothing there, and threads beyond about one per 256 cells only add
//! wake-ups. Measured at M2.15 (docs/climate/BENCHMARKS.md).

use std::sync::Arc;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};

use clap::ValueEnum;
use fms_grid::operators::splits;

/// Cells to a thread on a grid too small for its loops to be split.
const CELLS_PER_THREAD: usize = 256;

/// What the pool's threads do when they have no work.
#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
pub enum Idle {
    /// Keep looking for work: every core stays busy for the whole run.
    Spin,
    /// Yield and then sleep until woken.
    Sleep,
}

/// A run's threads.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct PoolPlan {
    pub threads: usize,
    /// How many of them spin when idle (all but one, or none).
    pub spinners: usize,
}

/// The threads for a run on `cells` cells of a machine with `cores` cores.
/// Unless `threads` is asked for: every core where the grid's loops are
/// split between threads, otherwise one for every [`CELLS_PER_THREAD`]
/// cells. Unless `idle` is asked for: they spin only on such a grid with
/// every core in use (fewer threads leave the machine room, which spinning
/// would take back).
pub fn plan(cells: usize, cores: usize, threads: Option<usize>, idle: Option<Idle>) -> PoolPlan {
    let large = splits(cells);
    let cores = cores.max(1);
    let automatic = if large {
        cores
    } else {
        cores.min((cells / CELLS_PER_THREAD).max(1))
    };
    let threads = threads.unwrap_or(automatic).max(1);
    let spin = match idle {
        Some(Idle::Spin) => true,
        Some(Idle::Sleep) => false,
        None => large && threads >= cores,
    };
    PoolPlan {
        threads,
        // One thread is left to sleep: with the thread that runs the loop,
        // a spinner on every thread would be one more than the cores.
        spinners: if spin { threads - 1 } else { 0 },
    }
}

/// The plan for a run on `cells` cells, with the global pool sized to it.
/// `threads` falls back to `RAYON_NUM_THREADS`. A pool that already exists
/// keeps its size.
pub fn configure(cells: usize, threads: Option<usize>, idle: Option<Idle>) -> PoolPlan {
    let from_environment = || {
        std::env::var("RAYON_NUM_THREADS")
            .ok()
            .and_then(|n| n.parse().ok())
            .filter(|&n: &usize| n > 0)
    };
    let cores = std::thread::available_parallelism().map_or(1, usize::from);
    let wanted = plan(cells, cores, threads.or_else(from_environment), idle);
    let built = rayon::ThreadPoolBuilder::new()
        .num_threads(wanted.threads)
        .build_global();
    if built.is_ok() {
        return wanted;
    }
    let threads = rayon::current_num_threads();
    PoolPlan {
        threads,
        spinners: wanted.spinners.min(threads - 1),
    }
}

/// While this lives, some of the global pool's threads spin when they have
/// no work instead of sleeping.
pub struct Spinners {
    stop: Arc<AtomicBool>,
    count: usize,
}

impl Spinners {
    /// Start `count` spinners, on the pool's first threads. Starts none
    /// from inside the pool: the calling thread's own would sit on top of
    /// the caller and hold it until the end.
    pub fn start(count: usize) -> Self {
        let stop = Arc::new(AtomicBool::new(false));
        let count = count.min(rayon::current_num_threads());
        if count == 0 || rayon::current_thread_index().is_some() {
            return Self { stop, count: 0 };
        }
        let started = Arc::new(AtomicUsize::new(0));
        let (flag, running) = (Arc::clone(&stop), Arc::clone(&started));
        rayon::spawn_broadcast(move |context| {
            if context.index() >= count {
                return;
            }
            running.fetch_add(1, Ordering::Relaxed);
            while !flag.load(Ordering::Relaxed) {
                // Run a waiting piece of work if there is one.
                if rayon::yield_now() != Some(rayon::Yield::Executed) {
                    std::hint::spin_loop();
                }
            }
        });
        // Each spinner must sit at the bottom of its thread's stack before
        // the pool is given work: one started on top of a piece of work
        // would hold that work's thread until the end.
        while started.load(Ordering::Relaxed) < count {
            std::thread::yield_now();
        }
        Self { stop, count }
    }

    /// How many threads spin.
    pub fn count(&self) -> usize {
        self.count
    }
}

impl Drop for Spinners {
    fn drop(&mut self) {
        self.stop.store(true, Ordering::Relaxed);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rayon::prelude::*;

    #[test]
    fn a_large_grid_takes_every_core_and_spins_all_but_one() {
        let p = plan(20_252, 15, None, None);
        assert_eq!((p.threads, p.spinners), (15, 14));
        let p = plan(10_242, 8, None, None);
        assert_eq!((p.threads, p.spinners), (8, 7));
    }

    #[test]
    fn a_small_grid_takes_a_thread_for_every_256_cells_and_sleeps() {
        let p = plan(2_562, 15, None, None);
        assert_eq!((p.threads, p.spinners), (10, 0));
        assert_eq!(plan(2_562, 4, None, None).threads, 4);
        assert_eq!(plan(642, 15, None, None).threads, 2);
        assert_eq!(plan(42, 15, None, None).threads, 1);
    }

    #[test]
    fn asking_for_fewer_threads_stops_the_spinning_unless_asked_for() {
        let p = plan(20_252, 15, Some(4), None);
        assert_eq!((p.threads, p.spinners), (4, 0));
        let p = plan(20_252, 15, Some(4), Some(Idle::Spin));
        assert_eq!((p.threads, p.spinners), (4, 3));
        let p = plan(20_252, 15, None, Some(Idle::Sleep));
        assert_eq!((p.threads, p.spinners), (15, 0));
        let p = plan(2_562, 15, None, Some(Idle::Spin));
        assert_eq!((p.threads, p.spinners), (10, 9));
        // One thread has nothing to spin.
        assert_eq!(plan(2_562, 15, Some(1), Some(Idle::Spin)).spinners, 0);
    }

    #[test]
    fn the_pool_works_while_its_threads_spin_and_after() {
        let sum = |n: u64| (0..n).into_par_iter().sum::<u64>();
        let spinners = Spinners::start(2);
        assert!(spinners.count() <= 2);
        assert_eq!(sum(100_000), 4_999_950_000);
        drop(spinners);
        assert_eq!(sum(100_000), 4_999_950_000);
    }

    #[test]
    fn no_spinner_starts_from_inside_the_pool() {
        let count = rayon::scope(|_| Spinners::start(2).count());
        assert_eq!(count, 0);
    }
}
