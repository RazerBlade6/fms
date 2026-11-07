//! World generators: plate layout, noise, heightmap import.
//!
//! So far: importing a global heightmap ([`heightmap::load`]) and resampling it
//! onto the terrain grid ([`resample()`]), plus the sub-grid terrain statistics
//! the climate grid derives from it ([`TerrainStats`]).
//!
//! See `docs/PLAN.md` for the design.

pub mod heightmap;
pub mod raster;
pub mod resample;
pub mod subgrid;

pub use heightmap::{HeightmapError, LoadOptions, load};
pub use raster::{Raster, RasterError};
pub use resample::resample;
pub use subgrid::TerrainStats;
