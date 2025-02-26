//! The geokit library allows you to work with geographical data expressed
//! in various coordinate reference systems and transform them to and from other CRS.
//!
//! IMPORTANT: Wherever coordinates, whether geodetic or cartesian, appear as raw floating point values
//! they are **always** expressed in radians (geodetic) or meters (cartesian).

pub mod math;

pub mod cs;
pub mod geodesy;
pub mod operations;
pub mod projections;
pub mod quantities;
pub mod transformations;
pub mod units;
