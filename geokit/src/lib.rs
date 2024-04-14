//! The geokit library allows you to work with geographical coordinates expressed
//! in various coordinate reference systems and transform them to and from other CRS.
//! IMPORTANT: Wherever coordinates, whether geodetic or cartesian, appear as raw floating point values
//! they are **always** expressed in radians (geodetic) or meters (cartesian).
#[macro_use]
extern crate lazy_static;

extern crate nalgebra as na;

pub mod units;
pub mod cs;
pub mod geodesy;
pub mod crs;
pub mod operation;
pub mod providers;
