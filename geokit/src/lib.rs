//! The geokit library allows you to work with geographical coordinates expressed
//! in various coordinate reference systems and transform them to and from other CRS.
//! IMPORTANT: Wherever coordinates, whether geodetic or cartesian, appear as raw floating point values
//! they are **always** expressed in radians (geodetic) or meters (cartesian).

#![warn(missing_docs)]

pub mod math;

pub mod crs;
pub mod cs;
pub mod geodesy;
pub mod operation;
pub mod providers;
pub mod quantity;
