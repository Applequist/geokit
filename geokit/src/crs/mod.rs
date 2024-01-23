use std::any::Any;

pub trait Crs: Any {}

pub mod geocentric;
pub mod geographic;

pub use geocentric::GeocentricCrs;
pub use geographic::{GeodeticAxes, GeographicCrs};
