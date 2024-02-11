use std::any::Any;

pub trait Crs: Any {
    /// Returns whether the given Crs uses *normalized* coordinates,
    /// i.e. well-known axes order, direction and units.
    fn is_normalized(&self) -> bool;
}

pub mod geocentric;
pub mod geographic;
pub mod projected;

pub use geocentric::GeocentricCrs;
pub use geographic::{GeodeticAxes, GeographicCrs};
