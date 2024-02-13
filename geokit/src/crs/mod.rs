// FIXME: Do we really need this. Should be instead go for an enum with Geocentric, Geographic and
// Projected variants ?
pub trait Crs {
    /// Returns whether the given Crs uses *normalized* coordinates,
    /// i.e. well-known axes order, direction and units.
    fn is_normalized(&self) -> bool;
}

pub mod geocentric;
pub mod geographic;
pub mod projected;

pub use geocentric::GeocentricCrs;
pub use geographic::{GeodeticAxes, GeographicCrs};
