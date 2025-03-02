//! This module defines [Coordinates Reference Systems](Crs).
//! A [Crs] ties a coordinates system to the Earth using a [GeodeticDatum](datum) and allows to
//! unambiguously assign coordinates to location on Earth.

use std::any::Any;

/// [Crs] is a marker trait that must be implemented by Coordinates reference systems.
/// It guarantees that CRS implement Any.
pub trait Crs: Any {
    /// Returns the unique id of the CRS.
    fn id(&self) -> &str;
}

pub mod geocentric;
pub mod geographic;
pub mod projected;
