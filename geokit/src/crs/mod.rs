//! This module defines [Coordinates Reference Systems](Crs).
//! A [Crs] ties a coordinates system to the Earth using a [GeodeticDatum](datum) and allows to
//! unambiguously assign coordinates to location on Earth.

/// [Crs] is the root trait for Coordinates reference systems.
pub trait Crs {
    /// Return the unique id this [Crs].
    fn id(&self) -> &str;
}

pub mod geocentric;
pub mod geographic;
pub mod projected;
