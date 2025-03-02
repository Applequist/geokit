//! This module defines [Coordinates Reference Systems](Crs).
//! A [Crs] ties a coordinates system to the Earth using a [GeodeticDatum](datum) and allows to
//! unambiguously assign coordinates to location on Earth.

use std::any::Any;

use crate::{math::fp::Float, quantities::length::Length};

/// [Crs] is a marker trait that must be implemented by Coordinates reference systems.
/// It guarantees that CRS implement Any.
pub trait Crs: Any {
    type Tolerance;

    /// Returns the unique id of the CRS.
    fn id(&self) -> &str;

    /// Returns the dimension of the CRS, either 1, 2 or 3.
    fn dim(&self) -> usize;

    /// Returns whether the `a` and `b` coordinates are all equal within the given error margins.
    ///
    /// This method checks that each ordinates are approximately equal and does NOT rely
    /// on the notion of distance.
    fn approx_eq(&self, a: &[Float], b: &[Float], err: Self::Tolerance) -> bool;

    /// Returns the distance between the 2 given coordinates
    ///
    /// The computed distance depends on the type and dimension of the CRS.
    /// See each implementation for details.
    fn dist(&self, a: &[Float], b: &[Float]) -> Length;
}

pub mod geocentric;
pub mod geographic;
pub mod projected;
