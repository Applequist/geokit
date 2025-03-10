//! This module defines [Coordinates Reference Systems](Crs).
//! A [Crs] ties a coordinates system to the Earth using a [GeodeticDatum](datum) and allows to
//! unambiguously assign coordinates to location on Earth.

use std::any::Any;

use approx::AbsDiffEq;

use crate::{math::fp::Float, quantities::length::Length};

/// [Crs] is a marker trait that must be implemented by Coordinates reference systems.
/// It guarantees that CRS implement Any.
pub trait Crs: Any {
    /// Returns the unique id of the CRS.
    fn id(&self) -> &str;

    /// Returns the dimension of the CRS, either 1, 2 or 3.
    fn dim(&self) -> usize;

    /// Returns whether two positions `a` and `b` are all equal within the given tolerance `tol`.
    ///
    /// `a` and `b` are approximately equal if the absolute difference of their coordinates on
    /// each axis is less or equal to the tolerance for that axis.
    ///
    /// # Parameters
    ///
    /// - `a`: first position's coordinates. Length must match this CRS's dimension.
    /// - `b`: second position's coordinates. Length must match this CRS's dimension.
    /// - `tol`: the tolerance for each axis. Length, order and units must match this CRS's axes
    fn approx_eq(&self, a: &[Float], b: &[Float], tol: &[Float]) -> bool {
        let in2d = a[0].abs_diff_eq(&b[0], tol[0]) && a[1].abs_diff_eq(&b[1], tol[1]);
        match self.dim() {
            2 => in2d,
            3 => in2d && a[2].abs_diff_eq(&b[2], tol[2]),
            _ => unreachable!(),
        }
    }

    /// Returns the distance between the 2 given coordinates
    ///
    /// The computed distance depends on the type and dimension of the CRS.
    /// See each implementation for details.
    fn dist(&self, a: &[Float], b: &[Float]) -> Result<Length, &'static str>;
}

pub mod geocentric;
pub mod geographic;
pub mod projected;
