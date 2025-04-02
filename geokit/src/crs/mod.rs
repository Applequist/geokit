//! This module defines [Coordinates Reference Systems](Crs).
//! A [Crs] ties a coordinates system to the Earth using a [GeodeticDatum](datum) and allows to
//! unambiguously assign coordinates to location on Earth.

use crate::{cs::Coord, math::fp::Float, quantities::length::Length};
use approx::AbsDiffEq;
use std::any::Any;
use std::fmt::Debug;

/// [Crs] is a marker trait that must be implemented by Coordinates reference systems.
/// It guarantees that CRS implement Any.
pub trait Crs: Any + Debug {
    /// Returns the unique id of the CRS.
    fn id(&self) -> &str;

    /// Returns the dimension of the CRS, either 1, 2 or 3.
    fn dim(&self) -> usize;

    /// Returns whether two coordinates `a` and `b` are approximately equal within the given tolerance `tol`.
    ///
    /// `a` and `b` are approximately equal if the absolute difference of their coordinates on
    /// each axis is less or equal to the tolerance for that axis.
    ///
    /// # Parameters
    ///
    /// - `a`: first position's coordinates. Length must match this CRS's dimension.
    /// - `b`: second position's coordinates. Length must match this CRS's dimension.
    /// - `tol`: the tolerance for each axis. Length, order and units must match this CRS's axes
    fn approx_eq(&self, a: &Coord, b: &Coord, tol: &[Float]) -> bool {
        assert_eq!(a.len(), self.dim());
        assert_eq!(b.len(), self.dim());
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
    ///
    /// # Parameters
    ///
    /// - `a`: first position's coordinates. Length must match this CRS's dimension.
    /// - `b`: second position's coordinates. Length must match this CRS's dimension.
    fn dist(&self, a: &Coord, b: &Coord) -> Result<Length, &'static str>;
}

pub mod geocentric;
pub mod geographic;
pub mod projected;
