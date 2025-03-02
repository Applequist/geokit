//! Cartesian coordinates systems are used to represent points from 2D and 3D cartesion space
//! using distances from hyperplanes along axes as coordinates.
//!
//! Cartesian coordinates systems are used for geocentric, projected and topocentric cartesian
//! spaces.
//!
//! This module Provide systems of axes for the following cartesian CS:
//! - [GeocentricAxes] for geocentric cartesian CS.
//! - [ProjectedAxes] for projected cartesian CS.

use crate::{
    math::fp::Float,
    quantities::length::Length,
    units::length::{LengthUnit, M},
};
use approx::AbsDiffEq;
use derive_more::derive::Display;

use super::CartesianErrors;

/// [GeocentricAxes] defines the possible set of axes used in **geocentric** 3D cartesion CS.
/// Geocentric CS uses [meter](crate::units::length::M) unit by default for all axes.
#[derive(Debug, Clone, Copy, PartialEq, Display)]
pub enum GeocentricAxes {
    /// Coordinates are:
    /// - x in meters. **The X-axis is the intersection of the equatirial
    /// plane and the Greenwich meridian plane, positive toward the Greenwich meridian.
    /// - y in meters, positive east
    /// - z in meters, positive north
    XYZ,
}

impl GeocentricAxes {
    pub fn normalize(&self, coords: &[Float]) -> XYZ {
        XYZ {
            x: Length::new(coords[0], M),
            y: Length::new(coords[1], M),
            z: Length::new(coords[2], M),
        }
    }

    pub fn denormalize(&self, xyz: XYZ, coords: &mut [Float]) {
        coords[0] = xyz.x.m();
        coords[1] = xyz.y.m();
        coords[2] = xyz.z.m();
    }

    pub fn dim(&self) -> usize {
        3
    }
}

/// [XYZ] represents **normalized** geocentric coordinates.
#[derive(Debug, Copy, Clone, PartialEq, Display)]
#[display("({}, {}, {})", x, y, z)]
pub struct XYZ {
    pub x: Length,
    pub y: Length,
    pub z: Length,
}

impl XYZ {
    /// Computes the distance between `self` and `other` in meters.
    pub fn dist_to(self, other: Self) -> Length {
        (self.x - other.x)
            .hypot(self.y - other.y)
            .hypot(self.z - other.z)
    }

    /// Checks whether this [XYZ] is approximately equal to the `other` [XYZ]
    /// within the given [CartesianErrors] bounds.
    ///
    /// `self` and `other` are approximately equal if the following conditions are
    /// all satisfied:
    /// - `self.dist_to(other) <= err.0`
    /// FIX: This is not the same as Crs::approx_eq.
    ///
    pub fn approx_eq(self, other: Self, err: CartesianErrors) -> bool {
        self.dist_to(other) <= err.length()
    }
}

/// Returns whether `res` is approximately equal to `exp` within the given [CartesianErrors] error
/// bounds, printing information about the coordinates when not equal.
///
/// Use only in tests.
pub fn approx_eq_xyz(res: XYZ, exp: XYZ, err: CartesianErrors) -> bool {
    let d = res.dist_to(exp);
    let is_approx_eq = d <= err.length();
    if !is_approx_eq {
        println!(
            "d({}, {}) = {:e} m > {:e} m",
            res,
            exp,
            d.m(),
            err.length().m()
        );
        if !res.x.abs_diff_eq(&exp.x, err.length()) {
            println!(
                "X error: | {} - {} | = {:e} m > {:e} m",
                res.x,
                exp.x,
                (res.x - exp.x).abs().m(),
                err.length().m()
            );
        } else {
            println!(
                "X ok      {} = {} +/- {:e} m",
                res.x,
                exp.x,
                err.length().m()
            );
        }
        if !res.y.abs_diff_eq(&exp.y, err.length()) {
            println!(
                "Y error: | {} - {} | = {:e} m > {:e} m",
                res.y,
                exp.y,
                (res.y - exp.y).abs().m(),
                err.length().m()
            );
        } else {
            println!(
                "Y ok      {} = {} +/- {:e} m",
                res.y,
                exp.y,
                err.length().m()
            );
        }
        if !res.z.abs_diff_eq(&exp.z, err.length()) {
            println!(
                "Z error: | {} - {} | = {:e} m > {:e} m",
                res.z,
                exp.z,
                (res.z - exp.z).abs().m(),
                err.length().m()
            );
        } else {
            println!(
                "Z ok      {} = {} +/- {:e} m",
                res.z,
                exp.z,
                err.length().m()
            );
        }
    }
    is_approx_eq
}
