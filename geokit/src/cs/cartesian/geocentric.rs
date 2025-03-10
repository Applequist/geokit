//! Cartesian coordinates systems are used to represent points from 2D and 3D cartesion space
//! using distances from hyperplanes along axes as coordinates.
//!
//! Cartesian coordinates systems are used for geocentric, projected and topocentric cartesian
//! spaces.
//!
//! This module Provide systems of axes for the following cartesian CS:
//! - [GeocentricAxes] for geocentric cartesian CS.
//! - [ProjectedAxes] for projected cartesian CS.

use super::CartesianTolerance;
use crate::{cs::Tolerance, math::fp::Float, quantities::length::Length, units::length::M};
use approx::AbsDiffEq;
use derive_more::derive::Display;
use smallvec::smallvec;

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
    pub fn dim(&self) -> usize {
        3
    }

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

    pub fn denormalize_tol(&self, tol: CartesianTolerance) -> Tolerance {
        smallvec![tol.all.m(); 3]
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
    pub fn approx_eq(self, other: Self, tol: CartesianTolerance) -> bool {
        self.x.abs_diff_eq(&other.x, tol.all)
            && self.y.abs_diff_eq(&other.y, tol.all)
            && self.z.abs_diff_eq(&other.z, tol.all)
    }
}

/// Returns whether `res` is approximately equal to `exp` within the given [CartesianErrors] error
/// bounds, printing information about the coordinates when not equal.
///
/// Use only in tests.
pub fn approx_eq_xyz(res: XYZ, exp: XYZ, tol: CartesianTolerance) -> bool {
    let is_approx_eq = res.approx_eq(exp, tol);
    if !is_approx_eq {
        println!("{} !~= {} (tol = {:e} m)", res, exp, tol.all.m());
        if !res.x.abs_diff_eq(&exp.x, tol.all) {
            println!(
                "X error: | {} - {} | = {:e} m > {:e} m",
                res.x,
                exp.x,
                (res.x - exp.x).abs().m(),
                tol.all.m()
            );
        } else {
            println!("X ok      {} = {} +/- {:e} m", res.x, exp.x, tol.all.m());
        }
        if !res.y.abs_diff_eq(&exp.y, tol.all) {
            println!(
                "Y error: | {} - {} | = {:e} m > {:e} m",
                res.y,
                exp.y,
                (res.y - exp.y).abs().m(),
                tol.all.m()
            );
        } else {
            println!("Y ok      {} = {} +/- {:e} m", res.y, exp.y, tol.all.m());
        }
        if !res.z.abs_diff_eq(&exp.z, tol.all) {
            println!(
                "Z error: | {} - {} | = {:e} m > {:e} m",
                res.z,
                exp.z,
                (res.z - exp.z).abs().m(),
                tol.all.m()
            );
        } else {
            println!("Z ok      {} = {} +/- {:e} m", res.z, exp.z, tol.all.m());
        }
    }
    is_approx_eq
}
