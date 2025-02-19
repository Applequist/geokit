//! Cartesian coordinates systems are used to represent points from 2D and 3D cartesion space
//! using distances from hyperplanes along axes as coordinates.
//!
//! Cartesian coordinates systems are used for geocentric, projected and topocentric cartesian
//! spaces.
//!
//! This module Provide systems of axes for the following cartesian CS:
//! - [GeocentricAxes] for geocentric cartesian CS.
//! - [ProjectedAxes] for projected cartesian CS.

use super::geodetic::Height;
use crate::{
    math::fp::Float,
    quantities::length::Length,
    units::length::{LengthUnit, M},
};
use approx::AbsDiffEq;
use derive_more::derive::Display;

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

/// Errors used to compare approximate equality between cartesian coordinates.
///
/// See [XYZ::approx_eq]
/// See [ENH::approx_eq]
pub struct CartesianErrors(Length);

impl CartesianErrors {
    pub fn tiny() -> Self {
        CartesianErrors(Length::tiny())
    }

    pub fn small() -> Self {
        CartesianErrors(Length::small())
    }
}

impl Default for CartesianErrors {
    fn default() -> Self {
        CartesianErrors::tiny()
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
    pub fn dist_to(&self, other: &Self) -> Length {
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
    pub fn approx_eq(&self, other: &Self, err: CartesianErrors) -> bool {
        self.dist_to(other) <= err.0
    }
}

/// Returns whether `res` is approximately equal to `exp` within the given [CartesianErrors] error
/// bounds, printing information about the coordinates when not equal.
///
/// Use only in tests.
pub fn approx_eq_xyz(res: &XYZ, exp: &XYZ, err: &CartesianErrors) -> bool {
    let d = res.dist_to(exp);
    let is_approx_eq = d <= err.0;
    if !is_approx_eq {
        println!("d({}, {}) = {:e} m > {:e} m", res, exp, d.m(), err.0.m());
        if !res.x.abs_diff_eq(&exp.x, err.0) {
            println!(
                "X error: | {} - {} | = {:e} m > {:e} m",
                res.x,
                exp.x,
                (res.x - exp.x).abs().m(),
                err.0.m()
            );
        } else {
            println!("X ok      {} = {} +/- {:e} m", res.x, exp.x, err.0.m());
        }
        if !res.y.abs_diff_eq(&exp.y, err.0) {
            println!(
                "Y error: | {} - {} | = {:e} m > {:e} m",
                res.y,
                exp.y,
                (res.y - exp.y).abs().m(),
                err.0.m()
            );
        } else {
            println!("Y ok      {} = {} +/- {:e} m", res.y, exp.y, err.0.m());
        }
        if !res.z.abs_diff_eq(&exp.z, err.0) {
            println!(
                "Z error: | {} - {} | = {:e} m > {:e} m",
                res.z,
                exp.z,
                (res.z - exp.z).abs().m(),
                err.0.m()
            );
        } else {
            println!("Z ok      {} = {} +/- {:e} m", res.z, exp.z, err.0.m());
        }
    }
    is_approx_eq
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

/// [ProjectedAxes] defines the possible set of axes used in 2D or 3D **projected** cartesian CS.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ProjectedAxes {
    /// The axes for a 3D projected CS. Coordinates are, in order:
    /// - easting positive eastward, using `horiz_unit` [LengthUnit],
    /// - northing positive northward, using `horiz_unit` [LengthUnit],
    /// - ellipsoidal height positive up. using `height_unit` [LengthUnit],
    EastNorthUp {
        /// the length unit for easting and northing, eg
        /// `M` or `US_FT`
        horiz_unit: LengthUnit,
        /// the length unit for easting and northing, eg
        /// `M` or `US_FT`
        height_unit: LengthUnit,
    },
    /// The axes for a 2D projected Crs. Coordinates are, in order:
    /// - easting positive eastward, using `horiz_unit` [LengthUnit],
    /// - northing positive northward, using `horiz_unit` [LengthUnit],
    EastNorth {
        /// the length unit for easting and northing, eg
        /// `M` or `US_FT`
        horiz_unit: LengthUnit,
    },
}

impl ProjectedAxes {
    pub fn normalize(&self, coords: &[Float]) -> ENH {
        match self {
            ProjectedAxes::EastNorthUp {
                horiz_unit,
                height_unit,
            } => ENH {
                easting: Length::new(coords[0], *horiz_unit),
                northing: Length::new(coords[1], *horiz_unit),
                height: Length::new(coords[2], *height_unit),
            },
            ProjectedAxes::EastNorth { horiz_unit } => ENH {
                easting: Length::new(coords[0], *horiz_unit),
                northing: Length::new(coords[1], *horiz_unit),
                height: Length::ZERO,
            },
        }
    }

    pub fn denormalize(&self, enh: ENH, coords: &mut [Float]) {
        match self {
            ProjectedAxes::EastNorthUp {
                horiz_unit,
                height_unit,
            } => {
                coords[0] = enh.easting.val(*horiz_unit);
                coords[1] = enh.northing.val(*horiz_unit);
                coords[2] = enh.height.val(*height_unit);
            }
            ProjectedAxes::EastNorth { horiz_unit } => {
                coords[0] = enh.easting.val(*horiz_unit);
                coords[1] = enh.northing.val(*horiz_unit);
            }
        }
    }

    pub fn dim(&self) -> usize {
        match self {
            ProjectedAxes::EastNorthUp {
                horiz_unit: _,
                height_unit: _,
            } => 3,
            ProjectedAxes::EastNorth { horiz_unit: _ } => 2,
        }
    }
}

/// [ENH] represents **normalized** Projected coordinates.
#[derive(Debug, Copy, Clone, PartialEq, Display)]
#[display("(e = {}, n = {}, h = {})", easting, northing, height)]
pub struct ENH {
    pub easting: Length,
    pub northing: Length,
    pub height: Height,
}
