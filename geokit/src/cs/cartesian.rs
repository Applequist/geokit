//! Cartesian coordinates spaces are used to represent points from 2D and 3D cartesion space.
//!
//! Points in those spaces are represented using distances from hyperplanes along axes.
//! A set of axes defines the direction, orientation and unit to use with these axes.
//!
//! A special set of axes is said to be **normalized**.
//!
//! This module Provide systems of axes for the following cartesian CS:
//! - [GeocentricAxes] for geocentric cartesian CS.
//! - [ProjectedAxes] for projected cartesian CS.
//!
//! And a set of value type to represent **normalized** coordinates.

use crate::{
    math::Float,
    quantities::length::Length,
    units::length::{LengthUnit, M},
};
use approx::AbsDiffEq;
use derive_more::derive::Display;

use super::geodetic::Height;

/// [GeocentricAxes] defines the possible set of axes used in **geocentric** 3D cartesion CS.
/// Geocentric CS uses [meter][crate::units::length::M] unit by default for all axes.
#[derive(Debug, Clone, Copy, PartialEq, Display)]
pub enum GeocentricAxes {
    /// Coordinates are x, y, and z in meters.
    XYZ,
}

/// [XYZ] represents **normalized** geocentric coordinates.
/// - x in meters. **The X-axis is the intersection of the equatirial
/// plane and the Greenwich meridian plane, positive toward the Greenwich meridian.
/// - y in meters, positive east
/// - z in meters, positive north
#[derive(Debug, Copy, Clone, PartialEq, Display)]
#[display("({}, {}, {})", x, y, z)]
pub struct XYZ {
    pub x: Length,
    pub y: Length,
    pub z: Length,
}

impl AbsDiffEq for XYZ {
    type Epsilon = Float;

    fn default_epsilon() -> Self::Epsilon {
        Float::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        self.x.abs_diff_eq(&other.x, epsilon)
            && self.y.abs_diff_eq(&other.y, epsilon)
            && self.z.abs_diff_eq(&other.z, epsilon)
    }
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
