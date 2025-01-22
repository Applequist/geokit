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

use crate::units::length::LengthUnit;

use super::r1::Length;

/// [GeocentricAxes] defines the possible set of axes used in **geocentric** 3D cartesion CS.
/// Geocentric CS uses [meter][crate::units::length::M] unit by default for all axes.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum GeocentricAxes {
    /// Coordinates are x, y, and z in meters.
    XYZ,
}

/// First normalized geocentric coordinates type.
pub type X = Length;
/// Second normalized geocentric coordinates type.
pub type Y = Length;
/// Third normalized geocentric coordinates type.
pub type Z = Length;

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

/// First normalized projected coordinates type
pub type Easting = Length;
/// Second normalized projected coordinates type
pub type Northing = Length;
