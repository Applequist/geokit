use approx::AbsDiffEq;
use derive_more::derive::Display;

use crate::{
    cs::geodetic::Height, math::fp::Float, quantities::length::Length, units::length::LengthUnit,
};

use super::CartesianErrors;

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

pub struct ProjectedErrors {
    pub horiz: (Float, LengthUnit),
    pub height: (Float, LengthUnit),
}

impl ProjectedErrors {
    pub fn horiz(&self) -> Length {
        self.horiz.0 * self.horiz.1
    }

    pub fn height(&self) -> Length {
        self.height.0 * self.height.1
    }

    pub(crate) fn convert_to(&self, horiz_unit: LengthUnit, height_unit: LengthUnit) -> [Float; 3] {
        todo!()
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

impl ENH {
    /// Computes the distance between `self` and `other` in meters.
    pub fn dist_to(self, other: Self) -> Length {
        (self.easting - other.easting)
            .hypot(self.northing - other.northing)
            .hypot(self.height - other.height)
    }

    /// Checks whether this [ENH] is approximately equal to the `other` [ENH]
    /// within the given [ProjectedErrors] bounds.
    pub fn approx_eq(self, other: Self, err: ProjectedErrors) -> bool {
        self.easting.abs_diff_eq(&other.easting, err.horiz())
            && self.northing.abs_diff_eq(&other.northing, err.horiz())
            && self.height.abs_diff_eq(&other.height, err.height())
    }
}
