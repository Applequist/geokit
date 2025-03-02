use crate::{
    cs::geodetic::Height,
    math::fp::Float,
    quantities::length::Length,
    units::length::{LengthUnit, M},
};
use approx::AbsDiffEq;
use derive_more::derive::Display;

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
}

#[derive(Copy, Clone, Debug)]
pub struct ProjectedTolerance {
    pub horiz: (Float, LengthUnit),
    pub height: (Float, LengthUnit),
}

impl ProjectedTolerance {
    /// [ProjectedTolerancee::tiny()] allows a `1e-4 m` error on both the
    /// horizontal and the vertical axes.
    /// This translates to a positional (distance) error less than a millimeter.
    pub fn tiny() -> Self {
        ProjectedTolerance {
            horiz: (1e-4, M),
            height: (1e-4, M),
        }
    }

    /// [ProjectedTolerancee::small()] allows a `1e-3 m` maximum error on both the
    /// horizontal and vertial axes.
    /// This translates to a positional (distance) error less than a centimeter.
    pub fn small() -> Self {
        ProjectedTolerance {
            horiz: (1e-3, M),
            height: (1e-3, M),
        }
    }

    /// Returns the maximum error on horizontal axes.
    pub fn horiz(&self) -> Length {
        self.horiz.0 * self.horiz.1
    }

    /// Returns the maximum error on the vertical axis.
    pub fn height(&self) -> Length {
        self.height.0 * self.height.1
    }

    /// Converts the given [ProjectedTolerancee] to the given units.
    pub(crate) fn convert_to(
        &self,
        horiz_unit: LengthUnit,
        height_unit: LengthUnit,
    ) -> ProjectedTolerance {
        let horiz = if self.horiz.1 == horiz_unit {
            self.horiz
        } else {
            (self.horiz().val(horiz_unit), horiz_unit)
        };
        let height = if self.height.1 == height_unit {
            self.height
        } else {
            (self.height().val(height_unit), height_unit)
        };
        ProjectedTolerance { horiz, height }
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
    pub fn approx_eq(self, other: Self, tol: ProjectedTolerance) -> bool {
        self.easting.abs_diff_eq(&other.easting, tol.horiz())
            && self.northing.abs_diff_eq(&other.northing, tol.horiz())
            && self.height.abs_diff_eq(&other.height, tol.height())
    }
}
