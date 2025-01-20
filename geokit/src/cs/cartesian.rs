use crate::quantity::length::units::LengthUnit;

/// The [GeocentricAxes] enum defines the *coordinates system* part of [GeocentricCrs].
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum GeocentricAxes {
    /// Coordinates are x, y, and z in meters.
    XYZ,
}

/// A [ProjectedAxes] value defines the **coordinates system** part of a [ProjectedCrs].
/// That is:
/// - the ordering and direction of the axes,
/// - the horizontal length unit used for easting and northing,
/// - the height length unit used for the ellipsoidal height.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ProjectedAxes {
    /// The axes for a 3D projected Crs. Coordinates are, in order:
    /// - easting positive eastward,
    /// - northing positive northward,
    /// - ellipsoidal height positive up.
    EastNorthUp {
        /// the length unit for easting and northing, eg
        /// `M` or `US_FT`
        horiz_unit: LengthUnit,
        /// the length unit for easting and northing, eg
        /// `M` or `US_FT`
        height_unit: LengthUnit,
    },
    /// The axes for a 2D projected Crs. Coordinates are, in order:
    /// - easting positive eastward,
    /// - northing positive northward,
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
