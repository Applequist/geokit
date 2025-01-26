//! This module defines [Crs](Coordinates Reference Systems).
//! A [Crs] ties a coordinates system to the Earth using a [GeodeticDatum](datum) and allows to
//! unambiguously assign coordinates to location on Earth.

use crate::cs::cartesian::{GeocentricAxes, ProjectedAxes};
use crate::cs::geodetic::GeodeticAxes;
use crate::geodesy::GeodeticDatum;
use crate::projections::ProjectionSpec;
use smol_str::SmolStr;

/// A [GeocentricCrs] is a **3D cartesian coordinates reference system** in which
/// coordinates are given by distance **in meters** along the following axes:
/// - X: axis from the center of the datum's ellipsoid in the equatorial and prime meridian plane,
/// - Y: axis from the center of the datum's ellipsoid in the equatorial plane and 90 degrees
/// meridian plane (east)
/// - Z: axis from the center of the datum's ellipsoid through the North Pole.
#[derive(Debug, Clone, PartialEq)]
pub struct GeocentricCrs {
    pub id: SmolStr,
    pub datum: GeodeticDatum,
    pub axes: GeocentricAxes,
}

/// A [GeographicCrs] is a **2D or 3D geodetic coordinates reference system** in which
/// coordinates are made up of longitude, latitude and optionally ellipsoidal height in various order, direction
/// and quantity.
#[derive(Debug, Clone, PartialEq)]
pub struct GeographicCrs {
    pub id: SmolStr,
    pub datum: GeodeticDatum,
    pub axes: GeodeticAxes,
}

/// A [ProjectedCrs] is a **2D/3D cartesian coordinates reference system** derived from a
/// 2D or 3D [Geographic] Crs using a projection.
/// Coordinates are made up of easting, northing and optionally ellipsoidal height in various order,
/// direction and quantity.
#[derive(Debug, Clone, PartialEq)]
pub struct ProjectedCrs {
    pub id: SmolStr,
    pub datum: GeodeticDatum,
    pub axes: ProjectedAxes,
    pub projection: ProjectionSpec,
}
