use crate::{coord::Coord3D, geodesy::GeodeticDatum};

/// A `GeocentricCrs` is a **3D cartesian coordinates reference system** in which
/// coordinates are given by distance **in meters** along the following axes:
/// - X: axis from the center of the datum's ellipsoid in the equatorial and prime meridian plane,
/// - Y: axis from the center of the datum's ellipsoid in the equatorial plane and 90 degrees
/// meridian plane (east)
/// - Z: axis from the center of the datum's ellipsoid through the north pole.
pub struct GeocentricCrs {
    /// The geodetic datum of this CRS.
    pub datum: GeodeticDatum,
}

pub enum Geodetic3DAxes {
    /// Coordinates are, in the given order:
    /// - longitude positive east of prime meridian,
    /// - latitude positive north of equatorial plane,
    /// - ellipsoidal height positive upward.
    EastNorthUp,
    /// Coordinates are, in the given order:
    /// - latitude positive north of equatorial plane,
    /// - longitude positive east of prime meridian,
    /// - ellipsoidal height positive upward.
    NorthEastUp,
}

/// A `Geodetic3DCrs` is a **3D geodetic coordinates reference system** in which
/// coordinates are given by longitude, latitude and ellipsoidal height in various order, direction
/// and units.
pub struct Geodetic3DCrs {
    /// The geodetici datum of this CRS.
    pub datum: GeodeticDatum,
    /// Specify the axes order and direction
    pub axes: Geodetic3DAxes,
    /// The angle unit used for longitude and latitude.
    pub angle_unit: f64,
    /// The length unit used for ellipsoidal height.
    pub height_unit: f64,
}

pub enum Geodetic2DAxes {
    /// Coordinates are, in the given order:
    /// - longitude positive east of prime meridian,
    /// - latitude positive north of equatorial plane.
    EastNorth,
    /// Coordinates are, in the given order:
    /// - latitude positive north of equatorial plane,
    /// - longitude positive east of prime meridian.
    NorthEast,
}

/// A `Geodetic2DCrs` is a **2D geodetic coordinates reference system** in which
/// coordinates are given by longitude and latitude in various order, direction
/// and units.
pub struct Geodetic2DCrs {
    /// The geodetici datum of this CRS.
    pub datum: GeodeticDatum,
    /// Specify the axes order and direction
    pub axes: Geodetic2DAxes,
    /// The angle unit used for longitude and latitude.
    pub angle_unit: f64,
}

/// A `TopocentricCrs` is a **3D cartesian coordinates reference system** whose origin is specified
/// as a geodetic location in a base 3D geodetic CRS and axes are derived from the base CRS.
pub struct TopocentricCrs {
    /// The base 3D geodetic CRS.
    base_crs: Geodetic3DCrs,
    /// The origin of the 3D cartesian frame. The axes are derive from the base CRS at the given
    /// location.
    origin: Coord3D,
    /// The length unit used for the coordinates in this CRS.
    unit: f64,
}

/// A dynamic CRS.
pub enum Crs {
    Geocentric(GeocentricCrs),
    Geodetic3D(Geodetic3DCrs),
    Geodetic2D(Geodetic2DCrs),
    Topocentric(TopocentricCrs),
}
