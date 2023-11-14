use super::geocentric::GeocentricCrs;
use super::{Crs, LowerTransformation};
use crate::coord::CoordSpace;
use crate::geodesy::GeodeticDatum;
use crate::id::Id;
use std::fmt::*;

/// Defines the ordering and direction of the axes of a 3D geodetic CRS.
#[derive(Debug, Copy, Clone, PartialEq)]
pub enum Geodetic3DAxes {
    /// Coordinates are given in the following order:
    /// - longitude positive east of prime meridian,
    /// - latitude positive north of equatorial plane,
    /// - ellipsoidal height positive upward.
    EastNorthUp,
    /// Coordinates are given in the following order:
    /// - latitude positive north of equatorial plane,
    /// - longitude positive east of prime meridian,
    /// - ellipsoidal height positive upward.
    NorthEastUp,
}

/// A `Geodetic3DCrs` is a **3D geodetic coordinates reference system** in which
/// coordinates are given by longitude, latitude and ellipsoidal height in various order, direction
/// and units.
#[derive(Debug, Clone, PartialEq)]
pub struct Geodetic3DCrs {
    id: Id,
    datum: GeodeticDatum,
    axes: Geodetic3DAxes,
    angle_unit: f64,
    height_unit: f64,
}

impl Geodetic3DCrs {
    pub fn new(
        id: Id,
        datum: GeodeticDatum,
        axes: Geodetic3DAxes,
        angle_unit: f64,
        height_unit: f64,
    ) -> Self {
        Self {
            id,
            datum,
            axes,
            angle_unit,
            height_unit,
        }
    }
}

impl Crs for Geodetic3DCrs {
    fn id(&self) -> &Id {
        &self.id
    }

    fn dim(&self) -> usize {
        3
    }

    fn kind(&self) -> CoordSpace {
        CoordSpace::Geodetic
    }

    fn datum(&self) -> &GeodeticDatum {
        &self.datum
    }

    fn lower(&self, transformation: LowerTransformation) -> Option<(Box<dyn Crs>, String)> {
        let t = match transformation {
            LowerTransformation::TO => String::from("Normalization + norm_geod3d to norm_geoc"),
            LowerTransformation::FROM => String::from("norm_geoc to norm_geod3d + Denormalization"),
        };
        Some((
            Box::new(GeocentricCrs::new(Id::name("n/a"), self.datum.clone())),
            t,
        ))
    }
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
    id: Id,
    datum: GeodeticDatum,
    axes: Geodetic2DAxes,
    angle_unit: f64,
}

impl Geodetic2DCrs {
    /// Create a new 2D geodetic CRS.
    pub fn new(id: Id, datum: GeodeticDatum, axes: Geodetic2DAxes, angle_unit: f64) -> Self {
        Self {
            id,
            datum,
            axes,
            angle_unit,
        }
    }
}

impl Crs for Geodetic2DCrs {
    fn id(&self) -> &Id {
        &self.id
    }

    fn dim(&self) -> usize {
        2
    }

    fn kind(&self) -> CoordSpace {
        CoordSpace::Geodetic
    }

    fn datum(&self) -> &GeodeticDatum {
        &self.datum
    }

    fn lower(&self, transformation: LowerTransformation) -> Option<(Box<dyn Crs>, String)> {
        let t = match transformation {
            LowerTransformation::TO => String::from("Normalization + norm_geod2d to norm_geoc"),
            LowerTransformation::FROM => String::from("norm_geoc to norm_geod2d + Denormalization"),
        };
        Some((
            Box::new(GeocentricCrs::new(Id::name("n/a"), self.datum.clone())),
            t,
        ))
    }
}
