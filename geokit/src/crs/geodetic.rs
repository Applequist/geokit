use super::geocentric::GeocentricCrs;
use super::Crs;
use crate::coord::CoordSpace;
use crate::geodesy::GeodeticDatum;
use crate::id::Id;
use crate::transformation::{Identity, InversibleTransformation};
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

    fn lower(&self) -> Option<(Box<dyn Crs>, Box<dyn InversibleTransformation>)> {
        Some((
            Box::new(GeocentricCrs::new(Id::name("n/a"), self.datum.clone())),
            // FIX: Replace with proper transformation
            Box::new(Identity),
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

    fn lower(&self) -> Option<(Box<dyn Crs>, Box<dyn InversibleTransformation>)> {
        Some((
            Box::new(GeocentricCrs::new(Id::name("n/a"), self.datum.clone())),
            // FIX: Replace with proper transformation
            Box::new(Identity),
        ))
    }
}

#[cfg(test)]
mod tests {
    mod geodetic3d {
        use crate::{
            crs::geodetic::{Geodetic3DAxes, Geodetic3DCrs},
            geodesy::{Ellipsoid, GeodeticDatum, PrimeMeridian},
            id::Id,
        };

        #[test]
        fn partial_eq() {
            let geod3d = Geodetic3DCrs::new(
                Id::name("WGS 84 (geodetic3d)"),
                GeodeticDatum::default(),
                Geodetic3DAxes::EastNorthUp,
                1.,
                1.,
            );

            let different_id = Geodetic3DCrs::new(
                Id::name("WGS 84.1 (geodetic3d)"),
                GeodeticDatum::default(),
                Geodetic3DAxes::EastNorthUp,
                1.,
                1.,
            );
            assert!(!geod3d.eq(&different_id));
            assert!(geod3d.ne(&different_id));
            assert!(!different_id.eq(&geod3d));
            assert!(different_id.ne(&geod3d));

            let different_datum = Geodetic3DCrs::new(
                Id::name("WGS 84 (geodetic3d)"),
                GeodeticDatum::new(
                    Id::name("WGS 84.1"),
                    Ellipsoid::default(),
                    PrimeMeridian::default(),
                    None,
                ),
                Geodetic3DAxes::EastNorthUp,
                1.,
                1.,
            );
            assert!(!geod3d.eq(&different_datum));
            assert!(geod3d.ne(&different_datum));
            assert!(!different_datum.eq(&geod3d));
            assert!(different_datum.ne(&geod3d));

            let different_axes = Geodetic3DCrs::new(
                Id::name("WGS 84 (geodetic3d)"),
                GeodeticDatum::default(),
                Geodetic3DAxes::NorthEastUp,
                1.,
                1.,
            );
            assert!(!geod3d.eq(&different_axes));
            assert!(geod3d.ne(&different_axes));
            assert!(!different_axes.eq(&geod3d));
            assert!(different_axes.ne(&geod3d));

            let different_angle_unit = Geodetic3DCrs::new(
                Id::name("WGS 84 (geodetic3d)"),
                GeodeticDatum::default(),
                Geodetic3DAxes::NorthEastUp,
                1.0_f64.to_radians(), // degres
                1.,
            );
            assert!(!geod3d.eq(&different_angle_unit));
            assert!(geod3d.ne(&different_angle_unit));
            assert!(!different_angle_unit.eq(&geod3d));
            assert!(different_angle_unit.ne(&geod3d));

            let different_height_unit = Geodetic3DCrs::new(
                Id::name("WGS 84 (geodetic3d)"),
                GeodeticDatum::default(),
                Geodetic3DAxes::NorthEastUp,
                1.,   // degres
                0.29, // approx feet
            );
            assert!(!geod3d.eq(&different_angle_unit));
            assert!(geod3d.ne(&different_angle_unit));
            assert!(!different_angle_unit.eq(&geod3d));
            assert!(different_angle_unit.ne(&geod3d));
        }
    }

    mod geodetic2d {}
}
