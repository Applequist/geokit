use super::geocentric::GeocentricCrs;
use super::{CoordSpace, Crs};
use crate::geodesy::GeodeticDatum;
use crate::id::Id;
use crate::transformation::{Identity, InvertibleTransformation, Transformation};
use std::fmt::*;

/// Defines the ordering and direction of the axes of a 3D geodetic CRS.
#[derive(Debug, Copy, Clone, PartialEq)]
pub enum GeodeticAxes {
    /// Coordinates are given in the following order:
    /// - longitude positive east of prime meridian,
    /// - latitude positive north of equatorial plane,
    /// - ellipsoidal height positive upward.
    EastNorthUp { angle_unit: f64, height_unit: f64 },
    /// Coordinates are given in the following order:
    /// - latitude positive north of equatorial plane,
    /// - longitude positive east of prime meridian,
    /// - ellipsoidal height positive upward.
    NorthEastUp { angle_unit: f64, height_unit: f64 },
    /// Coordinates are, in the given order:
    /// - longitude positive east of prime meridian,
    /// - latitude positive north of equatorial plane.
    EastNorth { angle_unit: f64 },
    /// Coordinates are, in the given order:
    /// - latitude positive north of equatorial plane,
    /// - longitude positive east of prime meridian.
    NorthEast { angle_unit: f64 },
}

/// A `GeodeticCrs` is a **2D/3D geodetic coordinates reference system** in which
/// coordinates are given by longitude, latitude and optionally ellipsoidal height in various order, direction
/// and units.
#[derive(Debug, Clone, PartialEq)]
pub struct GeodeticCrs {
    id: Id,
    datum: GeodeticDatum,
    axes: GeodeticAxes,
}

impl GeodeticCrs {
    pub fn new(id: Id, datum: GeodeticDatum, axes: GeodeticAxes) -> Self {
        Self { id, datum, axes }
    }
}

impl Crs for GeodeticCrs {
    fn id(&self) -> &Id {
        &self.id
    }

    fn dim(&self) -> usize {
        match self.axes {
            GeodeticAxes::EastNorthUp {
                angle_unit: _,
                height_unit: _,
            }
            | GeodeticAxes::NorthEastUp {
                angle_unit: _,
                height_unit: _,
            } => 3,
            GeodeticAxes::EastNorth { angle_unit: _ }
            | GeodeticAxes::NorthEast { angle_unit: _ } => 2,
        }
    }

    fn kind(&self) -> CoordSpace {
        CoordSpace::Geodetic
    }

    fn datum(&self) -> &GeodeticDatum {
        &self.datum
    }

    fn lower(&self) -> Option<(Box<dyn Crs>, Box<dyn InvertibleTransformation>)> {
        Some((
            Box::new(GeocentricCrs::new(
                Id::name(format!("Lowered from {}", self.id)),
                self.datum.clone(),
            )),
            // FIX: Replace with proper transformation
            Box::new(Identity::<3>),
        ))
    }

    fn normalization(&self) -> Box<dyn Transformation> {
        // FIX: Replace with proper transformation
        Box::new(Identity::<3>)
    }

    fn denormalization(&self) -> Box<dyn Transformation> {
        // FIX: Replace with proper transformation
        Box::new(Identity::<3>)
    }
}

#[cfg(test)]
mod tests {
    use crate::{
        crs::geodetic::{GeodeticAxes, GeodeticCrs},
        geodesy::{Ellipsoid, GeodeticDatum, PrimeMeridian},
        id::Id,
    };

    #[test]
    fn partial_eq() {
        let geod3d = GeodeticCrs::new(
            Id::name("WGS 84 (geodetic3d)"),
            GeodeticDatum::default(),
            GeodeticAxes::EastNorthUp {
                angle_unit: 1.0,
                height_unit: 1.0,
            },
        );

        let different_id = GeodeticCrs::new(
            Id::name("WGS 84.1 (geodetic3d)"),
            GeodeticDatum::default(),
            GeodeticAxes::EastNorthUp {
                angle_unit: 1.0,
                height_unit: 1.0,
            },
        );
        assert!(!geod3d.eq(&different_id));
        assert!(geod3d.ne(&different_id));
        assert!(!different_id.eq(&geod3d));
        assert!(different_id.ne(&geod3d));

        let different_datum = GeodeticCrs::new(
            Id::name("WGS 84 (geodetic3d)"),
            GeodeticDatum::new(
                Id::name("WGS 84.1"),
                Ellipsoid::default(),
                PrimeMeridian::default(),
                None,
            ),
            GeodeticAxes::EastNorthUp {
                angle_unit: 1.0,
                height_unit: 1.0,
            },
        );
        assert!(!geod3d.eq(&different_datum));
        assert!(geod3d.ne(&different_datum));
        assert!(!different_datum.eq(&geod3d));
        assert!(different_datum.ne(&geod3d));

        let different_axes = GeodeticCrs::new(
            Id::name("WGS 84 (geodetic3d)"),
            GeodeticDatum::default(),
            GeodeticAxes::NorthEastUp {
                angle_unit: 1.0,
                height_unit: 1.0,
            },
        );
        assert!(!geod3d.eq(&different_axes));
        assert!(geod3d.ne(&different_axes));
        assert!(!different_axes.eq(&geod3d));
        assert!(different_axes.ne(&geod3d));

        let different_angle_unit = GeodeticCrs::new(
            Id::name("WGS 84 (geodetic3d)"),
            GeodeticDatum::default(),
            GeodeticAxes::NorthEastUp {
                angle_unit: 1.0_f64.to_radians(),
                height_unit: 1.,
            },
        );
        assert!(!geod3d.eq(&different_angle_unit));
        assert!(geod3d.ne(&different_angle_unit));
        assert!(!different_angle_unit.eq(&geod3d));
        assert!(different_angle_unit.ne(&geod3d));

        let different_height_unit = GeodeticCrs::new(
            Id::name("WGS 84 (geodetic3d)"),
            GeodeticDatum::default(),
            GeodeticAxes::NorthEastUp {
                angle_unit: 1.,
                height_unit: 0.29,
            },
        );
        assert!(!geod3d.eq(&different_height_unit));
        assert!(geod3d.ne(&different_height_unit));
        assert!(!different_height_unit.eq(&geod3d));
        assert!(different_height_unit.ne(&geod3d));
    }

    #[test]
    fn clone() {
        let geod = GeodeticCrs::new(
            Id::name("WGS 84 (geodetic2d)"),
            GeodeticDatum::default(),
            GeodeticAxes::EastNorthUp {
                angle_unit: 1.,
                height_unit: 1.,
            },
        );

        let cpy = geod.clone();
        assert_eq!(geod, cpy);
    }
}
