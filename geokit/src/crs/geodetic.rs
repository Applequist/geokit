use super::geocentric::GeocentricCrs;
use super::{CoordSpace, Crs, LoweringTransformation};
use crate::geodesy::GeodeticDatum;
use crate::id::Id;
use crate::transformation::{inv, Chain, CoordScaling, GeodToGeoc, Inv, Transformation};
use std::fmt::*;

/// A [`GeodeticAxes`] value defines the *coordinates system* part of a [`GedoeticCrs`], that is:
/// - the ordering and direction of the axes,
/// - the angle unit used for longitude and latitude,
/// - the length unit used for the ellipsoidal height.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum GeodeticAxes {
    /// Coordinates are given in the following order:
    /// - longitude positive east of prime meridian,
    /// - latitude positive north of equatorial plane,
    /// - ellipsoidal height positive upward.
    EastNorthUp {
        angle_unit: f64,
        height_unit: f64,
    },
    /// Coordinates are given in the following order:
    /// - latitude positive north of equatorial plane,
    /// - longitude positive east of prime meridian,
    /// - ellipsoidal height positive upward.
    NorthEastUp {
        angle_unit: f64,
        height_unit: f64,
    },
    /// Coordinates are, in the given order:
    /// - longitude positive east of prime meridian,
    /// - latitude positive north of equatorial plane.
    EastNorth {
        angle_unit: f64,
    },
    /// Coordinates are, in the given order:
    /// - latitude positive north of equatorial plane,
    /// - longitude positive east of prime meridian.
    NorthEast {
        angle_unit: f64,
    },
    NorthWest {
        angle_unit: f64,
    },
}

impl GeodeticAxes {
    /// Return the dimension (2D or 3D) of the coordinates system.
    pub fn dim(&self) -> usize {
        match self {
            GeodeticAxes::EastNorthUp {
                angle_unit: _,
                height_unit: _,
            }
            | GeodeticAxes::NorthEastUp {
                angle_unit: _,
                height_unit: _,
            } => 3,
            GeodeticAxes::EastNorth { angle_unit: _ }
            | GeodeticAxes::NorthEast { angle_unit: _ }
            | GeodeticAxes::NorthWest { angle_unit: _ } => 2,
        }
    }

    /// Return the **coordinates normalization** conversion.
    pub fn normalization(&self) -> CoordScaling {
        let to_ord = match *self {
            GeodeticAxes::EastNorthUp {
                angle_unit,
                height_unit,
            } => vec![(0, angle_unit), (1, angle_unit), (2, height_unit)],
            GeodeticAxes::EastNorth { angle_unit } => vec![(0, angle_unit), (1, angle_unit)],
            GeodeticAxes::NorthEastUp {
                angle_unit,
                height_unit,
            } => vec![(1, angle_unit), (0, angle_unit), (2, height_unit)],
            GeodeticAxes::NorthEast { angle_unit } => vec![(1, angle_unit), (0, angle_unit)],
            GeodeticAxes::NorthWest { angle_unit } => vec![(1, angle_unit), (0, -angle_unit)],
        };
        CoordScaling(to_ord)
    }
}

impl Default for GeodeticAxes {
    /// Return the [`GeodeticAxes`] used in **normalized geodetic coordinates**.
    fn default() -> Self {
        Self::EastNorthUp {
            angle_unit: 1.0,
            height_unit: 1.0,
        }
    }
}

/// A `GeodeticCrs` is a **2D/3D geodetic coordinates reference system** in which
/// coordinates are made up of longitude, latitude and optionally ellipsoidal height in various order, direction
/// and units.
#[derive(Debug, Clone, PartialEq)]
pub struct GeodeticCrs {
    id: Id,
    datum: GeodeticDatum,
    axes: GeodeticAxes,
}

impl GeodeticCrs {
    /// Create a new [`GeodeticCrs`].
    pub fn new(id: Id, datum: GeodeticDatum, axes: GeodeticAxes) -> Self {
        Self { id, datum, axes }
    }

    /// Return a reference to this CRS [`datum`][GeodeticDatum].
    pub fn datum(&self) -> &GeodeticDatum {
        &self.datum
    }

    /// Return a reference to the [`axes`][GeodeticAxes] used by this CRS.
    pub fn axes(&self) -> &GeodeticAxes {
        &self.axes
    }

    /// Return the *lower* [`geocentric`][GeocentricCrs] CRS derived from this geodetic CRS
    /// and the coordinates conversion to and from this geocentric CRS.
    pub fn geoc_crs(
        &self,
    ) -> (
        GeocentricCrs,
        Chain<CoordScaling, GeodToGeoc>,
        Chain<Inv<GeodToGeoc>, Inv<CoordScaling>>,
    ) {
        let geoc_crs = GeocentricCrs::new(
            self.id.renamed(format!("{} (as geocentric)", self.id())),
            // FIX: Should we enforce Greewich prime meridian ? See [GeocentricCrs].
            self.datum.clone(),
        );
        let normalization = self.axes.normalization();
        let geod_to_geoc = GeodToGeoc::new(&self.datum);
        (
            geoc_crs,
            normalization.clone().and_then(geod_to_geoc),
            inv(geod_to_geoc).and_then(inv(normalization)),
        )
    }
}

impl Default for GeodeticCrs {
    /// Return EPSG:4326 as default the default geodetic CRS.
    fn default() -> Self {
        GeodeticCrs::new(
            Id::full("WGS 84", "EPSG", 4326),
            GeodeticDatum::default(),
            GeodeticAxes::NorthEast {
                angle_unit: 1.0_f64.to_radians(),
            },
        )
    }
}

impl Crs for GeodeticCrs {
    fn id(&self) -> &Id {
        &self.id
    }

    fn dim(&self) -> usize {
        self.axes.dim()
    }

    fn kind(&self) -> CoordSpace {
        CoordSpace::Geodetic
    }

    fn datum(&self) -> &GeodeticDatum {
        &self.datum
    }

    fn lowered(&self) -> Option<LoweringTransformation> {
        let (crs, to, from) = self.geoc_crs();
        Some((Box::new(crs), to.boxed(), from.boxed()))
    }
}

#[cfg(test)]
mod tests {
    use crate::{
        crs::{
            geodetic::{GeodeticAxes, GeodeticCrs},
            Crs,
        },
        geodesy::{Ellipsoid, GeodeticDatum, PrimeMeridian},
        id::Id,
    };

    #[test]
    fn deault() {
        let wgs84_2d = GeodeticCrs::default();
        assert_eq!(wgs84_2d.dim(), 2);
        assert_eq!(wgs84_2d.datum(), &GeodeticDatum::default());
    }

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

    #[test]
    fn lowered() {
        unimplemented!();
    }
}
