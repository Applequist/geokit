//! This module defines [Crs](Coordinates Reference Systems).
//! A [Crs] ties a coordinates system to the Earth using a [GeodeticDatum](datum) and allows to
//! unambiguously assign coordinates to location on Earth.

use smol_str::SmolStr;

use crate::cs::geodetic::Lon;
use crate::units::length::Meters;
use crate::{
    geodesy::{Ellipsoid, GeodeticDatum},
    operation::{
        conversion::{
            projection::cyl::{Mercator, TransverseMercator, WebMercator},
            GeogToGeoc, Normalization,
        },
        Inv, Operation,
    },
};

#[derive(Debug, Clone, PartialEq)]
pub enum Crs {
    /// A [Geocentric] Crs is a **3D cartesian coordinates reference system** in which
    /// coordinates are given by distance **in meters** along the following axes:
    /// - X: axis from the center of the datum's ellipsoid in the equatorial and prime meridian plane,
    /// - Y: axis from the center of the datum's ellipsoid in the equatorial plane and 90 degrees
    /// meridian plane (east)
    /// - Z: axis from the center of the datum's ellipsoid through the North Pole.
    Geocentric {
        id: SmolStr,
        datum: GeodeticDatum,
        axes: GeocentricAxes,
    },
    /// A `GeographicCrs` is a **2D/3D geodetic coordinates reference system** in which
    /// coordinates are made up of longitude, latitude and optionally ellipsoidal height in various order, direction
    /// and units.
    Geographic {
        id: SmolStr,
        datum: GeodeticDatum,
        axes: GeodeticAxes,
    },
    Projected {
        id: SmolStr,
        datum: GeodeticDatum,
        axes: ProjectedAxes,
        projection: ProjectionSpec,
    },
}

impl Crs {
    /// Returns the id of a [Crs].
    pub fn id(&self) -> &str {
        match self {
            Crs::Geocentric { id, .. } => id.as_str(),
            Crs::Geographic { id, .. } => id.as_str(),
            Crs::Projected { id, .. } => id.as_str(),
        }
    }

    /// Returns the coordinates dimension of this [Crs].
    pub fn dim(&self) -> usize {
        match self {
            Crs::Geocentric { axes: _, .. } => 3,
            Crs::Geographic { axes, .. } => axes.dim(),
            Crs::Projected { axes, .. } => axes.dim(),
        }
    }

    /// Returns the [GeodeticDatum] used by a [Crs].
    pub fn datum(&self) -> &GeodeticDatum {
        match self {
            Crs::Geocentric { datum, .. } => datum,
            Crs::Geographic { datum, .. } => datum,
            Crs::Projected { datum, .. } => datum,
        }
    }

    /// Returns the id of the [GeodeticDatum] used by this [Crs].
    pub fn datum_id(&self) -> &str {
        self.datum().id()
    }

    /// Returns the id of the reference [GeodeticDatum] of a [Crs].
    pub fn ref_datum_id(&self) -> &str {
        match self {
            Crs::Geocentric { datum, .. } => datum.ref_datum_id(),
            Crs::Geographic { datum, .. } => datum.ref_datum_id(),
            Crs::Projected { datum, .. } => datum.ref_datum_id(),
        }
    }

    /// Returns whether a [Crs] is normalized.
    /// A [Crs::Geocentric] CRS is normalized by default.
    /// A [Crs::Geographic] CRS is normalized if it is using the following
    /// axes `GeodeticAxes::EastNorthUp { angle_unit: 1.0, height_unit: 1.0 }`.
    /// And a [Crs::Projected] CRS is normalized if it is using the following axes:
    /// `ProjectedAxes::EastNorthUp { horiz_unit: 1.0, height_unit: 1.0 }`.
    pub fn is_normalized(&self) -> bool {
        match self {
            Crs::Geocentric { .. } => true, // TODO: Should we check for the prime meridian
            Crs::Geographic { axes, .. } => {
                if let GeodeticAxes::EastNorthUp {
                    angle_unit,
                    height_unit,
                } = axes
                {
                    *angle_unit == 1.0 && *height_unit == 1.0
                } else {
                    false
                }
            }
            Crs::Projected { axes, .. } => {
                if let ProjectedAxes::EastNorthUp {
                    horiz_unit,
                    height_unit,
                } = axes
                {
                    *horiz_unit == 1.0 && *height_unit == 1.0
                } else {
                    false
                }
            }
        }
    }

    /// Returns coordinates conversion to and from normalized geocentric coordinates in the
    /// reference datum.
    pub fn to_ref_geoc(&self) -> (Box<dyn Operation>, Box<dyn Operation>) {
        match self {
            Crs::Geocentric { datum, .. } => {
                let op = datum.to_ref_datum();
                (op.clone(), Inv(op).boxed())
            }
            Crs::Geographic { datum, axes, .. } => {
                let op = datum.to_ref_datum();
                let fwd = Normalization::from(*axes)
                    .and_then(GeogToGeoc::new(datum))
                    .and_then(op);
                let bwd = Inv(fwd.clone());
                (fwd.boxed(), bwd.boxed())
            }
            Crs::Projected {
                datum,
                axes,
                projection,
                ..
            } => {
                let op = datum.to_ref_datum();
                let fwd = Normalization::from(*axes)
                    .and_then(Inv(projection.projection(datum.ellipsoid())))
                    .and_then(GeogToGeoc::new(datum))
                    .and_then(op);
                let bwd = Inv(fwd.clone());
                (fwd.boxed(), bwd.boxed())
            }
        }
    }
}

/// The [GeocentricAxes] enum defines the *coordinates system* part of [Geocentric] CRS.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum GeocentricAxes {
    /// Coordinates are x, y, and z in meters.
    XYZ,
}

/// A [`GeodeticAxes`] value defines the *coordinates system* part of a [`GeodeticCrs`], that is:
/// - the ordering and direction of the axes,
/// - the angle unit used for longitude and latitude as a 'radian_per_unit' conversion factor,
/// - the length unit used for the ellipsoidal height as a 'meter_per_unit' conversion factor.
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
    /// Return the dimension (2D or 3D) of the coordinate system.
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

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ProjectedAxes {
    EastNorthUp { horiz_unit: f64, height_unit: f64 },
    EastNorth { horiz_unit: f64 },
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

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ProjectionSpec {
    Mercator1SP {
        lon0: Lon,
        k0: f64,
        false_easting: Meters,
        false_northing: Meters,
    },
    Mercator2SP {
        lon0: Lon,
        lat0: f64,
        false_easting: Meters,
        false_northing: Meters,
    },
    UTMNorth {
        zone: u8,
    },
    UTMSouth {
        zone: u8,
    },
    TransverseMercator {
        lon0: Lon,
        lat0: f64,
        k0: f64,
        false_easting: Meters,
        false_northing: Meters,
    },
    WebMercator {
        lon0: Lon,
        lat0: f64,
        false_easting: Meters,
        false_northing: Meters,
    },
}

impl ProjectionSpec {
    pub fn projection(&self, ellipsoid: &Ellipsoid) -> Box<dyn Operation> {
        match *self {
            Self::Mercator1SP {
                lon0,
                k0,
                false_easting,
                false_northing,
            } => Mercator::new_1_sp(
                ellipsoid,
                lon0.rad(),
                k0,
                false_easting.m(),
                false_northing.m(),
            )
            .boxed(),
            Self::Mercator2SP {
                lon0,
                lat0,
                false_easting,
                false_northing,
            } => Mercator::new_2_sp(
                ellipsoid,
                lon0.rad(),
                lat0,
                false_easting.m(),
                false_northing.m(),
            )
            .boxed(),
            Self::UTMNorth { zone } => TransverseMercator::new_utm_north(ellipsoid, zone).boxed(),
            Self::UTMSouth { zone } => TransverseMercator::new_utm_south(ellipsoid, zone).boxed(),
            Self::TransverseMercator {
                lon0,
                lat0,
                k0,
                false_easting,
                false_northing,
            } => TransverseMercator::new(
                ellipsoid,
                lon0.rad(),
                lat0,
                k0,
                false_easting.m(),
                false_northing.m(),
            )
            .boxed(),
            Self::WebMercator {
                lon0,
                lat0,
                false_easting,
                false_northing,
            } => WebMercator::new(
                ellipsoid,
                lon0.rad(),
                lat0,
                false_easting.m(),
                false_northing.m(),
            )
            .boxed(),
        }
    }
}

#[cfg(test)]
mod tests {

    use super::*;
    use crate::geodesy::*;

    #[test]
    fn clone() {
        let geoc = Crs::Geocentric {
            id: "WGS84".into(),
            datum: GeodeticDatum::new(
                "WGS84",
                ellipsoid::consts::WGS84,
                prime_meridian::consts::GREENWICH,
                None,
            ),
            axes: GeocentricAxes::XYZ,
        };
        let cpy = geoc.clone();
        assert_eq!(geoc, cpy);

        let geod = Crs::Geographic {
            id: "WGS 84 (geodetic2d)".into(),
            datum: GeodeticDatum::new(
                "WGS84",
                ellipsoid::consts::WGS84,
                prime_meridian::consts::GREENWICH,
                None,
            ),
            axes: GeodeticAxes::EastNorthUp {
                angle_unit: 1.,
                height_unit: 1.,
            },
        };

        let cpy = geod.clone();
        assert_eq!(geod, cpy);

        // TODO: Projected
    }

    #[test]
    fn partial_eq() {
        let geoc = Crs::Geocentric {
            id: "WGS84".into(),
            datum: GeodeticDatum::new(
                "WGS84",
                ellipsoid::consts::WGS84,
                prime_meridian::consts::GREENWICH,
                None,
            ),
            axes: GeocentricAxes::XYZ,
        };
        let cpy = geoc.clone();
        assert!(geoc.eq(&cpy));
        assert!(!geoc.ne(&cpy));

        let different_tag = Crs::Geocentric {
            id: "WGS84.1".into(),
            datum: GeodeticDatum::new(
                "WGS84",
                ellipsoid::consts::WGS84,
                prime_meridian::consts::GREENWICH,
                None,
            ),
            axes: GeocentricAxes::XYZ,
        };
        assert_ne!(geoc, different_tag);

        let different_datum = Crs::Geocentric {
            id: "WGS84".into(),
            datum: GeodeticDatum::new(
                "WGS 84.1",
                ellipsoid::consts::GRS80,
                prime_meridian::consts::GREENWICH,
                None,
            ),
            axes: GeocentricAxes::XYZ,
        };
        assert_ne!(geoc, different_datum);

        let geod3d = Crs::Geographic {
            id: "WGS 84 (geodetic3d)".into(),
            datum: GeodeticDatum::new(
                "WGS84",
                ellipsoid::consts::WGS84,
                prime_meridian::consts::GREENWICH,
                None,
            ),
            axes: GeodeticAxes::EastNorthUp {
                angle_unit: 1.0,
                height_unit: 1.0,
            },
        };

        let different_id = Crs::Geographic {
            id: "WGS 84.1 (geodetic3d)".into(),
            datum: GeodeticDatum::new(
                "WGS84",
                ellipsoid::consts::WGS84,
                prime_meridian::consts::GREENWICH,
                None,
            ),
            axes: GeodeticAxes::EastNorthUp {
                angle_unit: 1.0,
                height_unit: 1.0,
            },
        };
        assert!(!geod3d.eq(&different_id));
        assert!(geod3d.ne(&different_id));
        assert!(!different_id.eq(&geod3d));
        assert!(different_id.ne(&geod3d));

        let different_datum = Crs::Geographic {
            id: "WGS 84 (geodetic3d)".into(),
            datum: GeodeticDatum::new(
                "WGS84.1",
                ellipsoid::consts::WGS84,
                prime_meridian::consts::GREENWICH,
                None,
            ),
            axes: GeodeticAxes::EastNorthUp {
                angle_unit: 1.0,
                height_unit: 1.0,
            },
        };
        assert!(!geod3d.eq(&different_datum));
        assert!(geod3d.ne(&different_datum));
        assert!(!different_datum.eq(&geod3d));
        assert!(different_datum.ne(&geod3d));

        let different_axes = Crs::Geographic {
            id: "WGS 84 (geodetic3d)".into(),
            datum: GeodeticDatum::new(
                "WGS84",
                ellipsoid::consts::WGS84,
                prime_meridian::consts::GREENWICH,
                None,
            ),
            axes: GeodeticAxes::NorthEastUp {
                angle_unit: 1.0,
                height_unit: 1.0,
            },
        };
        assert!(!geod3d.eq(&different_axes));
        assert!(geod3d.ne(&different_axes));
        assert!(!different_axes.eq(&geod3d));
        assert!(different_axes.ne(&geod3d));

        let different_angle_unit = Crs::Geographic {
            id: "WGS 84 (geodetic3d)".into(),
            datum: GeodeticDatum::new(
                "WGS84",
                ellipsoid::consts::WGS84,
                prime_meridian::consts::GREENWICH,
                None,
            ),
            axes: GeodeticAxes::NorthEastUp {
                angle_unit: 1.0_f64.to_radians(),
                height_unit: 1.,
            },
        };
        assert!(!geod3d.eq(&different_angle_unit));
        assert!(geod3d.ne(&different_angle_unit));
        assert!(!different_angle_unit.eq(&geod3d));
        assert!(different_angle_unit.ne(&geod3d));

        let different_height_unit = Crs::Geographic {
            id: "WGS 84 (geodetic3d)".into(),
            datum: GeodeticDatum::new(
                "WGS84",
                ellipsoid::consts::WGS84,
                prime_meridian::consts::GREENWICH,
                None,
            ),
            axes: GeodeticAxes::NorthEastUp {
                angle_unit: 1.,
                height_unit: 0.29,
            },
        };
        assert!(!geod3d.eq(&different_height_unit));
        assert!(geod3d.ne(&different_height_unit));
        assert!(!different_height_unit.eq(&geod3d));
        assert!(different_height_unit.ne(&geod3d));
        // TODO: Geographic and projected.
    }
}
