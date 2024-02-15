use smol_str::SmolStr;

use crate::{
    geodesy::{Ellipsoid, GeodeticDatum},
    operation::{
        conversion::{projection::WebMercator, GeogToGeoc, Normalization},
        identity, Bwd, DynOperation, Fwd, Operation,
    },
};

#[derive(Debug, Clone, PartialEq)]
pub enum Crs {
    /// A [Geocentric] Crs is a **3D cartesian coordinates reference system** in which
    /// coordinates are given by distance **in meters** along the following axes:
    /// - X: axis from the center of the datum's ellipsoid in the equatorial and prime meridian plane,
    /// - Y: axis from the center of the datum's ellipsoid in the equatorial plane and 90 degrees
    /// meridian plane (east)
    /// - Z: axis from the center of the datum's ellipsoid through the north pole.
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
    pub fn id(&self) -> &str {
        match self {
            Crs::Geocentric {
                id,
                datum: _,
                axes: _,
            } => id.as_str(),
            Crs::Geographic {
                id,
                datum: _,
                axes: _,
            } => id.as_str(),
            Crs::Projected {
                id,
                datum: _,
                axes: _,
                projection: _,
            } => id.as_str(),
        }
    }

    pub fn dim(&self) -> usize {
        match self {
            Crs::Geocentric {
                id: _,
                datum: _,
                axes: _,
            } => 3,
            Crs::Geographic {
                id: _,
                datum: _,
                axes,
            } => axes.dim(),
            Crs::Projected {
                id: _,
                datum: _,
                axes,
                projection: _,
            } => axes.dim(),
        }
    }

    /// Return this CRS geodetic datum as a reference.
    pub fn datum(&self) -> &GeodeticDatum {
        match self {
            Crs::Geocentric {
                id: _,
                datum,
                axes: _,
            } => datum,
            Crs::Geographic {
                id: _,
                datum,
                axes: _,
            } => datum,
            Crs::Projected {
                id: _,
                datum,
                axes: _,
                projection: _,
            } => datum,
        }
    }

    pub fn ref_datum_id(&self) -> &str {
        match self {
            Crs::Geocentric {
                id: _,
                datum,
                axes: _,
            } => datum.ref_datum_id(),
            Crs::Geographic {
                id: _,
                datum,
                axes: _,
            } => datum.ref_datum_id(),
            Crs::Projected {
                id: _,
                datum,
                axes: _,
                projection: _,
            } => datum.ref_datum_id(),
        }
    }
    pub fn is_normalized(&self) -> bool {
        match self {
            Crs::Geocentric {
                id: _,
                datum: _,
                axes: _,
            } => true, // TODO: Should we check for the prime meridian
            Crs::Geographic {
                id: _,
                datum: _,
                axes,
            } => {
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
            Crs::Projected {
                id: _,
                datum: _,
                axes,
                projection: _,
            } => {
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

    pub fn normalization(&self) -> Normalization {
        match self {
            Crs::Geocentric {
                id: _,
                datum: _,
                axes,
            } => Normalization::from(*axes),
            Crs::Geographic {
                id: _,
                datum: _,
                axes,
            } => Normalization::from(*axes),
            Crs::Projected {
                id: _,
                datum: _,
                axes,
                projection: _,
            } => Normalization::from(*axes),
        }
    }

    pub fn to_geoc(&self) -> (Box<dyn Operation>, Box<dyn Operation>) {
        match self {
            Crs::Geocentric {
                id: _,
                datum: _,
                axes: _,
            } => (Box::new(identity::<3>()), Box::new(identity::<3>())),
            Crs::Geographic { id: _, datum, axes } => {
                let fwd = Fwd(Normalization::from(*axes)).and_then(Fwd(GeogToGeoc::new(datum)));
                let bwd = Bwd(GeogToGeoc::new(datum)).and_then(Bwd(Normalization::from(*axes)));
                (Box::new(fwd), Box::new(bwd))
            }
            Crs::Projected {
                id: _,
                datum,
                axes,
                projection,
            } => {
                let fwd = Fwd(Normalization::from(*axes))
                    .and_then(Bwd(projection.projection(datum.ellipsoid())))
                    .and_then(Fwd(GeogToGeoc::new(datum)));
                let bwd = Bwd(GeogToGeoc::new(datum))
                    .and_then(Fwd(projection.projection(datum.ellipsoid())))
                    .and_then(Bwd(Normalization::from(*axes)));
                (Box::new(fwd), Box::new(bwd))
            }
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum GeocentricAxes {
    XYZ,
}

/// A [`GeodeticAxes`] value defines the *coordinates system* part of a [`GeodeticCrs`], that is:
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
    WebMercator {
        lon0: f64,
        lat0: f64,
        false_easting: f64,
        false_northing: f64,
    },
}

impl ProjectionSpec {
    pub fn projection(&self, ellipsoid: &Ellipsoid) -> Box<dyn DynOperation> {
        match *self {
            Self::WebMercator {
                lon0,
                lat0,
                false_easting,
                false_northing,
            } => Box::new(WebMercator::new(
                ellipsoid.clone(),
                lon0,
                lat0,
                false_easting,
                false_northing,
            )),
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
