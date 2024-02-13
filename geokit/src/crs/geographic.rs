use smol_str::SmolStr;

use crate::{
    geodesy::GeodeticDatum,
    operation::{
        conversion::{GeogToGeoc, Normalization},
        Bwd, Fwd, Operation,
    },
};
use std::fmt::*;

use super::Crs;

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

/// A `GeographicCrs` is a **2D/3D geodetic coordinates reference system** in which
/// coordinates are made up of longitude, latitude and optionally ellipsoidal height in various order, direction
/// and units.
#[derive(Debug, Clone, PartialEq)]
pub struct GeographicCrs {
    id: SmolStr,
    datum: GeodeticDatum,
    axes: GeodeticAxes,
}

impl GeographicCrs {
    /// Create a new [`GeodeticCrs`].
    pub fn new<T: AsRef<str>>(id: T, datum: GeodeticDatum, axes: GeodeticAxes) -> Self {
        Self {
            id: SmolStr::new(id),
            datum,
            axes,
        }
    }

    pub fn id(&self) -> &str {
        self.id.as_str()
    }

    /// Return the coordinates space dimension.
    #[inline]
    pub fn dim(&self) -> usize {
        self.axes.dim()
    }

    /// Return a reference to this CRS [`datum`][GeodeticDatum].
    #[inline]
    pub fn datum(&self) -> &GeodeticDatum {
        &self.datum
    }

    /// Return a reference to the [`axes`][GeodeticAxes] used by this CRS.
    #[inline]
    pub fn axes(&self) -> GeodeticAxes {
        self.axes
    }

    pub fn to_geoc(&self) -> (impl Operation + Clone, impl Operation + Clone) {
        let fwd =
            Fwd(Normalization::from(self.axes())).and_then(Fwd(GeogToGeoc::new(self.datum())));
        let bwd =
            Bwd(GeogToGeoc::new(self.datum())).and_then(Bwd(Normalization::from(self.axes())));
        (fwd, bwd)
    }
}

impl Crs for GeographicCrs {
    fn is_normalized(&self) -> bool {
        if let GeodeticAxes::EastNorthUp {
            angle_unit,
            height_unit,
        } = self.axes()
        {
            angle_unit == 1.0 && height_unit == 1.0
        } else {
            false
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::{
        crs::geographic::{GeodeticAxes, GeographicCrs},
        geodesy::{ellipsoid, prime_meridian, GeodeticDatum},
    };

    #[test]
    fn partial_eq() {
        let geod3d = GeographicCrs::new(
            "WGS 84 (geodetic3d)",
            GeodeticDatum::new(
                "WGS84",
                ellipsoid::consts::WGS84,
                prime_meridian::consts::GREENWICH,
                None,
            ),
            GeodeticAxes::EastNorthUp {
                angle_unit: 1.0,
                height_unit: 1.0,
            },
        );

        let different_id = GeographicCrs::new(
            "WGS 84.1 (geodetic3d)",
            GeodeticDatum::new(
                "WGS84",
                ellipsoid::consts::WGS84,
                prime_meridian::consts::GREENWICH,
                None,
            ),
            GeodeticAxes::EastNorthUp {
                angle_unit: 1.0,
                height_unit: 1.0,
            },
        );
        assert!(!geod3d.eq(&different_id));
        assert!(geod3d.ne(&different_id));
        assert!(!different_id.eq(&geod3d));
        assert!(different_id.ne(&geod3d));

        let different_datum = GeographicCrs::new(
            "WGS 84 (geodetic3d)",
            GeodeticDatum::new(
                "WGS84.1",
                ellipsoid::consts::WGS84,
                prime_meridian::consts::GREENWICH,
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

        let different_axes = GeographicCrs::new(
            "WGS 84 (geodetic3d)",
            GeodeticDatum::new(
                "WGS84",
                ellipsoid::consts::WGS84,
                prime_meridian::consts::GREENWICH,
                None,
            ),
            GeodeticAxes::NorthEastUp {
                angle_unit: 1.0,
                height_unit: 1.0,
            },
        );
        assert!(!geod3d.eq(&different_axes));
        assert!(geod3d.ne(&different_axes));
        assert!(!different_axes.eq(&geod3d));
        assert!(different_axes.ne(&geod3d));

        let different_angle_unit = GeographicCrs::new(
            "WGS 84 (geodetic3d)",
            GeodeticDatum::new(
                "WGS84",
                ellipsoid::consts::WGS84,
                prime_meridian::consts::GREENWICH,
                None,
            ),
            GeodeticAxes::NorthEastUp {
                angle_unit: 1.0_f64.to_radians(),
                height_unit: 1.,
            },
        );
        assert!(!geod3d.eq(&different_angle_unit));
        assert!(geod3d.ne(&different_angle_unit));
        assert!(!different_angle_unit.eq(&geod3d));
        assert!(different_angle_unit.ne(&geod3d));

        let different_height_unit = GeographicCrs::new(
            "WGS 84 (geodetic3d)",
            GeodeticDatum::new(
                "WGS84",
                ellipsoid::consts::WGS84,
                prime_meridian::consts::GREENWICH,
                None,
            ),
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
        let geod = GeographicCrs::new(
            "WGS 84 (geodetic2d)",
            GeodeticDatum::new(
                "WGS84",
                ellipsoid::consts::WGS84,
                prime_meridian::consts::GREENWICH,
                None,
            ),
            GeodeticAxes::EastNorthUp {
                angle_unit: 1.,
                height_unit: 1.,
            },
        );

        let cpy = geod.clone();
        assert_eq!(geod, cpy);
    }
}
