use super::geocentric::GeocentricCrs;
use super::{CoordSpace, Crs};
use crate::geodesy::{Ellipsoid, GeodeticDatum, PrimeMeridian};
use crate::id::Id;
use crate::transformation::{Chain, CoordScaling, Result, Transformation};
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

impl GeodeticAxes {
    /// Return the dimension (2D or 3D) of the coordinates system.
    fn dim(&self) -> usize {
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
            | GeodeticAxes::NorthEast { angle_unit: _ } => 2,
        }
    }

    /// Return the **coordinates normalization** conversion.
    fn normalization(&self) -> CoordScaling {
        let ord = match *self {
            GeodeticAxes::EastNorthUp {
                angle_unit,
                height_unit,
            } => {
                vec![(0, angle_unit), (1, angle_unit), (2, height_unit)]
            }
            GeodeticAxes::EastNorth { angle_unit } => {
                vec![(0, angle_unit), (1, angle_unit)]
            }
            GeodeticAxes::NorthEastUp {
                angle_unit,
                height_unit,
            } => {
                vec![(1, angle_unit), (0, angle_unit), (2, height_unit)]
            }
            GeodeticAxes::NorthEast { angle_unit } => {
                vec![(1, angle_unit), (0, angle_unit)]
            }
        };
        CoordScaling(ord)
    }

    /// Return the **coordinates denormalization** conversion.
    fn denormalization(&self) -> CoordScaling {
        let ord = match *self {
            GeodeticAxes::EastNorthUp {
                angle_unit,
                height_unit,
            } => {
                let from_rad = 1. / angle_unit;
                let from_meters = 1. / height_unit;
                vec![(0, from_rad), (1, from_rad), (2, from_meters)]
            }
            GeodeticAxes::EastNorth { angle_unit } => {
                let from_rad = 1. / angle_unit;
                vec![(0, from_rad), (1, from_rad)]
            }
            GeodeticAxes::NorthEastUp {
                angle_unit,
                height_unit,
            } => {
                let from_rad = 1. / angle_unit;
                let from_meters = 1. / height_unit;
                vec![(1, from_rad), (0, from_rad), (2, from_meters)]
            }
            GeodeticAxes::NorthEast { angle_unit } => {
                let from_rad = 1. / angle_unit;
                vec![(1, from_rad), (0, from_rad)]
            }
        };
        CoordScaling(ord)
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

/// [`GeodeticToGeocentric`] converts **normalized geodetic coordinates** to **normalized
/// geocentric coordinates**.
///
/// **IMPORTANT**: As mentioned in the EPSG guidance note 7 part 2, paragraph 4.1.1, this Transformation
/// first transform **normalized geodetic coordinates** to greenwich based normalized geodetic coordinates.
#[derive(Debug, Clone, Copy)]
pub struct GeodeticToGeocentric {
    ellipsoid: Ellipsoid,
    prime_meridian: PrimeMeridian,
}

impl GeodeticToGeocentric {
    /// Create a new [`GeodeticToGeocentric`] conversion from the given [`datum`][GeodeticDatum].
    fn from(datum: &GeodeticDatum) -> Self {
        Self {
            ellipsoid: datum.ellipsoid(),
            prime_meridian: datum.prime_meridian(),
        }
    }
}

impl Transformation for GeodeticToGeocentric {
    fn in_dim(&self) -> usize {
        3
    }

    fn out_dim(&self) -> usize {
        3
    }

    fn apply(&self, input: &[f64], output: &mut [f64]) -> Result<()> {
        // convert longitude to Greenwich-base longitude
        // FIX: what if we cross the antimeridian ?
        let llh_gw = [input[0] + self.prime_meridian.lon(), input[1], input[2]];
        self.ellipsoid.llh_to_xyz(&llh_gw, output);
        Ok(())
    }

    fn boxed(&self) -> Box<dyn Transformation> {
        Box::new(*self)
    }
}

/// [`GeocentricToGeodetic`] converts **normalized geocentric coordinates** to **normalized
/// geodetic coordinates**.
///
/// **IMPORTANT**: As mentioned in the EPSG guidance note 7 part 2, paragraph 4.1.1, this Transformation
/// also transform the normalized geodetic coordinates to non-greenwich based normalized geodetic coordinates.
#[derive(Debug, Clone, Copy)]
pub struct GeocentricToGeodetic {
    ellipsoid: Ellipsoid,
    prime_meridian: PrimeMeridian,
}

impl GeocentricToGeodetic {
    pub fn from(datum: &GeodeticDatum) -> Self {
        Self {
            ellipsoid: datum.ellipsoid(),
            prime_meridian: datum.prime_meridian(),
        }
    }
}

impl Transformation for GeocentricToGeodetic {
    fn in_dim(&self) -> usize {
        3
    }

    fn out_dim(&self) -> usize {
        3
    }

    fn apply(&self, input: &[f64], output: &mut [f64]) -> Result<()> {
        self.ellipsoid.xyz_to_llh(input, output);
        // Convert longitude from Greenwich-based longitude to this prime_meridian base longitude.
        // FIX: What if we cross the anti-meridian
        output[0] -= self.prime_meridian.lon();
        Ok(())
    }

    fn boxed(&self) -> Box<dyn Transformation> {
        Box::new(*self)
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
        Chain<CoordScaling, GeodeticToGeocentric>,
        Chain<GeocentricToGeodetic, CoordScaling>,
    ) {
        let geoc_crs = GeocentricCrs::new(
            self.id.renamed(format!("{} (as geocentric)", self.id())),
            // FIX: Should we enforce Greewich prime meridian ? See [GeocentricCrs].
            self.datum.clone(),
        );
        (
            geoc_crs,
            self.axes
                .normalization()
                .and_then(GeodeticToGeocentric::from(&self.datum)),
            GeocentricToGeodetic::from(&self.datum).and_then(self.axes.denormalization()),
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

    fn normalized(
        &self,
    ) -> (
        Box<dyn Crs>,
        Box<dyn Transformation>,
        Box<dyn Transformation>,
    ) {
        (
            Box::new(GeodeticCrs::new(
                Id::name(format!("lowered from {}", self.id)),
                self.datum.clone(),
                GeodeticAxes::default(),
            )),
            self.axes.normalization().boxed(),
            self.axes.denormalization().boxed(),
        )
    }

    fn lowered(
        &self,
    ) -> Option<(
        Box<dyn Crs>,
        Box<dyn Transformation>,
        Box<dyn Transformation>,
    )> {
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
    fn normalized() {
        let geod3d = GeodeticCrs::new(
            Id::name("WGS 84 (geodetic 3d)"),
            GeodeticDatum::default(),
            GeodeticAxes::NorthEastUp {
                angle_unit: 1.0_f64.to_radians(), // degrees
                height_unit: 1.,                  // meters
            },
        );

        let (ngeod3d, n, d) = geod3d.normalized();
        assert_eq!(ngeod3d.datum(), geod3d.datum());

        let p = [45., 90., 100.]; // lat(deg), lon(deg), h(m)
        let mut norm_p = [0.; 3]; // lon(rad), lat(rad), h(m)
        n.apply(&p, &mut norm_p).unwrap();
        assert!(
            (norm_p[0] - std::f64::consts::FRAC_PI_2).abs() < 1e-5
                && (norm_p[1] - std::f64::consts::FRAC_PI_4).abs() < 1e-5
                && norm_p[2] == 100.
        );

        let mut dnorm_p = [0.; 3];
        d.apply(&norm_p, &mut dnorm_p).unwrap();
        assert!(
            (dnorm_p[0] - p[0]).abs() < 1e-5
                && (dnorm_p[1] - p[1]).abs() < 1e-5
                && (dnorm_p[2] - p[2]).abs() < 1e-3
        );
    }

    #[test]
    fn lowered() {
        unimplemented!();
    }
}
