use super::geocentric::GeocentricCrs;
use super::{CoordSpace, Crs};
use crate::geodesy::{Ellipsoid, GeodeticDatum, PrimeMeridian};
use crate::id::Id;
use crate::transformation::{Result, Transformation};
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

impl GeodeticAxes {
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

type ToOrd = (usize, f64);

#[derive(Debug, Clone)]
struct GeodeticNormalization {
    ord: Vec<ToOrd>,
}

impl GeodeticNormalization {
    /// Create a new [`GeodeticNormalization`] conversion base on [`GeodeticAxes`]
    fn from(value: GeodeticAxes) -> Self {
        let ord = match value {
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
        Self { ord }
    }
}

impl Transformation for GeodeticNormalization {
    fn in_dim(&self) -> usize {
        self.ord.len()
    }

    fn out_dim(&self) -> usize {
        3
    }

    fn apply(&self, input: &[f64], output: &mut [f64]) -> Result<()> {
        output.copy_from_slice(&[0.; 3]);
        for (ox, (ix, f)) in self.ord.iter().enumerate() {
            output[ox] = input[*ix] * f;
        }
        Ok(())
    }

    fn boxed(&self) -> Box<dyn Transformation> {
        Box::new(self.clone())
    }
}

#[derive(Debug, Clone)]
struct GeodeticDenormalization {
    ord: Vec<ToOrd>,
}

impl GeodeticDenormalization {
    /// Create a new [`GeodeticDenormalization`] based on [`GeodeticAxes`].
    fn from(axes: GeodeticAxes) -> Self {
        let ord = match axes {
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
        Self { ord }
    }
}

impl Transformation for GeodeticDenormalization {
    fn in_dim(&self) -> usize {
        3
    }

    fn out_dim(&self) -> usize {
        self.ord.len()
    }

    fn apply(&self, input: &[f64], output: &mut [f64]) -> Result<()> {
        for (ox, (ix, f)) in self.ord.iter().enumerate() {
            output[ox] = input[*ix] * f;
        }
        Ok(())
    }

    fn boxed(&self) -> Box<dyn Transformation> {
        Box::new(self.clone())
    }
}

/// [`GeodeticToGeocentric`] converts **normalized geodetic coordinates** to **normalized
/// geocentric coordinates**.
///
/// As mentioned in the EPSG guidance note 7 part 2, paragraph 4.1.1, this Transformation
/// first transform to/from non-greenwich based geodetic coordinates into greenwich based
/// ones.
#[derive(Debug, Clone, Copy)]
struct GeodeticToGeocentric {
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
        let llh_gw = [
            self.prime_meridian.to_greenwich(input[0]),
            input[1],
            input[2],
        ];

        self.ellipsoid.llh_to_xyz(&llh_gw, output);
        Ok(())
    }

    fn boxed(&self) -> Box<dyn Transformation> {
        Box::new(*self)
    }
}

#[derive(Debug, Clone, Copy)]
struct GeocentricToGeodetic {
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
        output[0] = self.prime_meridian.from_greenwich(output[0]);
        Ok(())
    }

    fn boxed(&self) -> Box<dyn Transformation> {
        Box::new(*self)
    }
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

impl Default for GeodeticCrs {
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
            GeodeticNormalization::from(self.axes).boxed(),
            GeodeticDenormalization::from(self.axes).boxed(),
        )
    }

    fn lowered(
        &self,
    ) -> Option<(
        Box<dyn Crs>,
        Box<dyn Transformation>,
        Box<dyn Transformation>,
    )> {
        Some((
            Box::new(GeocentricCrs::new(
                Id::name(format!("Lowered from {}", self.id)),
                self.datum.clone(),
            )),
            GeodeticNormalization::from(self.axes)
                .and_then(GeodeticToGeocentric::from(&self.datum))
                .boxed(),
            GeocentricToGeodetic::from(&self.datum)
                .and_then(GeodeticDenormalization::from(self.axes))
                .boxed(),
        ))
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
