#![allow(non_snake_case)]

use super::{Ellipsoid, PrimeMeridian};
use crate::cs::cartesian::XYZ;
use crate::cs::geodetic::{Lat, Lon, LLH};
use crate::math::fp::PI_2;
use crate::units::angle::RAD;
use crate::units::length::M;
use derive_more::derive::Display;
use nalgebra::ComplexField;
use smol_str::SmolStr;
use std::fmt::Debug;

/// A `datum` is the information required to attach a coordinate system to an object.
/// A [GeodeticDatum] is a `datum` describing the relationship of an ellipsoidal model of the Earth
/// with the real Earth.
/// It is defined by an id, an [Ellipsoid] and a [PrimeMeridian].
///
/// # Equality
///
/// 2 datums are considered equal if they have **the same id**, the same ellipsoid and prime meridian.
/// That is: it is not enough to use the same ellipsoid and prime meridian for 2 datums to be considered
/// the same.
#[derive(Debug, Clone, Display)]
#[display("(id = {}, ellps = {}, pm = {})", id, ellipsoid, prime_meridian)]
pub struct GeodeticDatum {
    id: SmolStr,
    ellipsoid: Ellipsoid,
    prime_meridian: PrimeMeridian,
}

impl GeodeticDatum {
    /// Creates a new [`GeodeticDatum`].
    /// Used to define well known datum as constant.
    pub const fn new_static(
        id: &'static str,
        ellipsoid: Ellipsoid,
        prime_meridian: PrimeMeridian,
    ) -> Self {
        Self {
            id: SmolStr::new_static(id),
            ellipsoid,
            prime_meridian,
        }
    }

    /// Creates a new [`GeodeticDatum`] from the given parts.
    pub fn new<T: AsRef<str>>(id: T, ellipsoid: Ellipsoid, prime_meridian: PrimeMeridian) -> Self {
        Self {
            id: SmolStr::new(id),
            ellipsoid,
            prime_meridian,
        }
    }

    /// Return this datum's id.
    pub fn id(&self) -> &str {
        self.id.as_str()
    }

    /// Return a reference to this datum's ellipsoid.
    pub fn ellipsoid(&self) -> &Ellipsoid {
        &self.ellipsoid
    }

    /// Return a reference to this datum's prime meridian.
    pub fn prime_meridian(&self) -> &PrimeMeridian {
        &self.prime_meridian
    }

    /// Convert [normalized geodetic coordinates][LLH] into [normalized geocentric coordinates][XYZ].
    /// NOTE: As mentioned in the EPSG guidance note 7 part 2, paragraph 4.1.1, this transformation
    /// first changes the longitude origin from this datum's prime meridian to the Greenwich
    /// meridian.
    pub fn llh_to_xyz(&self, llh: LLH) -> XYZ {
        // change of longitude origin from this datum's prime meridian to the Greenwich meridian.
        let lon = llh.lon + self.prime_meridian.lon();

        if self.ellipsoid.is_spherical() {
            let (sin_lon, cos_lon) = llh.lon.sin_cos();
            let (sin_lat, cos_lat) = llh.lat.sin_cos();
            let d = self.ellipsoid.a() + llh.height;
            let z = d * sin_lat;
            let r = d * cos_lat;
            let x = r * cos_lon;
            let y = r * sin_lon;
            XYZ { x, y, z }
        } else {
            let v = self.ellipsoid.prime_vertical_radius(llh.lat);
            let (sin_lon, cos_lon) = lon.sin_cos();
            let (sin_lat, cos_lat) = llh.lat.sin_cos();

            XYZ {
                x: (v + llh.height) * cos_lat * cos_lon,
                y: (v + llh.height) * cos_lat * sin_lon,
                z: (v * (1.0 - self.ellipsoid.e_sq()) + llh.height) * sin_lat,
            }
        }
    }

    /// Convert [normalized geocentric coordinates][XYZ] into [normalized geodetic coordinates][LLH].
    /// We use the method described in H. Vermeille 2004 - Computing geodetic coordinates from
    /// geocentric coordinates.
    /// NOTE: As mentioned in the EPSG guidance note 7 part 2, paragraph 4.1.1, this transformation
    /// changes the calculated longitude origin (Greenwich meridian) to this Datum's prime
    /// meridian.
    pub fn xyz_to_llh(&self, xyz: XYZ) -> LLH {
        let x = xyz.x.m();
        let y = xyz.y.m();
        let z = xyz.z.m();
        let a = self.ellipsoid.a().m();

        if self.ellipsoid.is_spherical() {
            let lon = y.atan2(x);
            let r = x.hypot(y);
            let lat = z.atan2(r);
            let h = r.hypot(z) - a;
            LLH {
                lon: Lon::new(lon * RAD),
                lat: Lat::new(lat * RAD),
                height: h * M,
            }
        } else {
            let e2 = self.ellipsoid.e_sq();
            let e4 = e2 * e2;

            // Eq (1)
            let _r = x.hypot(y); // sqrt(x^2 + y^2)
            let _p = _r / a;
            let p = _p * _p; // (x^2 + y^2) / a^2

            // Eq (2)
            let z_a = z / a;
            let q = (1. - e2) * z_a * z_a; // (1 - e^2) * z^2 / a^2

            // Eq (3)
            let r = (p + q - e4) / 6.;

            // Eq (4)
            let s = e4 * p * q / (4. * r.powi(3));

            let t = (1. + s + (s * (2. + s)).sqrt()).cbrt();
            let u = r * (1. + t + 1. / t);
            let v = (u * u + e4 * q).sqrt();
            let w = e2 * (u + v - q) / (2. * v);
            let k = (u + v + w * w).sqrt() - w;
            let D = k * _r / (k + e2);

            let _h = D.hypot(z);
            let h = (k + e2 - 1.) / k * _h;
            let lat_rad = 2. * z.atan2(D + _h);

            let lon_rad = if x == 0. && y == 0. {
                0. // on polar axis, we choose lon = 0.
            } else if y >= 0. {
                PI_2 - 2. * x.atan2(_r + y)
            } else {
                -PI_2 + 2. * x.atan2(_r - y)
            };

            let mut lon = Lon::new(lon_rad * RAD);

            // Change longitude origin from Greenwich meridian to this datum's prime meridian
            lon = lon - self.prime_meridian().lon();

            LLH {
                lon,
                lat: Lat::new(lat_rad * RAD),
                height: h * M,
            }
        }
    }
}

impl PartialEq for GeodeticDatum {
    fn eq(&self, other: &Self) -> bool {
        self.id == other.id
            && self.ellipsoid == other.ellipsoid
            && self.prime_meridian == other.prime_meridian
    }
}

/// Well known datum definition.
pub mod consts {
    use super::GeodeticDatum;
    use crate::geodesy::{ellipsoid, prime_meridian};

    pub const WGS84: GeodeticDatum = GeodeticDatum::new_static(
        "WGS84",
        ellipsoid::consts::WGS84,
        prime_meridian::consts::GREENWICH,
    );

    pub const GGRS87: GeodeticDatum = GeodeticDatum::new_static(
        "GGRS87",
        ellipsoid::consts::GRS80,
        prime_meridian::consts::GREENWICH,
    );

    pub const NAD83: GeodeticDatum = GeodeticDatum::new_static(
        "NAD83",
        ellipsoid::consts::GRS80,
        prime_meridian::consts::GREENWICH,
    );

    pub const RNB72: GeodeticDatum = GeodeticDatum::new_static(
        "RNB72",
        ellipsoid::consts::INTL,
        prime_meridian::consts::GREENWICH,
    );
}

#[cfg(test)]
mod tests {
    use smol_str::SmolStr;

    use super::GeodeticDatum;
    use crate::cs::cartesian::{approx_eq_xyz, CartesianErrors, XYZ};
    use crate::cs::geodetic::{approx_eq_llh, GeodeticErrors, Lat, Lon, LLH};
    use crate::geodesy::ellipsoid::consts::WGS84;
    use crate::geodesy::prime_meridian::consts::GREENWICH;
    use crate::geodesy::{ellipsoid, geodetic_datum, prime_meridian, Ellipsoid, PrimeMeridian};
    use crate::units::angle::DEG;
    use crate::units::length::M;

    #[test]
    fn clone() {
        let d = GeodeticDatum::new(
            "WGS 84",
            Ellipsoid::from_ainvf("WGS84", 6_378_137.0 * M, 298.257_223_563),
            PrimeMeridian::new("Greenwich", Lon::new(0.0 * DEG)),
        );
        let cpy = d.clone();
        assert_eq!(d, cpy);
        let _e = d.ellipsoid;
    }

    #[test]
    fn partial_eq() {
        let d = GeodeticDatum::new(
            "WGS 84",
            ellipsoid::consts::WGS84,
            prime_meridian::consts::GREENWICH,
        );
        let same = GeodeticDatum::new(
            "WGS 84",
            ellipsoid::consts::WGS84,
            prime_meridian::consts::GREENWICH,
        );
        let different_id = GeodeticDatum::new(
            "WGS 84.1",
            ellipsoid::consts::WGS84,
            prime_meridian::consts::GREENWICH,
        );
        let different_ellipsoid = GeodeticDatum::new(
            "WGS 84",
            ellipsoid::consts::GRS80,
            prime_meridian::consts::GREENWICH,
        );
        let different_pm = GeodeticDatum::new(
            "WGS 84",
            ellipsoid::consts::WGS84,
            prime_meridian::consts::PARIS,
        );

        assert!(d.eq(&same));
        assert!(!d.ne(&same));

        assert_ne!(d, different_id);
        assert_ne!(d, different_ellipsoid);
        assert_ne!(d, different_pm);
    }

    #[test]
    fn llh_to_xyz() {
        let datum = geodetic_datum::consts::WGS84;

        let computed = datum.llh_to_xyz(LLH {
            lon: Lon::new(2.344 * DEG),
            lat: Lat::new(44.0 * DEG),
            height: 100.0 * M,
        });

        let expected = XYZ {
            x: 4591703.1092 * M,
            y: 187953.8205 * M,
            z: 4408161.0783 * M,
        };
        assert!(approx_eq_xyz(
            &computed,
            &expected,
            &CartesianErrors::default()
        ));
    }

    #[test]
    fn llh_to_xyz_on_sphere() {
        let datum = GeodeticDatum {
            id: SmolStr::new_static("Spherical WGS84"),
            prime_meridian: GREENWICH,
            ellipsoid: Ellipsoid::from_ab("Spherical WGS84", WGS84.a(), WGS84.a()),
        };

        let computed = datum.llh_to_xyz(LLH {
            lon: Lon::new(17.7562015132 * DEG),
            lat: Lat::new(45.3935192042 * DEG),
            height: 133.12 * M,
        });

        let expected = XYZ {
            x: 4265666.7773 * M,
            y: 1365959.7291 * M,
            z: 4540987.8537 * M,
        };
        assert!(approx_eq_xyz(
            &computed,
            &expected,
            &CartesianErrors::default()
        ));
    }

    #[test]
    fn xyz_to_llh() {
        let datum = geodetic_datum::consts::WGS84;

        let computed = datum.xyz_to_llh(XYZ {
            x: 4591703.1092 * M,
            y: 187953.8205 * M,
            z: 4408161.0783 * M,
        });

        let expected = LLH {
            lon: Lon::new(2.344 * DEG),
            lat: Lat::new(44.0 * DEG),
            height: 100. * M,
        };
        assert!(approx_eq_llh(
            &computed,
            &expected,
            &GeodeticErrors::small()
        ));
    }

    #[test]
    fn xyz_to_llh_on_sphere() {
        let datum = GeodeticDatum {
            id: SmolStr::new_static("Spherical WGS84"),
            prime_meridian: GREENWICH,
            ellipsoid: Ellipsoid::from_ab("Spherical WGS84", WGS84.a(), WGS84.a()),
        };

        let computed = datum.xyz_to_llh(XYZ {
            x: 4265666.7773 * M,
            y: 1365959.7291 * M,
            z: 4540987.8537 * M,
        });

        let expected = LLH {
            lon: Lon::new(17.7562015132 * DEG),
            lat: Lat::new(45.3935192042 * DEG),
            height: 133.12 * M,
        };
        assert!(approx_eq_llh(
            &computed,
            &expected,
            &GeodeticErrors::small()
        ));
    }
}
