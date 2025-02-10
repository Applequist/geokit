use super::{Ellipsoid, PrimeMeridian};
use crate::cs::cartesian::XYZ;
use crate::cs::geodetic::{Lat, Lon, LLH};
use crate::units::length::M;
use derive_more::derive::Display;
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
        let lon = llh.lon + self.prime_meridian.lon().angle();

        let v = self.ellipsoid.prime_vertical_radius(llh.lat);
        let (sin_lon, cos_lon) = lon.sin_cos();
        let (sin_lat, cos_lat) = llh.lat.sin_cos();

        XYZ {
            x: (v + llh.height) * cos_lat * cos_lon,
            y: (v + llh.height) * cos_lat * sin_lon,
            z: (v * (1.0 - self.ellipsoid.e_sq()) + llh.height) * sin_lat,
        }
    }

    /// Convert [normalized geocentric coordinates][XYZ] into [normalized geodetic coordinates][LLH].
    /// NOTE: As mentioned in the EPSG guidance note 7 part 2, paragraph 4.1.1, this transformation
    /// changes the calculated longitude origin (Greenwich meridian) to this Datum's prime
    /// meridian.
    pub fn xyz_to_llh(&self, xyz: XYZ) -> LLH {
        let x = xyz.x;
        let y = xyz.y;
        let z = xyz.z;

        let a2 = self.ellipsoid.a_sq();
        let b2 = self.ellipsoid.b_sq();
        let e2 = self.ellipsoid.e_sq();

        let greenwich_lon = Lon::new(y.atan2(x));
        // Change longitude origin from Greenwich meridian to this datum's prime meridian
        let lon = greenwich_lon - self.prime_meridian().lon();

        let p = x.hypot(y);
        let mut lat = z.atan2(p * (1.0 - e2));
        let (sin_lat, cos_lat) = lat.sin_cos();
        let n = a2 / (a2 * cos_lat * cos_lat + b2 * sin_lat * sin_lat).sqrt() * M;

        let mut h = p / cos_lat - n;
        loop {
            let next_lat = z.atan2(p * (1.0 - e2 * n / (n + h)));
            let (sin_nlat, cos_nlat) = next_lat.sin_cos();
            let next_n = a2 / ((a2 * cos_nlat * cos_nlat) + b2 * sin_nlat * sin_nlat).sqrt() * M;
            let next_h = p / cos_nlat - next_n;
            let delta_lat = (lat - next_lat).abs();
            let delta_h = (h - next_h).abs();
            lat = next_lat;
            h = next_h;
            if delta_lat.rad() < 0.5e-5 && delta_h.m() < 0.5e-3 {
                break;
            }
        }

        LLH {
            lon,
            lat: Lat::new(lat),
            height: h,
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
    use approx::assert_abs_diff_eq;
    use regex::Regex;

    use crate::cs::cartesian::{CartesianErrors, XYZ};
    use crate::cs::geodetic::{GeodeticErrors, Lat, Lon, LLH};
    use crate::geodesy::{ellipsoid, geodetic_datum, prime_meridian, Ellipsoid, PrimeMeridian};
    use crate::math::Float;
    use crate::units::angle::DEG;
    use crate::units::length::M;

    use super::consts::WGS84;
    use super::GeodeticDatum;

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

        let xyz = datum.llh_to_xyz(LLH {
            lon: Lon::new(2.344 * DEG),
            lat: Lat::new(44.0 * DEG),
            height: 100.0 * M,
        });

        xyz.approx_eq(
            &XYZ {
                x: 4591703.1092 * M,
                y: 187953.8205 * M,
                z: 4408161.0783 * M,
            },
            CartesianErrors::default(),
        );
    }

    #[test]
    fn xyz_to_llh() {
        let datum = geodetic_datum::consts::WGS84;

        let llh = datum.xyz_to_llh(XYZ {
            x: 4591703.1092 * M,
            y: 187953.8205 * M,
            z: 4408161.0783 * M,
        });

        llh.approx_eq(
            &LLH {
                lon: Lon::new(2.344 * DEG),
                lat: Lat::new(44.0 * DEG),
                height: 100. * M,
            },
            GeodeticErrors::default(),
        );
    }

    const XYZLLH: &'static str = include_str!("xyz_to_llh_100.txt");

    fn read_conversions() -> impl Iterator<Item = (XYZ, LLH)> {
        XYZLLH
            .lines()
            .filter(|&l| !l.is_empty() && !l.starts_with("#"))
            .map(move |l| {
                let xyzllh = l
                    .split(',')
                    .map(|s| s.trim().parse::<Float>().unwrap())
                    .collect::<Vec<_>>();
                let xyz = XYZ {
                    x: xyzllh[0] * M,
                    y: xyzllh[1] * M,
                    z: xyzllh[2] * M,
                };
                let llh = LLH {
                    lon: Lon::new(xyzllh[3] * DEG),
                    lat: Lat::new(xyzllh[4] * DEG),
                    height: xyzllh[5] * M,
                };
                (xyz, llh)
            })
    }

    #[test]
    fn xyz_to_lly_100() {
        let datum = WGS84;

        for (xyz, expected_llh) in read_conversions() {
            let llh = datum.xyz_to_llh(xyz);
            llh.approx_eq(&expected_llh, GeodeticErrors::default());
        }
    }
}
