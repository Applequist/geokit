use std::fmt::Debug;

use smol_str::SmolStr;

use super::{Ellipsoid, PrimeMeridian};
use crate::operation::{
    identity,
    transformation::{GeocentricTranslation, Helmert7Params, RotationConvention},
    Operation,
};

/// Coordinates can be transformed between different datum.
/// A [`ToWGS84`] specifies the transformation to transform normalized geocentric coordinates from
/// a datum to geocentric WGS84 coordinates.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ToWGS84 {
    /// A simple translation of source normalized geocentric coordinates by adding an offset **in
    /// meters**.
    GeocentricTranslation { tx: f64, ty: f64, tz: f64 },
    /// The Helmert 7-parameters transformation transforms normalized geocentric coordinates using
    /// small rotations around x, y and z axes, a translation and a small scaling in ppm.
    /// Rotations angle are in radians, translation along axes are in meters and scale is in ppm.
    Helmert7Params {
        conv: RotationConvention,
        rx: f64,
        ry: f64,
        rz: f64,
        tx: f64,
        ty: f64,
        tz: f64,
        scale: f64,
    },
}

impl ToWGS84 {
    /// Returns the actual transformation.
    pub fn transformation(&self) -> Box<dyn Operation> {
        match *self {
            ToWGS84::GeocentricTranslation { tx, ty, tz } => {
                GeocentricTranslation::new(tx, ty, tz).boxed()
            }
            ToWGS84::Helmert7Params {
                conv,
                rx,
                ry,
                rz,
                tx,
                ty,
                tz,
                scale,
            } => Helmert7Params::new(conv, rx, ry, rz, tx, ty, tz, scale).boxed(),
        }
    }
}

/// A `datum` is the information required to fix a coordinate system to an object.
/// A `GeodeticDatum` is a `datum` describing the relationship of an ellipsoidal model of the Earth
/// with the real Earth.
/// It is defined by an [Ellipsoid] and a [PrimeMeridian].
#[derive(Debug, Clone)]
pub struct GeodeticDatum {
    id: SmolStr,
    ellipsoid: Ellipsoid,
    prime_meridian: PrimeMeridian,
    to_wgs84: Option<ToWGS84>,
}

impl GeodeticDatum {
    /// Creates a new [`GeodeticDatum`].
    pub const fn new_static(
        id: &'static str,
        ellipsoid: Ellipsoid,
        prime_meridian: PrimeMeridian,
        to_wgs84: Option<ToWGS84>,
    ) -> Self {
        Self {
            id: SmolStr::new_static(id),
            ellipsoid,
            prime_meridian,
            to_wgs84,
        }
    }

    /// Creates a new [`GeodeticDatum`].
    pub fn new<T: AsRef<str>>(
        id: T,
        ellipsoid: Ellipsoid,
        prime_meridian: PrimeMeridian,
        to_wgs84: Option<ToWGS84>,
    ) -> Self {
        Self {
            id: SmolStr::new(id),
            ellipsoid,
            prime_meridian,
            to_wgs84,
        }
    }

    /// Return this datum's id as a reference.
    pub fn id(&self) -> &str {
        self.id.as_str()
    }

    /// Return a copy of this datum's ellipsoid.
    pub fn ellipsoid(&self) -> &Ellipsoid {
        &self.ellipsoid
    }

    /// Return a copy of this datum's prime meridian.
    pub fn prime_meridian(&self) -> &PrimeMeridian {
        &self.prime_meridian
    }

    /// Returns the id of the reference datum used by this datum.
    pub fn ref_datum_id(&self) -> &str {
        self.to_wgs84.map_or(self.id(), |_| "WGS84")
    }

    pub fn to_wgs84(&self) -> Box<dyn Operation> {
        self.to_wgs84
            .map_or(identity::<3>().boxed(), |spec| spec.transformation())
    }
}

impl PartialEq for GeodeticDatum {
    fn eq(&self, other: &Self) -> bool {
        self.id == other.id
            && self.ellipsoid == other.ellipsoid
            && self.prime_meridian == other.prime_meridian
    }
}

pub mod consts {

    use crate::{
        geodesy::{ellipsoid, prime_meridian},
        operation::transformation::RotationConvention,
    };

    use super::{GeodeticDatum, ToWGS84};

    pub const WGS84: GeodeticDatum = GeodeticDatum::new_static(
        "WGS84",
        ellipsoid::consts::WGS84,
        prime_meridian::consts::GREENWICH,
        None,
    );

    pub const GGRS87: GeodeticDatum = GeodeticDatum::new_static(
        "GGRS87",
        ellipsoid::consts::GRS80,
        prime_meridian::consts::GREENWICH,
        Some(ToWGS84::GeocentricTranslation {
            tx: -199.87,
            ty: 74.79,
            tz: 246.64,
        }),
    );

    pub const NAD83: GeodeticDatum = GeodeticDatum::new_static(
        "NAD83",
        ellipsoid::consts::GRS80,
        prime_meridian::consts::GREENWICH,
        None,
    );

    pub const RNB72: GeodeticDatum = GeodeticDatum::new_static(
        "RNB72",
        ellipsoid::consts::INTL,
        prime_meridian::consts::GREENWICH,
        Some(ToWGS84::Helmert7Params {
            conv: RotationConvention::CoordinateFrame,
            // TODO: Convert arcsec to radians
            rx: -0.33657,
            ry: 0.456955,
            rz: -1.84218,
            tx: 106.869,
            ty: -52.2978,
            tz: 103.724,
            scale: 0.,
        }),
    );
}

#[cfg(test)]
mod tests {

    use super::GeodeticDatum;
    use crate::geodesy::{ellipsoid, geodetic_datum, prime_meridian, Ellipsoid, PrimeMeridian};

    #[test]
    fn clone() {
        let d = GeodeticDatum::new(
            "WGS 84",
            Ellipsoid::from_ainvf("WGS84", 6_378_137.0, 298.257_223_563),
            PrimeMeridian::new("Greenwich", 0.0),
            None,
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
            None,
        );
        let same = GeodeticDatum::new(
            "WGS 84",
            ellipsoid::consts::WGS84,
            prime_meridian::consts::GREENWICH,
            None,
        );
        let different_id = GeodeticDatum::new(
            "WGS 84.1",
            ellipsoid::consts::WGS84,
            prime_meridian::consts::GREENWICH,
            None,
        );
        let different_ellipsoid = GeodeticDatum::new(
            "WGS 84",
            ellipsoid::consts::GRS80,
            prime_meridian::consts::GREENWICH,
            None,
        );
        let different_pm = GeodeticDatum::new(
            "WGS 84",
            ellipsoid::consts::WGS84,
            prime_meridian::consts::PARIS,
            None,
        );

        assert!(d.eq(&same));
        assert!(!d.ne(&same));

        assert_ne!(d, different_id);
        assert_ne!(d, different_ellipsoid);
        assert_ne!(d, different_pm);
    }

    #[test]
    fn ref_datum_id() {
        assert_eq!(geodetic_datum::consts::WGS84.ref_datum_id(), "WGS84");
        assert_eq!(geodetic_datum::consts::RNB72.ref_datum_id(), "WGS84");
    }
}
