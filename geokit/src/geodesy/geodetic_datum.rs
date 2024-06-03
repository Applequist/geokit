use std::fmt::Debug;

use smol_str::SmolStr;

use crate::operation::{
    identity,
    transformation::{GeocentricTranslation, Helmert7Params, RotationConvention},
    Operation,
};
use crate::quantity::scale::PPM;

use super::{Ellipsoid, PrimeMeridian};

/// Coordinates can be transformed between different datum.
/// A [`DatumTransformation`] specifies how to transform normalized geocentric coordinates from
/// a datum to geocentric coordinates in a reference datum, usually WGS84.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum DatumTransformation {
    /// A simple translation of source normalized geocentric coordinates by adding an offset **in
    /// meters**.
    GeocentricTranslation { tx: f64, ty: f64, tz: f64 },
    /// The Helmert 7-parameters transformation transforms normalized geocentric coordinates using
    /// small rotations around x, y and z axes, a translation and a small scaling in ppm.
    /// Rotations angle are in radians, translation along axes are in meters and scale is in ppm.
    Helmert7Params {
        conv: RotationConvention,
        rotation: [f64; 3],
        translation: [f64; 3],
        scale: PPM,
    },
}

impl DatumTransformation {
    /// Returns the actual transformation.
    pub fn transformation(&self) -> Box<dyn Operation> {
        match *self {
            DatumTransformation::GeocentricTranslation { tx, ty, tz } => {
                GeocentricTranslation::new(tx, ty, tz).boxed()
            }
            DatumTransformation::Helmert7Params {
                conv,
                rotation,
                translation,
                scale,
            } => Helmert7Params::new(conv, rotation, translation, scale).boxed(),
        }
    }
}

/// A `datum` is the information required to fix a coordinate system to an object.
/// A [GeodeticDatum] is a `datum` describing the relationship of an ellipsoidal model of the Earth
/// with the real Earth.
/// It is defined by an id, an [Ellipsoid] and a [PrimeMeridian].
/// It can also come with an optional transformation to a reference datum, usually `WGS84`.
///
/// # Equality
///
/// 2 datums are considered equal if they have **the same id**, the same ellipsoid and prime meridian.
/// That is: it is not enough to use the same ellipsoid and prime meridian for 2 datums to be considered
/// the same.
#[derive(Debug, Clone)]
pub struct GeodeticDatum {
    id: SmolStr,
    ellipsoid: Ellipsoid,
    prime_meridian: PrimeMeridian,
    ref_datum: Option<(SmolStr, DatumTransformation)>,
}

impl GeodeticDatum {
    /// Creates a new [`GeodeticDatum`].
    /// Used to define well known datum as constant.
    pub const fn new_static(
        id: &'static str,
        ellipsoid: Ellipsoid,
        prime_meridian: PrimeMeridian,
        ref_datum: Option<(&'static str, DatumTransformation)>,
    ) -> Self {
        Self {
            id: SmolStr::new_static(id),
            ellipsoid,
            prime_meridian,
            ref_datum: match ref_datum {
                Some((id, tx)) => Some((SmolStr::new_static(id), tx)),
                _ => None,
            },
        }
    }

    /// Creates a new [`GeodeticDatum`] from the given parts.
    pub fn new<T: AsRef<str>>(
        id: T,
        ellipsoid: Ellipsoid,
        prime_meridian: PrimeMeridian,
        ref_datum: Option<(T, DatumTransformation)>,
    ) -> Self {
        Self {
            id: SmolStr::new(id),
            ellipsoid,
            prime_meridian,
            ref_datum: ref_datum.map(|(id, tx)| (SmolStr::new(id), tx)),
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

    /// Returns the id of the reference datum used by this datum.
    pub fn ref_datum_id(&self) -> &str {
        match &self.ref_datum {
            Some((id, _)) => id,
            None => self.id(),
        }
    }

    /// Returns the [Operation] to transform **geocentric** coordinates in this datum into
    /// **geocentric** coordinates in the reference datum.
    pub fn to_ref_datum(&self) -> Box<dyn Operation> {
        match &self.ref_datum {
            Some((_, tx)) => tx.transformation(),
            None => identity(3, 3).boxed(),
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
    use crate::quantity::scale::PPM;
    use crate::{
        geodesy::{ellipsoid, prime_meridian},
        operation::transformation::RotationConvention,
    };

    use super::{DatumTransformation, GeodeticDatum};

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
        Some((
            "WGS84",
            DatumTransformation::GeocentricTranslation {
                tx: -199.87,
                ty: 74.79,
                tz: 246.64,
            },
        )),
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
        Some((
            "WGS84",
            DatumTransformation::Helmert7Params {
                conv: RotationConvention::CoordinateFrame,
                rotation: [-1.63172286E-6, 2.21538036E-6, -8.9311407E-6],
                translation: [106.869, -52.2978, 103.724],
                scale: PPM(0.),
            },
        )),
    );
}

#[cfg(test)]
mod tests {
    use crate::geodesy::{ellipsoid, geodetic_datum, prime_meridian, Ellipsoid, PrimeMeridian};
    use crate::quantity::angle::units::DEG;

    use super::GeodeticDatum;

    #[test]
    fn clone() {
        let d = GeodeticDatum::new(
            "WGS 84",
            Ellipsoid::from_ainvf("WGS84", 6_378_137.0, 298.257_223_563),
            PrimeMeridian::new("Greenwich", 0.0 * DEG),
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
