use std::fmt::Debug;

use smol_str::SmolStr;

use super::{Ellipsoid, PrimeMeridian};
use crate::operation::transformation::RotationConvention;

#[derive(Debug, Clone)]
pub struct ReferenceDatum {
    id: SmolStr,
    to_ref: DatumTransformation,
}

/// Coordinates can be transformed between different datum.
/// A [`DatumTransformation`] specifies the supported transformations.
#[derive(Debug, Clone)]
pub enum DatumTransformation {
    /// A simple datum shift transforms source normalized geocentric coordinates by adding an offset **in
    /// meters**.
    GeocentricTranslation { tx: f64, ty: f64, tz: f64 },
    /// The Helmert 7-parameters transformation transforms normalized geocentric coordinates using
    /// small rotations around x, y and z axes, a translation and a small scaling in ppm.
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

/// A `datum` is the information required to fix a coordinate system to an object.
/// A `GeodeticDatum` is a `datum` describing the relationship of an ellipsoidal model of the Earth
/// with the real Earth.
/// It is defined by an [Ellipsoid] and a [PrimeMeridian].
#[derive(Clone)]
pub struct GeodeticDatum {
    id: SmolStr,
    ellipsoid: Ellipsoid,
    prime_meridian: PrimeMeridian,
    ref_datum: Option<ReferenceDatum>,
}

impl GeodeticDatum {
    /// Creates a new [`GeodeticDatum`].
    pub const fn new_static(
        id: &'static str,
        ellipsoid: Ellipsoid,
        prime_meridian: PrimeMeridian,
        ref_datum: Option<ReferenceDatum>,
    ) -> Self {
        Self {
            id: SmolStr::new_static(id),
            ellipsoid,
            prime_meridian,
            ref_datum,
        }
    }

    /// Creates a new [`GeodeticDatum`].
    pub fn new<T: AsRef<str>>(
        id: T,
        ellipsoid: Ellipsoid,
        prime_meridian: PrimeMeridian,
        ref_datum: Option<ReferenceDatum>,
    ) -> Self {
        Self {
            id: SmolStr::new(id),
            ellipsoid,
            prime_meridian,
            ref_datum,
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

    /// Return the option [DatumTransformation] to a reference datum.
    pub fn ref_datum(&self) -> Option<&ReferenceDatum> {
        self.ref_datum.as_ref()
    }
}

impl PartialEq for GeodeticDatum {
    fn eq(&self, other: &Self) -> bool {
        self.id == other.id
            && self.ellipsoid == other.ellipsoid
            && self.prime_meridian == other.prime_meridian
    }
}

impl Debug for GeodeticDatum {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("GeodeticDatum")
            .field("id", &self.id)
            .field("ellipsoid", &self.ellipsoid)
            .field("prime_meridian", &self.prime_meridian)
            .field("ref_datum", &self.ref_datum)
            .finish()
    }
}

pub mod consts {
    use crate::geodesy::{ellipsoid, prime_meridian};

    use super::GeodeticDatum;

    pub const WGS84: GeodeticDatum = GeodeticDatum::new_static(
        "EPSG:6326",
        ellipsoid::consts::WGS84,
        prime_meridian::consts::GREENWICH,
        None,
    );
}

#[cfg(test)]
mod tests {

    use super::GeodeticDatum;
    use crate::geodesy::{ellipsoid, prime_meridian, Ellipsoid, PrimeMeridian};

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

    #[ignore = "unimplemented"]
    #[test]
    fn debug() {
        unimplemented!()
    }

    #[ignore = "unimplemented"]
    #[test]
    fn display() {
        unimplemented!()
    }
}
