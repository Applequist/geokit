use std::fmt::Debug;

use super::{Ellipsoid, PrimeMeridian};
use crate::operation::transformation::RotationConvention;
use crate::tag::Tag;

#[derive(Debug, Clone)]
pub struct ReferenceDatum {
    tag: Tag,
    to_ref: DatumTransformation,
}

/// Coordinates can be transformed between different datum.
/// A [`DatumTransformation`] specifies the supported transformations.
#[derive(Debug, Clone, Copy)]
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
    tag: Tag,
    ellipsoid: Ellipsoid,
    prime_meridian: PrimeMeridian,
    ref_datum: Option<ReferenceDatum>,
}

impl GeodeticDatum {
    /// Creates a new [`GeodeticDatum`].
    pub fn new<T: Into<Tag>>(
        tag: T,
        ellipsoid: Ellipsoid,
        prime_meridian: PrimeMeridian,
        ref_datum: Option<ReferenceDatum>,
    ) -> Self {
        Self {
            tag: tag.into(),
            ellipsoid,
            prime_meridian,
            ref_datum,
        }
    }

    /// Return this datum's id as a reference.
    pub fn tag(&self) -> &Tag {
        &self.tag
    }

    /// Return a copy of this datum's ellipsoid.
    pub fn ellipsoid(&self) -> Ellipsoid {
        self.ellipsoid
    }

    /// Return a copy of this datum's prime meridian.
    pub fn prime_meridian(&self) -> PrimeMeridian {
        self.prime_meridian
    }

    /// Return the option [DatumTransformation] to a reference datum.
    pub fn ref_datum(&self) -> Option<&ReferenceDatum> {
        self.ref_datum.as_ref()
    }
}

impl PartialEq for GeodeticDatum {
    fn eq(&self, other: &Self) -> bool {
        self.tag == other.tag
            && self.ellipsoid == other.ellipsoid
            && self.prime_meridian == other.prime_meridian
    }
}

impl Default for GeodeticDatum {
    /// Return the WGS 84 datum (EPSG:6326) as default GeodeticDatum.
    fn default() -> Self {
        GeodeticDatum::new(
            Tag::full("WGS 1984", "EPSG", 6326),
            Ellipsoid::default(),
            PrimeMeridian::default(),
            None,
        )
    }
}

impl Debug for GeodeticDatum {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("GeodeticDatum")
            .field("tag", &self.tag)
            .field("ellipsoid", &self.ellipsoid)
            .field("prime_meridian", &self.prime_meridian)
            .field("ref_datum", &self.ref_datum)
            .finish()
    }
}

#[cfg(test)]
mod tests {

    use super::GeodeticDatum;
    use crate::geodesy::{Ellipsoid, PrimeMeridian};
    use crate::tag::Tag;

    #[test]
    fn clone() {
        let d = GeodeticDatum::new(
            Tag::name("WGS 84"),
            Ellipsoid::from_ainvf(6_378_137.0, 298.257_223_563),
            PrimeMeridian::new(0.0),
            None,
        );
        let cpy = d.clone();
        assert_eq!(d, cpy);
        let _e = d.ellipsoid;
    }

    #[test]
    fn partial_eq() {
        let d = GeodeticDatum::new(
            Tag::name("WGS 84"),
            Ellipsoid::from_ainvf(6_378_137.0, 298.257_223_563),
            PrimeMeridian::new(0.0),
            None,
        );
        let same = GeodeticDatum::new(
            Tag::name("WGS 84"),
            Ellipsoid::from_ainvf(6_378_137.0, 298.257_223_563),
            PrimeMeridian::new(0.0),
            None,
        );
        let different_id = GeodeticDatum::new(
            Tag::name("WGS 84.1"),
            Ellipsoid::from_ainvf(6_378_137.0, 298.257_223_563),
            PrimeMeridian::new(0.0),
            None,
        );
        let different_ellipsoid = GeodeticDatum::new(
            Tag::name("WGS 84"),
            Ellipsoid::from_ab(6_378_137.0, 6_378_137.0),
            PrimeMeridian::new(0.0),
            None,
        );
        let different_pm = GeodeticDatum::new(
            Tag::name("WGS 84"),
            Ellipsoid::from_ainvf(6_378_137.0, 298.257_223_563),
            PrimeMeridian::new(2.32_f64.to_radians()),
            None,
        );

        assert!(d.eq(&same));
        assert!(!d.ne(&same));

        assert_ne!(d, different_id);
        assert_ne!(d, different_ellipsoid);
        assert_ne!(d, different_pm);
    }

    #[test]
    fn default() {
        let wgs84 = GeodeticDatum::default();
        assert_eq!(
            wgs84.tag,
            Tag::full("WGS_1984", "EPSG", 6326),
            "Expected 'WGS 84'. Got {}",
            wgs84.tag
        );
        assert_eq!(wgs84.ellipsoid, Ellipsoid::default());
        assert_eq!(wgs84.prime_meridian, PrimeMeridian::default());
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
