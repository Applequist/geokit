use super::{ellipsoid::Ellipsoid, prime_meridian::PrimeMeridian};
use crate::crs::CoordSpace;
use crate::id::Id;
use crate::transformation::InvertibleTransformation;

/// A `datum` is the information required to fix a coordinate system to an object.
/// A `GeodeticDatum` is a `datum` describing the relationship of an ellipsoidal model of the Earth
/// with the real Earth.
/// It is defined by an [Ellipsoid] and a [PrimeMeridian].
#[derive(Debug, Clone)]
pub struct GeodeticDatum {
    id: Id,
    ellipsoid: Ellipsoid,
    prime_meridian: PrimeMeridian,
    to_ref: Option<(CoordSpace, Id, Box<dyn InvertibleTransformation>)>,
}

impl GeodeticDatum {
    /// Creates a new [`GeodeticDatum`].
    pub fn new(
        id: Id,
        ellipsoid: Ellipsoid,
        prime_meridian: PrimeMeridian,
        to_ref: Option<(CoordSpace, Id, Box<dyn InvertibleTransformation>)>,
    ) -> Self {
        Self {
            id,
            ellipsoid,
            prime_meridian,
            to_ref,
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

impl Default for GeodeticDatum {
    /// Returns the WGS 84 datum (epsg:6326) as default GeodeticDatum.
    fn default() -> Self {
        GeodeticDatum::new(
            Id::name("WGS 84"),
            Ellipsoid::default(),
            PrimeMeridian::default(),
            None,
        )
    }
}

#[cfg(test)]
mod tests {

    use super::GeodeticDatum;
    use crate::geodesy::{Ellipsoid, PrimeMeridian};
    use crate::id::Id;

    #[test]
    fn clone() {
        let d = GeodeticDatum::new(
            Id::name("WGS 84"),
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
            Id::name("WGS 84"),
            Ellipsoid::from_ainvf(6_378_137.0, 298.257_223_563),
            PrimeMeridian::new(0.0),
            None,
        );
        let same = GeodeticDatum::new(
            Id::name("WGS 84"),
            Ellipsoid::from_ainvf(6_378_137.0, 298.257_223_563),
            PrimeMeridian::new(0.0),
            None,
        );
        let different_id = GeodeticDatum::new(
            Id::name("WGS 84.1"),
            Ellipsoid::from_ainvf(6_378_137.0, 298.257_223_563),
            PrimeMeridian::new(0.0),
            None,
        );
        let different_ellipsoid = GeodeticDatum::new(
            Id::name("WGS 84"),
            Ellipsoid::from_ab(6_378_137.0, 6_378_137.0),
            PrimeMeridian::new(0.0),
            None,
        );
        let different_pm = GeodeticDatum::new(
            Id::name("WGS 84"),
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
        assert_eq!(wgs84.id, Id::name("WGS 84"));
        assert_eq!(wgs84.ellipsoid, Ellipsoid::default());
        assert_eq!(wgs84.prime_meridian, PrimeMeridian::default());
    }
}
