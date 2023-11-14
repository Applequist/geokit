use super::{ellipsoid::Ellipsoid, prime_meridian::PrimeMeridian};

/// A `datum` is the information required to fix a coordinate system to an object.
/// A `GeodeticDatum` is a `datum` describing the relationship of an ellipsoidal model of the Earth
/// with the real Earth.
/// It is defined by an [Ellipsoid] and a [PrimeMeridian].
#[derive(Debug, Clone, PartialEq)]
pub struct GeodeticDatum {
    id: String,
    ellipsoid: Ellipsoid,
    prime_meridian: PrimeMeridian,
}

impl GeodeticDatum {
    /// Creates a new [`GeodeticDatum`].
    pub fn new<S: Into<String>>(
        id: S,
        ellipsoid: Ellipsoid,
        prime_meridian: PrimeMeridian,
    ) -> Self {
        Self {
            id: id.into(),
            ellipsoid,
            prime_meridian,
        }
    }
}

impl Default for GeodeticDatum {
    /// Returns the WGS 84 datum (epsg:6326) as default GeodeticDatum.
    fn default() -> Self {
        GeodeticDatum::new("WGS 84", Ellipsoid::default(), PrimeMeridian::default())
    }
}

#[cfg(test)]
mod tests {
    use crate::geodesy::{Ellipsoid, PrimeMeridian};

    use super::GeodeticDatum;

    #[test]
    fn clone() {
        let d = GeodeticDatum::new(
            "WGS84",
            Ellipsoid::from_ainvf(6_378_137.0, 298.257_223_563),
            PrimeMeridian::new(0.0),
        );
        let cpy = d.clone();
        assert_eq!(d, cpy);
        let _s = d.ellipsoid;
    }

    #[test]
    fn partial_eq() {
        let d = GeodeticDatum::new(
            "WGS84",
            Ellipsoid::from_ainvf(6_378_137.0, 298.257_223_563),
            PrimeMeridian::new(0.0),
        );
        let same = GeodeticDatum::new(
            "WGS84",
            Ellipsoid::from_ainvf(6_378_137.0, 298.257_223_563),
            PrimeMeridian::new(0.0),
        );
        let different_id = GeodeticDatum::new(
            "blah",
            Ellipsoid::from_ainvf(6_378_137.0, 298.257_223_563),
            PrimeMeridian::new(0.0),
        );
        let different_ellipsoid = GeodeticDatum::new(
            "WGS84",
            Ellipsoid::from_ab(6_378_137.0, 6_378_137.0),
            PrimeMeridian::new(0.0),
        );
        let different_pm = GeodeticDatum::new(
            "WGS84",
            Ellipsoid::from_ainvf(6_378_137.0, 298.257_223_563),
            PrimeMeridian::new(2.32_f64.to_radians()),
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
        assert_eq!(wgs84.id, String::from("WGS 84"));
        assert_eq!(wgs84.ellipsoid, Ellipsoid::default());
        assert_eq!(wgs84.prime_meridian, PrimeMeridian::default());
    }
}
