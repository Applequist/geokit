use smol_str::SmolStr;

use crate::geodesy::GeodeticDatum;

use super::Crs;

/// A [GeocentricCrs] is a **3D cartesian coordinates reference system** in which
/// coordinates are given by distance **in meters** along the following axes:
/// - X: axis from the center of the datum's ellipsoid in the equatorial and prime meridian plane,
/// - Y: axis from the center of the datum's ellipsoid in the equatorial plane and 90 degrees
/// meridian plane (east)
/// - Z: axis from the center of the datum's ellipsoid through the north pole.
#[derive(Debug, Clone, PartialEq)]
pub struct GeocentricCrs {
    id: SmolStr,
    datum: GeodeticDatum,
}

impl GeocentricCrs {
    /// Creates a new [`GeocentricCrs`].
    pub fn new<T: AsRef<str>>(id: T, datum: GeodeticDatum) -> Self {
        debug_assert!(
            datum.prime_meridian().lon() == 0.0,
            "Expected Greenwich prime meridian"
        );
        Self {
            id: SmolStr::new(id),
            datum,
        }
    }

    pub fn id(&self) -> &str {
        self.id.as_str()
    }

    #[inline]
    pub fn dim(&self) -> usize {
        3
    }

    /// Return this CRS geodetic datum as a reference.
    #[inline]
    pub fn datum(&self) -> &GeodeticDatum {
        &self.datum
    }
}

impl Crs for GeocentricCrs {
    fn is_normalized(&self) -> bool {
        true
    }
}

#[cfg(test)]
mod tests {

    use super::*;
    use crate::geodesy::*;

    #[test]
    fn clone() {
        let geoc = GeocentricCrs::new(
            "WGS84",
            GeodeticDatum::new(
                "WGS84",
                ellipsoid::consts::WGS84,
                prime_meridian::consts::GREENWICH,
                None,
            ),
        );
        let cpy = geoc.clone();
        assert_eq!(geoc, cpy);
    }

    #[test]
    fn partial_eq() {
        let geoc = GeocentricCrs::new(
            "WGS84",
            GeodeticDatum::new(
                "WGS84",
                ellipsoid::consts::WGS84,
                prime_meridian::consts::GREENWICH,
                None,
            ),
        );
        let cpy = geoc.clone();
        assert!(geoc.eq(&cpy));
        assert!(!geoc.ne(&cpy));

        let mut different_tag = geoc.clone();
        different_tag.id = "WGS 84.1".into();
        assert_ne!(geoc, different_tag);

        let mut different_datum = geoc.clone();
        different_datum.datum = GeodeticDatum::new(
            "WGS 84.1",
            ellipsoid::consts::GRS80,
            prime_meridian::consts::GREENWICH,
            None,
        );
        assert_ne!(geoc, different_datum);
    }
}
