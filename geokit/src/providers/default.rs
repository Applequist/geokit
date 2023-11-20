use crate::{
    geodesy::{Ellipsoid, GeodeticDatum, PrimeMeridian},
    id::Id,
};

use super::GeodesyProvider;

/// The default [`GeodesyProvider`] provides commonly used geodesy elements.
pub struct DefaultGeodesyProvider;

impl GeodesyProvider for DefaultGeodesyProvider {
    fn ellipsoid_ids(&self) -> Vec<Id> {
        vec![Id::name("WGS 84")]
    }

    fn ellipsoid<I: Into<Id>>(&self, id: I) -> Option<Ellipsoid> {
        if id.into() == Id::name("WGS 84") {
            Some(Ellipsoid::default())
        } else {
            None
        }
    }

    fn prime_meridian_ids(&self) -> Vec<Id> {
        vec![Id::name("Greenwich")]
    }

    fn prime_meridian<I: Into<Id>>(&self, id: I) -> Option<PrimeMeridian> {
        if id.into() == Id::name("Greenwich") {
            Some(PrimeMeridian::default())
        } else {
            None
        }
    }

    fn datum_ids(&self) -> Vec<Id> {
        vec![Id::name("WGS_1984")]
    }

    fn datum<I: Into<Id>>(&self, id: I) -> Option<GeodeticDatum> {
        if id.into() == Id::name("WGS_1984") {
            Some(GeodeticDatum::default())
        } else {
            None
        }
    }
}

#[cfg(test)]
mod tests {
    use super::DefaultGeodesyProvider;
    use super::GeodesyProvider;

    #[test]
    fn ellipsoid() {
        let provider = DefaultGeodesyProvider;
        assert_eq!(provider.ellipsoid_ids().len(), 1);
        let wgs84 = provider.ellipsoid("WGS 84");
        assert!(wgs84.is_some());
        assert!((wgs84.unwrap().a - 6378137.0).abs() < 1e-5);
        assert!((wgs84.unwrap().invf - 298.257223563).abs() < 1e-5);
    }

    #[test]
    fn prime_meridian() {
        let provider = DefaultGeodesyProvider;
        assert_eq!(provider.prime_meridian_ids().len(), 1);
        let gw = provider.prime_meridian("Greenwich");
        assert!(gw.is_some());
        assert_eq!(gw.unwrap().greenwich_lon(), 0.0);
    }
}
