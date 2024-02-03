use std::iter::empty;

use phf::phf_ordered_map;

use crate::geodesy::{Ellipsoid, GeodeticDatum, PrimeMeridian};

use super::GeodesyProvider;

macro_rules! ellipsoid_def {
    ( (a = $a:literal, b = $b:literal ) ) => {
        Ellipsoid::new(
            $a,
            $b,
            if $a != $b {
                $a / ($a - $b)
            } else {
                f64::INFINITY
            },
        )
    };
    ( (a = $a:literal, invf = $invf:literal ) ) => {
        Ellipsoid::new($a, $a * (1. - 1. / $invf), $invf)
    };
}

macro_rules! ellipsoids {
    ( $( $name:ident : $id:literal => $def:tt ),+ ) => {

        /// Ids of ellipsoid supported by the default provider.
        pub mod ellipsoids {
            $(pub const $name: &'static str = $id;)+
        }

        static ELLIPSOIDS: phf::OrderedMap<&str, Ellipsoid> = phf_ordered_map! {
            $( $id => ellipsoid_def!($def),)+
        };
    }
}

ellipsoids! {
    AIRY : "Airy 1830" =>  (a = 6_377_563.396, invf = 299.3249646),
    ANDRAE : "Andrae 1876 (Den., Iclnd.)" => (a = 6_377_104.43, invf = 300.0),
    DANISH : "Andrae 1876 (Denmark, Iceland)" => (a = 6_377_019.256_3, invf = 300.0),
    APL_4_9 : "Appl. Physics. 1965" => (a = 6_378_137., invf = 298.25),
    AUST_SA : "Australian Natl & S. Amer. 1969" => (a = 6_378_160., invf = 298.25),
    BESSEL : "Bessel 1841" => (a = 6_377_397.155, invf = 299.1528128),
    BESS_NAM : "Bessel 1841 (Namibia)" => (a = 6_377_483.865, invf = 299.1528128),
    CLRK66 : "Clarke 1866" => (a = 6_378_206.4, b = 6_356_583.8),
    CLRK80 : "Clarke 1880 mod." => (a = 6_378_249.145, invf = 293.4663),
    CLRK80IGN : "Clarke 1880 (IGN)." => (a = 6_378_249.2, invf = 293.4660212936269),
    CPM : "Comm. des Poids et Mesures 1799" => (a = 6_375_738.7, invf = 334.29),
    DELMBR : "Delambre 1810 (Belgium)" => (a = 6_376_428., invf = 311.5),
    ENGELIS : "Engelis 1985" => (a = 6_378_136.05, invf = 298.2566),
    EVRST30 : "Everest 1830" => (a = 6_377_276.345, invf = 300.8017),
    EVRST48 : "Everest 1948" => (a = 6_377_304.063, invf = 300.8017),
    EVRST56 : "Everest 1956" => (a = 6_377_301.243, invf = 300.8017),
    EVRST69 : "Everest 1969" => (a = 6_377_295.664, invf = 300.8017),
    EVRSTSS : "Everest (Sabah & Sarawak)" => (a = 6_377_298.556, invf = 300.8017),
    FSCHR60 : "Fischer (Mercury Datum) 1960" => (a = 6_378_166., invf = 298.3),
    FSCHR68 : "Fischer 1968" => (a = 6_378_150., invf = 298.3),
    GSK2011 : "GSK-2011" => (a = 6_378_136.5, invf = 298.2564151),
    GRS67 : "GRS 1967(IUGG, 1967)" => (a = 6_378_160., invf = 298.2471674270),
    GRS80 : "GRS 1980(IUGG, 1980)" => (a = 6_378_137., invf = 298.257222101),
    IAU76 : "IAU 1976" => (a = 6_378_140., invf = 298.257),
    HELMERT : "Helmert 1906" => (a = 6_378_200., invf = 298.3),
    HOUGH : "Hough" => (a = 6_378_270., invf = 297.),
    INTL : "International 1924 (Hayford 1909, 1910)" => (a = 6_378_388., invf = 297.),
    KAULA : "Kaula 1961" => (a = 6_378_163., invf = 298.24),
    KRASS: "Krassovsky, 1942" => (a = 6_378_245., invf = 298.3),
    LERCH : "Lerch 1979" => (a = 6_378_139., invf = 298.257),
    MPRTS : "Maupertuis 1738" => (a = 6_397_300., invf = 191.),
    MERIT : "MERIT 1983" => (a = 6_378_137., invf = 298.257),
    MOD_AIRY : "Modified Airy" => (a = 6_377_340.189, b = 6_356_034.446),
    FSCHR60M : "Modified Fischer 1960" => (a = 6_378_155., invf = 298.3),
    NWL9D : "Naval Weapons Lab., 1965" => (a = 6_378_145., invf = 298.25),
    NEW_INTL : "New International 1967" => (a = 6_378_157.5, b = 6_356_772.2),
    PLESSIS : "Plessis 1817 (France)" => (a = 6_376_523., b = 6_355_863.),
    PZ90 : "PZ-90" => (a = 6_378_136., invf = 298.25784),
    SEASIA : "Southeast Asia" => (a = 6_378_155., b = 6_356_773.320_5),
    SGS85 : "Soviet Geodetic System 85" => (a = 6_378_136., invf = 298.257),
    WALBECK : "Walbeck" => (a = 6_376_896., b = 6_355_834.846_7),
    WGS60 : "WGS 60" => (a = 6_378_165., invf = 298.3),
    WGS66 : "WGS 66" => (a = 6_378_145., invf = 298.25),
    WGS72 : "WGS 72" => (a = 6_378_135., invf = 298.26),
    WGS84 : "WGS 84" => (a = 6_378_137., invf = 298.257_223_563),
    SPHERE : "Normal Sphere (r=6370997)" => (a = 6_370_997., invf = 6_370_997.)
}

macro_rules! deg {
    ($d:expr) => {
        $d * std::f64::consts::PI / 180.0
    };
}

macro_rules! prime_meridians {
    ( $( $name:ident : $id:literal => $gw_lon_rad:expr ),+ ) => {

        pub mod prime_meridians {
            $(pub const $name: &str = $id;)+
        }

        static PRIME_MERIDIANS: phf::OrderedMap<&str, PrimeMeridian> = phf_ordered_map! {
            $( $id => PrimeMeridian::new($gw_lon_rad),)+
        };
    };
}

prime_meridians! {
    GREENWICH : "Greenwich"  => 0.,
    LISBON : "Lisbon"     => deg!(-9.131906111111),
    PARIS : "Paris"      => deg!(2.337229166667),
    BOGOTA : "Bogota"     => deg!(-74.080916666667),
    MADRID : "Madrid"     => deg!(-3.687938888889),
    ROME : "Rome"       => deg!(12.452333333333),
    BERN : "Bern"       => deg!(7.439583333333),
    JAKARTA : "Jakarta"    => deg!(106.807719444444),
    FERRO : "Ferro"      => deg!(-17.666666666667),
    BRUSSELS : "Brussels"   => deg!(4.367975),
    STOCKHOLM : "Stockholm"  => deg!(18.058277777778),
    ATHENS : "Athens"     => deg!(23.7163375),
    OSLO : "Oslo"       => deg!(10.722916666667),
    COPENHAGEN : "Copenhagen" => deg!(12.57788)
}

static DATUMS: phf::OrderedMap<&str, fn() -> GeodeticDatum> = phf_ordered_map! {
    "WGS 84" => GeodeticDatum::default,
};

/// The default [`GeodesyProvider`] provides commonly used geodesy elements.
pub struct DefaultGeodesyProvider;

impl GeodesyProvider for DefaultGeodesyProvider {
    fn ellipsoid_ids(&self) -> impl Iterator<Item = &str> {
        ELLIPSOIDS.keys().cloned()
    }

    fn ellipsoid(&self, id: &str) -> Option<Ellipsoid> {
        ELLIPSOIDS.get(id).cloned()
    }

    fn prime_meridian_ids(&self) -> impl Iterator<Item = &str> {
        PRIME_MERIDIANS.keys().cloned()
    }

    fn prime_meridian(&self, id: &str) -> Option<PrimeMeridian> {
        PRIME_MERIDIANS.get(id).cloned()
    }

    fn datum_ids(&self) -> impl Iterator<Item = &str> {
        // No pre-defined datum definition yet.
        empty()
    }

    fn datum(&self, id: &str) -> Option<GeodeticDatum> {
        DATUMS.get(id).map(|f| (*f)())
    }
}

#[cfg(test)]
mod tests {
    use crate::providers::default::ellipsoids::WGS84;
    use crate::providers::default::prime_meridians::GREENWICH;

    use super::DefaultGeodesyProvider;
    use super::GeodesyProvider;

    #[test]
    fn ellipsoid() {
        let provider = DefaultGeodesyProvider;
        assert_eq!(provider.ellipsoid_ids().collect::<Vec<_>>().len(), 46);
        let wgs84 = provider.ellipsoid(WGS84);
        assert!(wgs84.is_some());
        assert!((wgs84.unwrap().a() - 6378137.0).abs() < 1e-5);
        assert!((wgs84.unwrap().invf() - 298.257223563).abs() < 1e-5);
    }

    #[test]
    fn prime_meridian() {
        let provider = DefaultGeodesyProvider;
        assert_eq!(provider.prime_meridian_ids().collect::<Vec<_>>().len(), 14);
        let gw = provider.prime_meridian(GREENWICH);
        assert!(gw.is_some());
        assert_eq!(gw.unwrap().lon(), 0.0);
    }
}
