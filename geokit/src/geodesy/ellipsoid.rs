/// An `Ellipsoid` is a mathematical surface used as a model of the Earth surface.
#[derive(Debug, Copy, Clone, PartialEq)]
pub struct Ellipsoid {
    /// The semi major axis length **in meters**.
    a: f64,
    /// The semi minor axis length **in meters**.
    b: f64,
    /// The inverse flattening: `a / (a - b)``. INFINITY if a == b.
    invf: f64,
}

impl Ellipsoid {
    /// Create a new ellipsoid using semi-major and semi-minor axes **in meters**.
    ///
    /// # Panics
    ///
    /// If `b` is negative or zero or if `a` is less then `b`.
    pub fn from_ab(a: f64, b: f64) -> Self {
        assert!(b > 0., "Expected b > 0. Got {}", b);
        assert!(a >= b, "Expected a >= b. Got a = {}, b = {}", a, b);
        Self::new(a, b, if a > b { a / (a - b) } else { f64::INFINITY })
    }

    /// Create a new ellipsoid using semi-major axis **in meters** and inverse flattening.
    ///
    /// # Panics
    ///
    /// if `a` is negative or zero or if `invf` is not greater than 1.
    pub fn from_ainvf(a: f64, invf: f64) -> Self {
        assert!(a > 0., "Expected a > 0. Got {}", a);
        assert!(invf > 1., "Expected invf > 1. Got {}", invf);
        Self::new(a, a * (1. - 1. / invf), invf)
    }

    pub const fn new(a: f64, b: f64, invf: f64) -> Self {
        Self { a, b, invf }
    }

    /// Return whether the ellipsoid is actually a sphere or not.
    pub fn is_spherical(&self) -> bool {
        self.invf.is_infinite()
    }

    /// Return the semi-major axis length in meters.
    #[inline]
    pub fn a(&self) -> f64 {
        self.a
    }

    /// Return the semi-major axis squared.
    #[inline]
    pub fn a_sq(&self) -> f64 {
        self.a * self.a
    }

    /// Return the semi-minor axis length in meters.
    #[inline]
    pub fn b(&self) -> f64 {
        self.b
    }

    /// Return the semi minor axis squared.
    #[inline]
    pub fn b_sq(&self) -> f64 {
        self.b * self.b
    }

    /// Return the inverse flattening: `a / (a - b)`.
    #[inline]
    pub fn invf(&self) -> f64 {
        self.invf
    }

    /// Return the first eccentricity squared: `(a^2 -b^2) / a^2`
    #[inline]
    pub fn e_sq(&self) -> f64 {
        (self.a_sq() - self.b_sq()) / self.a_sq()
    }

    /// Return the radius of curvature **in meters** in the east-west direction
    /// at the given latitude **in radians**.
    pub fn prime_vertical_radius(&self, lat: f64) -> f64 {
        self.a / (1.0 - self.e_sq() * lat.sin().powi(2)).sqrt()
    }

    /// Return the radius of curvature **in meters** in the north-south direction
    /// at the given latitude **in radians**.
    pub fn prime_meridional_radius(&self, lat: f64) -> f64 {
        self.a * (1.0 - self.e_sq()) / (1.0 - self.e_sq() * lat.sin().powi(2)).powf(1.5)
    }

    /// Convert **normalized geodetic coordinates** (lon in rad, lat in rad, height in meters)
    /// into **normalized geocentric coordinates** (x, y, z) all in meters.
    pub fn llh_to_xyz(&self, llh: &[f64], xyz: &mut [f64]) {
        let lon = llh[0];
        let lat = llh[1];
        let h = llh[2];

        let v = self.prime_vertical_radius(lat);
        let (sin_lon, cos_lon) = lon.sin_cos();
        let (sin_lat, cos_lat) = lat.sin_cos();

        xyz[0] = (v + h) * cos_lat * cos_lon;
        xyz[1] = (v + h) * cos_lat * sin_lon;
        xyz[2] = (v * (1.0 - self.e_sq()) + h) * sin_lat;
    }

    /// Convert **normalized geocentric coordinates** (x, y, z) in meters
    /// into **normalized geodetic coordinates** (lon in rad, lat in rad, height in meters)
    /// using Heiskanen and Moritz iterative method.
    pub fn xyz_to_llh(&self, xyz: &[f64], llh: &mut [f64]) {
        let x = xyz[0];
        let y = xyz[1];
        let z = xyz[2];

        let a2 = self.a_sq();
        let b2 = self.b_sq();
        let e2 = self.e_sq();

        let lon = y.atan2(x);

        let p = x.hypot(y);
        let mut lat = z.atan2(p * (1.0 - e2));
        let (sin_lat, cos_lat) = lat.sin_cos();
        let n = a2 / (a2 * cos_lat * cos_lat + b2 * sin_lat * sin_lat).sqrt();
        let mut h = p / cos_lat - n;
        loop {
            let next_lat = z.atan2(p * (1.0 - e2 * n / (n + h)));
            let (sin_nlat, cos_nlat) = next_lat.sin_cos();
            let next_n = a2 / ((a2 * cos_nlat * cos_nlat) + b2 * sin_nlat * sin_nlat).sqrt();
            let next_h = p / cos_nlat - next_n;
            let delta_lat = (lat - next_lat).abs();
            let delta_h = (h - next_h).abs();
            lat = next_lat;
            h = next_h;
            if delta_lat < 0.5e-5 && delta_h < 0.5e-3 {
                break;
            }
        }

        llh[0] = lon;
        llh[1] = lat;
        llh[2] = h;
    }
}

impl Default for Ellipsoid {
    /// Return the WGS84 (EPSG:7030) ellipsoid as default.
    fn default() -> Self {
        Ellipsoid::from_ainvf(6_378_137.0, 298.257_223_563)
    }
}

#[rustfmt::skip]
pub mod consts {
    use super::Ellipsoid;

    macro_rules! ellipsoid {
        ($name:ident, a = $a:expr, b = $b:expr, $desc:expr) => {
            pub const $name: Ellipsoid = Ellipsoid::new( $a, $b, if $a > $b { $a / ($a - $b) } else { f64::INFINITY });
        };
        ($name:ident, a = $a:expr, invf = $invf:expr, $desc:expr) => {
            pub const $name: Ellipsoid = Ellipsoid::new($a, $a * (1. - 1. / $invf), $invf);
        };
    }

    ellipsoid!(MERIT,     a = 6_378_137.,      invf = 298.257,           "MERIT 1983");
    ellipsoid!(SGS85,     a = 6_378_136.,      invf = 298.257,           "Soviet Geodetic System 85");
    ellipsoid!(GRS80,     a = 6_378_137.,      invf = 298.257222101,     "GRS 1980(IUGG, 1980)");
    ellipsoid!(IAU76,     a = 6_378_140.,      invf = 298.257,           "IAU 1976");
    ellipsoid!(AIRY,      a = 6_377_563.396,   invf = 299.3249646,       "Airy 1830");
    ellipsoid!(APL4_9,    a = 6_378_137.,      invf = 298.25,            "Appl. Physics. 1965");
    ellipsoid!(NWL9D,     a = 6_378_145.,      invf = 298.25,            "Naval Weapons Lab., 1965");
    ellipsoid!(MOD_AIRY,  a = 6_377_340.189,   b = 6_356_034.446,        "Modified Airy");
    ellipsoid!(ANDRAE,    a = 6_377_104.43,    invf = 300.0,             "Andrae 1876 (Den., Iclnd.)");
    ellipsoid!(DANISH,    a = 6_377_019.256_3, invf = 300.0,             "Andrae 1876 (Denmark, Iceland)");
    ellipsoid!(AUST_SA,   a = 6_378_160.,      invf = 298.25,            "Australian Natl & S. Amer. 1969");
    ellipsoid!(GRS67,     a = 6_378_160.,      invf = 298.2471674270,    "GRS 67(IUGG 1967)");
    ellipsoid!(GSK2011,   a = 6_378_136.5,     invf = 298.2564151,       "GSK-2011");
    ellipsoid!(BESSEL,    a = 6_377_397.155,   invf = 299.1528128,       "Bessel 1841");
    ellipsoid!(BESS_NAM,  a = 6_377_483.865,   invf = 299.1528128,       "Bessel 1841 (Namibia)");
    ellipsoid!(CLRK66,    a = 6_378_206.4,     b = 6_356_583.8,          "Clarke 1866");
    ellipsoid!(CLRK80,    a = 6_378_249.145,   invf = 293.4663,          "Clarke 1880 mod.");
    ellipsoid!(CLRK80IGN, a = 6_378_249.2,     invf = 293.4660212936269, "Clarke 1880 (IGN).");
    ellipsoid!(CPM,       a = 6_375_738.7,     invf = 334.29,            "Comm. des Poids et Mesures 1799");
    ellipsoid!(DELMBR,    a = 6_376_428.,      invf = 311.5,             "Delambre 1810 (Belgium)");
    ellipsoid!(ENGELIS,   a = 6_378_136.05,    invf = 298.2566,          "Engelis 1985");
    ellipsoid!(EVRST30,   a = 6_377_276.345,   invf = 300.8017,          "Everest 1830");
    ellipsoid!(EVRST48,   a = 6_377_304.063,   invf = 300.8017,          "Everest 1948");
    ellipsoid!(EVRST56,   a = 6_377_301.243,   invf = 300.8017,          "Everest 1956");
    ellipsoid!(EVRST69,   a = 6_377_295.664,   invf = 300.8017,          "Everest 1969");
    ellipsoid!(EVRSTSS,   a = 6_377_298.556,   invf = 300.8017,          "Everest (Sabah & Sarawak)");
    ellipsoid!(FSCHR60,   a = 6_378_166.,      invf = 298.3,             "Fischer (Mercury Datum) 1960");
    ellipsoid!(FSCHR60M,  a = 6_378_155.,      invf = 298.3,             "Modified Fischer 1960");
    ellipsoid!(FSCHR68,   a = 6_378_150.,      invf = 298.3,             "Fischer 1968");
    ellipsoid!(HELMERT,   a = 6_378_200.,      invf = 298.3,             "Helmert 1906");
    ellipsoid!(HOUGH,     a = 6_378_270.,      invf = 297.,              "Hough");
    ellipsoid!(INTL,      a = 6_378_388.,      invf = 297.,              "International 1924 (Hayford 1909, 1910)");
    ellipsoid!(KRASS,     a = 6_378_245.,      invf = 298.3,             "Krassovsky, 1942");
    ellipsoid!(KAULA,     a = 6_378_163.,      invf = 298.24,            "Kaula 1961");
    ellipsoid!(LERCH,     a = 6_378_139.,      invf = 298.257,           "Lerch 1979");
    ellipsoid!(MPRTS,     a = 6_397_300.,      invf = 191.,              "Maupertius 1738");
    ellipsoid!(NEW_INTL,  a = 6_378_157.5,     b = 6_356_772.2,          "New International 1967");
    ellipsoid!(PLESSIS,   a = 6_376_523.,      b = 6_355_863.,           "Plessis 1817 (France)");
    ellipsoid!(PZ90,      a = 6_378_136.,      invf = 298.25784,         "PZ-90");
    ellipsoid!(SEASIA,    a = 6_378_155.,      b = 6_356_773.320_5,      "Southeast Asia");
    ellipsoid!(WALBECK,   a = 6_376_896.,      b = 6_355_834.846_7,      "Walbeck");
    ellipsoid!(WGS60,     a = 6_378_165.,      invf = 298.3,             "WGS 60");
    ellipsoid!(WGS66,     a = 6_378_145.,      invf = 298.25,            "WGS 66");
    ellipsoid!(WGS72,     a = 6_378_135.,      invf = 298.26,            "WGS 72");
    ellipsoid!(WGS84,     a = 6_378_137.,      invf = 298.257_223_563,   "WGS 84");
    ellipsoid!(SPHERE,    a = 6_370_997.,      invf = 6_370_997.,        "Normal Sphere (r=6370997)");

}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    #[should_panic = "Expected b > 0"]
    fn negative_b() {
        let _ellipsoid = Ellipsoid::from_ab(1.0, -0.5);
    }

    #[test]
    #[should_panic = "Expected b > 0"]
    fn zero_b() {
        let _e = Ellipsoid::from_ab(1.0, 0.0);
    }

    #[test]
    #[should_panic = "Expected a >= b"]
    fn a_less_then_b() {
        let _ellipsoid = Ellipsoid::from_ab(1.0, 2.0);
    }

    #[test]
    #[should_panic = "Expected a > 0"]
    fn negative_a() {
        let _ellipsoid = Ellipsoid::from_ainvf(-0.1, 297.0);
    }

    #[test]
    #[should_panic = "Expected invf > 1"]
    fn small_invf() {
        let _ellipsoid = Ellipsoid::from_ainvf(1., 0.5);
    }

    #[test]
    #[should_panic = "Expected invf > 1"]
    fn unit_invf() {
        let _ellipsoid = Ellipsoid::from_ainvf(1., 1.);
    }

    #[test]
    fn copy() {
        let e = Ellipsoid::from_ab(1., 0.9);
        let cpy = e;
        assert_eq!(e, cpy);
        let _a = e.a;
    }

    #[test]
    fn partial_eq() {
        let e = Ellipsoid::from_ab(1., 0.9);
        let eq = e;
        assert!(e.eq(&eq));
        assert!(!e.ne(&eq));

        let e2 = Ellipsoid::from_ab(1.01, 0.9);
        assert!(!e.eq(&e2));
        assert!(e.ne(&e2));

        let e3 = Ellipsoid::from_ainvf(1., 10.0);
        assert!(!e.eq(&e3));
        assert!(e.ne(&e3));
    }

    #[test]
    fn default() {
        let wgs84 = Ellipsoid::default();
        assert_eq!(wgs84.a, 6_378_137.0);
        assert_eq!(wgs84.invf, 298.257_223_563);
    }
}
