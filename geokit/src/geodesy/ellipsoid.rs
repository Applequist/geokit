use std::fmt::Debug;

use smol_str::SmolStr;

/// An `Ellipsoid` is a mathematical surface defined by rotating an ellipse around
/// it semi-minor axis.
/// It is used in a [GeodeticDatum] as a model of the Earth surface.
#[derive(Debug, Clone)]
pub struct Ellipsoid {
    /// The name of this ellipsoid.
    name: SmolStr,
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
    pub fn from_ab(name: &str, a: f64, b: f64) -> Self {
        assert!(b > 0., "Expected b > 0. Got {}", b);
        assert!(a >= b, "Expected a >= b. Got a = {}, b = {}", a, b);
        Self {
            name: SmolStr::new(name),
            a,
            b,
            invf: if a > b { a / (a - b) } else { f64::INFINITY },
        }
    }

    /// Create a new ellipsoid using semi-major axis **in meters** and inverse flattening.
    ///
    /// # Panics
    ///
    /// if `a` is negative or zero or if `invf` is not greater than 1.
    pub fn from_ainvf(name: &str, a: f64, invf: f64) -> Self {
        assert!(a > 0., "Expected a > 0. Got {}", a);
        assert!(invf > 1., "Expected invf > 1. Got {}", invf);
        Self {
            name: SmolStr::new(name),
            a,
            b: a * (1. - 1. / invf),
            invf,
        }
    }

    /// Creates a new [Ellipsoid] with the given precomputed elements:
    /// - `a`: semi-major axis **in meters**,
    /// - `b`: semi-minor axis **in meters**. **Must be less then or equal to `a`**
    /// - `invf`: inverse-flattening. **Must be greater than 1, possibly infinite if `a == b`**
    #[inline(always)]
    pub(crate) const fn new_static(name: &'static str, a: f64, b: f64, invf: f64) -> Self {
        Self {
            name: SmolStr::new_static(name),
            a,
            b,
            invf,
        }
    }

    pub fn name(&self) -> &str {
        self.name.as_str()
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

impl PartialEq for Ellipsoid {
    fn eq(&self, other: &Self) -> bool {
        self.a == other.a && self.b == other.b && self.invf == other.invf
    }
}

pub mod consts {
    use super::Ellipsoid;

    macro_rules! def {
        ( ($name:literal, a = $a:literal, b = $b:literal ) ) => {
            Ellipsoid::new_static(
                $name,
                $a,
                $b,
                if $a != $b {
                    $a / ($a - $b)
                } else {
                    f64::INFINITY
                },
            )
        };
        ( ($name:literal, a = $a:literal, invf = $invf:literal ) ) => {
            Ellipsoid::new_static($name, $a, $a * (1. - 1. / $invf), $invf)
        };
    }

    macro_rules! ellipsoids {
        ( $( $name:ident => $def:tt),+ ) => {
                $(pub const $name: Ellipsoid = def!($def);)+
        }
    }

    ellipsoids! {
        AIRY => ("Airy 1830", a = 6_377_563.396, invf = 299.3249646),
        ANDRAE => ("Andrae 1876 (Den., Iclnd.)", a = 6_377_104.43, invf = 300.0),
        DANISH => ("Andrae 1876 (Denmark, Iceland)", a = 6_377_019.256_3, invf = 300.0),
        APL_4_9 => ("Appl. Physics. 1965", a = 6_378_137., invf = 298.25),
        AUST_SA => ("Australian Natl & S. Amer. 1969", a = 6_378_160., invf = 298.25),
        BESSEL => ("Bessel 1841", a = 6_377_397.155, invf = 299.1528128),
        BESS_NAM => ("Bessel 1841 (Namibia)", a = 6_377_483.865, invf = 299.1528128),
        CLRK66 => ("Clarke 1866", a = 6_378_206.4, b = 6_356_583.8),
        CLRK80 => ("Clarke 1880 mod.", a = 6_378_249.145, invf = 293.4663),
        CLRK80IGN => ("Clarke 1880 (IGN).", a = 6_378_249.2, invf = 293.4660212936269),
        CPM => ("Comm. des Poids et Mesures 1799", a = 6_375_738.7, invf = 334.29),
        DELMBR => ("Delambre 1810 (Belgium)", a = 6_376_428., invf = 311.5),
        ENGELIS => ("Engelis 1985", a = 6_378_136.05, invf = 298.2566),
        EVRST30 => ("Everest 1830", a = 6_377_276.345, invf = 300.8017),
        EVRST48 => ("Everest 1948", a = 6_377_304.063, invf = 300.8017),
        EVRST56 => ("Everest 1956", a = 6_377_301.243, invf = 300.8017),
        EVRST69 => ("Everest 1969", a = 6_377_295.664, invf = 300.8017),
        EVRSTSS => ("Everest (Sabah & Sarawak)", a = 6_377_298.556, invf = 300.8017),
        FSCHR60 => ("Fischer (Mercury Datum) 1960", a = 6_378_166., invf = 298.3),
        FSCHR68 => ("Fischer 1968", a = 6_378_150., invf = 298.3),
        GSK2011 => ("GSK-2011", a = 6_378_136.5, invf = 298.2564151),
        GRS67 => ("GRS 1967(IUGG, 1967)", a = 6_378_160., invf = 298.2471674270),
        GRS80 => ("GRS 1980(IUGG, 1980)", a = 6_378_137., invf = 298.257222101),
        IAU76 => ("IAU 1976", a = 6_378_140., invf = 298.257),
        HELMERT => ("Helmert 1906", a = 6_378_200., invf = 298.3),
        HOUGH => ("Hough", a = 6_378_270., invf = 297.),
        INTL => ("International 1924 (Hayford 1909, 1910)", a = 6_378_388., invf = 297.),
        KAULA => ("Kaula 1961", a = 6_378_163., invf = 298.24),
        KRASS => ("Krassovsky, 1942", a = 6_378_245., invf = 298.3),
        LERCH => ("Lerch 1979", a = 6_378_139., invf = 298.257),
        MPRTS => ("Maupertuis 1738", a = 6_397_300., invf = 191.),
        MERIT => ("MERIT 1983", a = 6_378_137., invf = 298.257),
        MOD_AIRY => ("Modified Airy", a = 6_377_340.189, b = 6_356_034.446),
        FSCHR60M => ("Modified Fischer 1960", a = 6_378_155., invf = 298.3),
        NWL9D => ("Naval Weapons Lab., 1965", a = 6_378_145., invf = 298.25),
        NEW_INTL => ("New International 1967", a = 6_378_157.5, b = 6_356_772.2),
        PLESSIS => ("Plessis 1817 (France)", a = 6_376_523., b = 6_355_863.),
        PZ90 => ("PZ-90", a = 6_378_136., invf = 298.25784),
        SEASIA => ("Southeast Asia", a = 6_378_155., b = 6_356_773.320_5),
        SGS85 => ("Soviet Geodetic System 85", a = 6_378_136., invf = 298.257),
        WALBECK => ("Walbeck", a = 6_376_896., b = 6_355_834.846_7),
        WGS60 => ("WGS 60", a = 6_378_165., invf = 298.3),
        WGS66 => ("WGS 66", a = 6_378_145., invf = 298.25),
        WGS72 => ("WGS 72", a = 6_378_135., invf = 298.26),
        WGS84 => ("WGS 84", a = 6_378_137., invf = 298.257_223_563),
        SPHERE => ("Normal Sphere (r=6370997)", a = 6_370_997., invf = 6_370_997.)
    }
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    #[should_panic = "Expected b > 0"]
    fn negative_b() {
        let _ellipsoid = Ellipsoid::from_ab("Negative b", 1.0, -0.5);
    }

    #[test]
    #[should_panic = "Expected b > 0"]
    fn zero_b() {
        let _e = Ellipsoid::from_ab("Zero b", 1.0, 0.0);
    }

    #[test]
    #[should_panic = "Expected a >= b"]
    fn a_less_then_b() {
        let _ellipsoid = Ellipsoid::from_ab("a < b", 1.0, 2.0);
    }

    #[test]
    #[should_panic = "Expected a > 0"]
    fn negative_a() {
        let _ellipsoid = Ellipsoid::from_ainvf("Negative a", -0.1, 297.0);
    }

    #[test]
    #[should_panic = "Expected invf > 1"]
    fn small_invf() {
        let _ellipsoid = Ellipsoid::from_ainvf("Small invf", 1., 0.5);
    }

    #[test]
    #[should_panic = "Expected invf > 1"]
    fn unit_invf() {
        let _ellipsoid = Ellipsoid::from_ainvf("invf = 1", 1., 1.);
    }

    #[test]
    fn clone() {
        let e = Ellipsoid::from_ab("Cloned", 1., 0.9);
        let cpy = e.clone();
        assert_eq!(e, cpy);
        let _a = e.a;
    }

    #[test]
    fn partial_eq() {
        let e = Ellipsoid::from_ab("E", 1., 0.9);
        let eq = Ellipsoid::from_ab("E'", 1., 0.9);
        assert!(e.eq(&eq));
        assert!(!e.ne(&eq));

        let e2 = Ellipsoid::from_ab("E2", 1.01, 0.9);
        assert!(!e.eq(&e2));
        assert!(e.ne(&e2));

        let e3 = Ellipsoid::from_ainvf("E3", 1., 10.0);
        assert!(!e.eq(&e3));
        assert!(e.ne(&e3));
    }
}
