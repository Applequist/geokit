use crate::{
    cs::geodetic::Lat,
    math::Float,
    quantities::length::Length,
    units::{angle::RAD, length::M},
};
use derive_more::derive::Display;
use smol_str::SmolStr;
use std::fmt::Debug;

/// An `Ellipsoid` is a mathematical surface defined by rotating an ellipse around
/// it semi-minor axis.
/// It is used in a [GeodeticDatum] as a model of the Earth surface.
#[derive(Debug, Clone, Display)]
#[display("(name = {}, a = {}, b = {}, inv_f = {})", name, a, b, invf)]
pub struct Ellipsoid {
    /// The name of this ellipsoid.
    name: SmolStr,
    /// The semi major axis length,
    a: Length,
    /// The semi minor axis length.
    b: Length,
    /// The inverse flattening: `a / (a - b)``. INFINITY if a == b.
    invf: Float,
}

impl Ellipsoid {
    /// Create a new ellipsoid using semi-major and semi-minor axes **in meters**.
    ///
    /// # Panics
    ///
    /// If `a` is negative or zero or if `a` is less than `b`.
    pub fn from_ab(name: &str, a: Length, b: Length) -> Self {
        assert!(b.m() > 0., "Expected semi_minor_axis ({}) > 0 m", b);
        assert!(
            a >= b,
            "Expected semi_major_axis ({}) >= semi_minor_axis ({}).",
            a,
            b
        );
        Self {
            name: SmolStr::new(name),
            a,
            b,
            invf: if a > b { a / (a - b) } else { Float::INFINITY },
        }
    }

    /// Create a new ellipsoid using semi-major axis in meters and inverse flattening.
    ///
    /// # Panics
    ///
    /// if `a` is negative or zero or if `invf` is not greater than 1.
    pub fn from_ainvf(name: &str, a: Length, invf: Float) -> Self {
        assert!(a > Length::ZERO, "Expected semi_major_axis ({}) > 0.", a);
        assert!(invf > 1., "Expected invf ({}) > 1.", invf);
        Self {
            name: SmolStr::new(name),
            a,
            b: a * (1. - 1. / invf),
            invf,
        }
    }

    /// Creates a new [Ellipsoid] with the given precomputed elements:
    /// - `a` semi-major axis length,
    /// - `b`: semi-minor axis length. **Must be less than or equal to `a`**
    /// - `invf`: inverse-flattening. **Must be greater than 1, possibly infinite if `a == b`**
    #[inline(always)]
    pub(crate) const fn new_static(name: &'static str, a: Length, b: Length, invf: Float) -> Self {
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

    /// Return the semi-major axis length **in meters**.
    #[inline]
    pub fn a(&self) -> Length {
        self.a
    }

    /// Return the semi-major axis squared.
    #[inline]
    pub fn a_sq(&self) -> Float {
        self.a.m() * self.a.m()
    }

    /// Return the semi-minor axis length **in meters**.
    #[inline]
    pub fn b(&self) -> Length {
        self.b
    }

    /// Return the semi minor axis squared.
    #[inline]
    pub fn b_sq(&self) -> Float {
        self.b.m() * self.b.m()
    }

    /// Return the flattening.
    #[inline]
    pub fn f(&self) -> Float {
        1. / self.invf
    }
    /// Return the inverse flattening: `a / (a - b)`.
    #[inline]
    pub fn invf(&self) -> Float {
        self.invf
    }

    /// Return the 3rd flattening: `n = (a - b) / (a + b)`
    #[inline]
    pub fn n(&self) -> Float {
        (self.a - self.b) / (self.a + self.b)
    }

    /// Return the first eccentricity squared: `(a^2 -b^2) / a^2`.
    #[inline]
    pub fn e_sq(&self) -> Float {
        (self.a_sq() - self.b_sq()) / self.a_sq()
    }

    /// Return the first eccentricity.
    /// See [e_sq]
    pub fn e(&self) -> Float {
        self.e_sq().sqrt()
    }

    /// Return the second eccentricity squared: `(a^2 - b^2) / b^2`.
    pub fn e_prime_sq(&self) -> Float {
        (self.a_sq() - self.b_sq()) / self.b_sq()
    }

    /// Return the second eccentricity.
    /// See [e_prime_sq]
    pub fn e_prime(&self) -> Float {
        self.e_prime_sq().sqrt()
    }

    /// Return the reduced latitude, aka parametric latitude, `beta`: `tan(beta) = (1 - f)tan(lat)`
    /// where `lat` is the **geodetic latitude** in radians.
    ///
    /// The reduced latitude of a point P on the ellipsoid is the angle at the center of a sphere
    /// tangent to the ellipsoid on the equator (radius = a), between the equatorial plane and a
    /// radius to a point Q on the sphere whose projection along the z-axis intersect the ellipsoid
    /// at P.
    pub fn reduced_latitude(&self, lat: Lat) -> Lat {
        Lat::new(((1. - self.f()) * lat.tan()).atan() * RAD)
    }

    /// Return the geocentric latitude, `psi`: `tan(psi) = (1 - e^2)tan(lat)`
    /// where `lat` is the **geodetic latitude** in radians.
    ///
    /// The geocentric latitude of a point P on the ellipsoid is the angle at the centre of the
    /// ellipsoid between the equatorial plane and a line to the point P.
    pub fn geocentric_latitude(&self, lat: Lat) -> Lat {
        Lat::new(((1. - self.e_sq()) * lat.tan()).atan() * RAD)
    }

    /// Return the radius of curvature `N` **in meters** at the given **geodetic latitude in radians**
    /// of the prime vertical normal section, i.e. the normal section perpendicular to the meridional
    /// normal section.
    /// See [prime_meridional_radius]
    pub fn prime_vertical_radius(&self, lat: Lat) -> Length {
        self.a / (1.0 - self.e_sq() * lat.sin().powi(2)).sqrt()
    }

    /// Return the radius of curvature `M` **in meters** at the given **geodetic latitude in radians**
    /// of the meridional normal section, i.e. the normal section passing through the poles.
    pub fn prime_meridional_radius(&self, lat: Lat) -> Length {
        self.a * (1.0 - self.e_sq()) / (1.0 - self.e_sq() * lat.sin().powi(2)).powf(1.5)
    }

    pub fn conformal_sphere_radius(&self, lat: Lat) -> Length {
        (self.prime_meridional_radius(lat).m() * self.prime_vertical_radius(lat).m()).sqrt() * M
    }
}

impl PartialEq for Ellipsoid {
    fn eq(&self, other: &Self) -> bool {
        self.a == other.a && self.b == other.b && self.invf == other.invf
    }
}

/// Well known ellipsoid definitions.
pub mod consts {
    use super::Ellipsoid;
    use crate::math::Float;
    use crate::quantities::length::Length;
    use crate::units::length::M;
    macro_rules! def {
        ( ($name:literal, a = $a:expr, b = $b:expr ) ) => {
            Ellipsoid::new_static(
                $name,
                Length::new($a, M),
                Length::new($b, M),
                if $a != $b {
                    $a / ($a - $b)
                } else {
                    Float::INFINITY
                },
            )
        };
        ( ($name:literal, a = $a:expr, invf = $invf:literal ) ) => {
            Ellipsoid::new_static(
                $name,
                Length::new($a, M),
                Length::new($a * (1. - 1. / $invf), M),
                $invf,
            )
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
    use crate::cs::geodetic::Lat;
    use crate::geodesy::ellipsoid::consts;
    use crate::geodesy::Ellipsoid;
    use crate::units::length::M;
    use approx::assert_abs_diff_eq;

    #[test]
    #[should_panic(expected = "Expected semi_minor_axis (-0.5 m) > 0")]
    fn negative_b() {
        let _ellipsoid = Ellipsoid::from_ab("Negative b", 1.0 * M, -0.5 * M);
    }

    #[test]
    #[should_panic(expected = "Expected semi_minor_axis (0 m) > 0")]
    fn zero_b() {
        let _e = Ellipsoid::from_ab("Zero b", 1.0 * M, 0.0 * M);
    }

    #[test]
    #[should_panic(expected = "Expected semi_major_axis (1 m) >= semi_minor_axis (2 m)")]
    fn a_less_then_b() {
        let _ellipsoid = Ellipsoid::from_ab("a < b", 1.0 * M, 2.0 * M);
    }

    #[test]
    #[should_panic(expected = "Expected semi_major_axis (-0.1 m) > 0")]
    fn negative_a() {
        let _ellipsoid = Ellipsoid::from_ainvf("Negative a", -0.1 * M, 297.0);
    }

    #[test]
    #[should_panic(expected = "Expected invf (0.5) > 1")]
    fn small_invf() {
        let _ellipsoid = Ellipsoid::from_ainvf("Small invf", 1. * M, 0.5);
    }

    #[test]
    #[should_panic(expected = "Expected invf (1) > 1")]
    fn unit_invf() {
        let _ellipsoid = Ellipsoid::from_ainvf("invf = 1", 1. * M, 1.);
    }

    #[test]
    fn test_prime_meridional_radius() {
        let wgs84 = consts::WGS84;
        assert_abs_diff_eq!(
            wgs84.prime_meridional_radius(Lat::ZERO),
            wgs84.a() * (1. - wgs84.e_sq())
        );
        assert_abs_diff_eq!(
            wgs84.prime_meridional_radius(Lat::MAX),
            wgs84.a_sq() / wgs84.b().m() * M,
        );
    }

    #[test]
    fn test_prime_vertical_radius() {
        let wgs84 = consts::WGS84;
        assert_abs_diff_eq!(wgs84.prime_vertical_radius(Lat::ZERO), wgs84.a(),);
        assert_abs_diff_eq!(
            wgs84.prime_vertical_radius(Lat::MAX),
            wgs84.a_sq() / wgs84.b().m() * M,
        );
    }

    #[test]
    fn clone() {
        let e = Ellipsoid::from_ab("Cloned", 1. * M, 0.9 * M);
        let cpy = e.clone();
        assert_eq!(e, cpy);
        let _a = e.a;
    }

    #[test]
    fn partial_eq() {
        let e = Ellipsoid::from_ab("E", 1. * M, 0.9 * M);
        let eq = Ellipsoid::from_ab("E'", 1. * M, 0.9 * M);
        assert!(e.eq(&eq));
        assert!(!e.ne(&eq));

        let e2 = Ellipsoid::from_ab("E2", 1.01 * M, 0.9 * M);
        assert!(!e.eq(&e2));
        assert!(e.ne(&e2));

        let e3 = Ellipsoid::from_ainvf("E3", 1. * M, 10.0);
        assert!(!e.eq(&e3));
        assert!(e.ne(&e3));
    }
}
