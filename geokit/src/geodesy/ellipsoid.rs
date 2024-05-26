use std::fmt::{Debug, Display};

use smol_str::SmolStr;
use crate::cs::geodetic::Lon;
use crate::math::complex::Complex;
use crate::math::polynomial::Polynomial;
use crate::units::angle::{Angle, Radians};

use crate::units::length::Length;

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
    /// Create a new ellipsoid using semi-major and semi-minor axes.
    ///
    /// # Panics
    ///
    /// If `semi_minor_axis` is negative or zero or if `semi_major_axis` is less than `semi_minor_axis`.
    pub fn from_ab<U: Length + Display>(
        name: &str,
        semi_major_axis: U,
        semi_minor_axis: U,
    ) -> Self {
        let a = semi_major_axis.m();
        let b = semi_minor_axis.m();
        assert!(b > 0., "Expected semi_minor_axis ({} m) > 0 m", b);
        assert!(
            a >= b,
            "Expected semi_major_axis ({} m) >= semi_minor_axis ({} m).",
            a,
            b
        );
        Self {
            name: SmolStr::new(name),
            a,
            b,
            invf: if a > b { a / (a - b) } else { f64::INFINITY },
        }
    }

    /// Create a new ellipsoid using semi-major axis and inverse flattening.
    ///
    /// # Panics
    ///
    /// if `semi_major_axis` is negative or zero or if `invf` is not greater than 1.
    pub fn from_ainvf<U: Length + Display>(name: &str, semi_major_axis: U, invf: f64) -> Self {
        let a = semi_major_axis.m();
        assert!(a > 0., "Expected semi_major_axis ({} m) > 0.", a);
        assert!(invf > 1., "Expected invf ({}) > 1.", invf);
        Self {
            name: SmolStr::new(name),
            a,
            b: a * (1. - 1. / invf),
            invf,
        }
    }

    /// Creates a new [Ellipsoid] with the given precomputed elements:
    /// - `semi_major_axis`: semi-major axis length **in meters**,
    /// - `semi_minor_axis`: semi-minor axis length **in meters**. **Must be less than or equal to `semi_major_axis`**
    /// - `invf`: inverse-flattening. **Must be greater than 1, possibly infinite if `semi_major_axis == semi_minor_axis`**
    #[inline(always)]
    pub(crate) const fn new_static(
        name: &'static str,
        semi_major_axis: f64,
        semi_minor_axis: f64,
        invf: f64,
    ) -> Self {
        // TODO: const assert ?
        Self {
            name: SmolStr::new_static(name),
            a: semi_major_axis,
            b: semi_minor_axis,
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
    pub fn a(&self) -> f64 {
        self.a
    }

    /// Return the semi-major axis squared.
    #[inline]
    pub fn a_sq(&self) -> f64 {
        self.a * self.a
    }

    /// Return the semi-minor axis length **in meters**.
    #[inline]
    pub fn b(&self) -> f64 {
        self.b
    }

    /// Return the semi minor axis squared.
    #[inline]
    pub fn b_sq(&self) -> f64 {
        self.b * self.b
    }

    /// Return the flattening.
    #[inline]
    pub fn f(&self) -> f64 {
        1. / self.invf
    }
    /// Return the inverse flattening: `a / (a - b)`.
    #[inline]
    pub fn invf(&self) -> f64 {
        self.invf
    }

    /// Return the 3rd flattening: `n = (a - b) / (a + b)`
    #[inline]
    pub fn n(&self) -> f64 {
        (self.a - self.b) / (self.a + self.b)
    }

    /// Return the first eccentricity squared: `(a^2 -b^2) / a^2`
    #[inline]
    pub fn e_sq(&self) -> f64 {
        (self.a_sq() - self.b_sq()) / self.a_sq()
    }

    pub fn e(&self) -> f64 {
        self.e_sq().sqrt()
    }

    /// Return the second eccentricity squared.
    /// See [e_prime]
    pub fn e_prime_sq(&self) -> f64 {
        (self.a_sq() - self.b_sq()) / self.b_sq()
    }

    /// Return the second eccentricity: `(a^2 - b^2) / b^2`
    pub fn e_prime(&self) -> f64 {
        self.e_prime_sq().sqrt()
    }

    /// Return the reduce latitude, aka parametric latitude, `beta`: `tan(beta) = (1 - f)tan(lat)`
    /// where `lat` is the geodetic latitude.
    ///
    /// The reduced latitude of a point P on the ellipsoid is the angle at the center of a sphere
    /// tangent to the ellipsoid on the equator (radius = a), between the equatorial plane and a
    /// radius to a point Q on the sphere whose projection along the z-axis intersect the ellipsoid
    /// at P.
    pub fn reduced_latitude(&self, lat: f64) -> f64 {
        ((1. - self.f()) * lat.tan()).atan()
    }

    /// Return the geocentric latitude, `psi`: `tan(psi) = (1 - e^2)tan(lat)`
    /// where `lat` is the geodetic latitude.
    ///
    /// The geocentric latitude of a point P on the ellipsoid is the angle at the centre of the
    /// ellipsoid between the equatorial plane and a line to the point P.
    pub fn geocentric_latitude(&self, lat: f64) -> f64 {
        ((1. - self.e_sq()) * lat.tan()).atan()
    }

    /// Return the radius of curvature `N` **in meters** at the given geodetic latitude **in radians**
    /// of the prime vertical normal section, i.e. the normal section perpendicular to the meridional
    /// normal section.
    /// See [prime_meridional_radius]
    pub fn prime_vertical_radius(&self, lat: f64) -> f64 {
        self.a / (1.0 - self.e_sq() * lat.sin().powi(2)).sqrt()
    }

    /// Return the radius of curvature `M` **in meters** at the given geodetic latitude **in radians**
    /// of the meridional normal section, i.e. the normal section passing through the poles.
    pub fn prime_meridional_radius(&self, lat: f64) -> f64 {
        self.a * (1.0 - self.e_sq()) / (1.0 - self.e_sq() * lat.sin().powi(2)).powf(1.5)
    }

    pub fn conformal_sphere_radius(&self, lat: f64) -> f64 {
        (self.prime_meridional_radius(lat) * self.prime_vertical_radius(lat)).sqrt()
    }

    /// Compute the coordinates and forward azimuth of the point at `s12` meters away
    /// from `p1` following the geodesic in the azimuth `alpha1` at `p1`.
    ///
    /// # Parameters
    ///
    /// - ellipsoid: the ellipsoid on which the geodesic is solved
    /// - p1: the **normalized geodetic** coordinates of the starting point
    /// - alpha1: the azimuth **in radians** of the geodesic at `p1`,
    /// - s12: the distance **in meters** along the geodesic from `p1` of the returned point coordinates.
    pub fn solve_direct<A: Angle, S: Length>(&self, p1: &[f64], alpha1: A, s12: S) -> ([f64; 2], f64) {
        // Step 1: solve NEP_1 to give alpha0, sigma1 and omega1
        let lat1 = p1[1];
        let beta1 = self.reduced_latitude(lat1);
        let (sin_beta1, cos_beta1) = beta1.sin_cos();
        let (sin_alpha1, cos_alpha1) = alpha1.rad().sin_cos();
        // alpha0 is the azimuth at the equator of the spherical triangle NEP_1, N: north pole,
        // E: intersection of the geodesic and the equator
        // Karney - Algorithms for geodesics eqn 10
        let alpha0 = Complex::new(
            Complex::new(cos_alpha1, sin_alpha1 * sin_beta1).abs(),
            sin_alpha1 * cos_beta1).theta();
        let (sin_alpha0, cos_alpha0) = alpha0.sin_cos();
        // Karney - Algorithms for geodesics eqn 11
        let sigma1 = Complex::new(cos_alpha1 * cos_beta1, sin_beta1).theta();
        let (sin_sigma1, cos_sigma1) = sigma1.sin_cos();
        // Karney - Algorithms for geodesics eqn 12
        let omega1 = Complex::new(cos_sigma1, sin_alpha0 * sin_sigma1).theta();

        // Step 2: Determine sigma2
        // Karney - Algorithms for geodesics eqn 9
        let k = self.e_prime() * cos_alpha0;
        // Karney - Algorithms for geodesics eqn 16
        let t = (1. + k * k).sqrt();
        let epsilon = (t - 1.) / (t + 1.);

        let a1 = Self::a1(epsilon);
        // geodesic arc length from E to p1 in meters
        let s1 = self.b() * Self::i1(a1, epsilon, sigma1);

        // geodesic arc length from E to p2 in meters (p2 is what we are looking for)
        let s2 = s1 + s12.m();

        let tau2 = s2 / (self.b() * a1);
        // arc length on the auxiliary sphere from E to p2 in radians
        let sigma2 = Self::sigma(epsilon, tau2);
        let (sin_sigma2, cos_sigma2) = sigma2.sin_cos();

        // Step 3: Solve NEB to give alpha2, beta2 (hence p2.lat) and omega2
        // Karney - Algorithms for geodesics eqn 14
        let alpha2 = Complex::new(cos_alpha0 * cos_sigma2, sin_alpha0).theta();
        // Karney - Algorithms for geodesics eqn 13
        let beta2 = Complex::new(Complex::new(cos_alpha0 * cos_sigma2, sin_alpha0).abs(), cos_alpha0 * sin_sigma2).theta();
        // Karney - Algorithms for geodesics eqn 12
        let omega2 = Complex::new(cos_sigma2, sin_alpha0 * sin_sigma2).theta();

        // Step 4: finally determine lambda2 and lambda12
        let a3 = Self::a3(self.n(), epsilon);
        // Karney - Algorithms for geodesics eqn 8
        let lambda1 = omega1 - self.f() * sin_alpha0 * Self::i3(a3, self.n(), epsilon, sigma1);
        let lambda2 = omega2 - self.f() * sin_alpha0 * Self::i3(a3, self.n(), epsilon, sigma2);

        let lambda12 = lambda2 - lambda1;

        // from: tan(beta) = (1 - f)*tan(lat)
        let lat2 = (beta2.tan() / (1. - self.f())).atan();

        // wrap the longitude in [-pi, pi]
        ([Lon::new(Radians(p1[0] + lambda12)).rad(), lat2], alpha2)

    }

    /// From Karney - Algorithms for geodesics eqn 17:
    /// A1 = (1 + 1/4 eps^2 + 1/64 eps^4 + 1/256 eps^6 + ...) / (1 - eps)
    fn a1(epsilon: f64) -> f64 {
        Polynomial::new([1., 0., 1. / 4., 0., 1. / 64., 0., 1. / 256.]).eval_at(epsilon) / (1. - epsilon)
    }

    /// From Karney - Algorithms for geodesics eqn 15:
    /// I1(sigma) = A1 * (sigma + sum(1, inf, C1l * sin(2l * sigma)))
    fn i1(a1: f64, epsilon: f64, sigma: f64) -> f64 {
        let c1xs = [
            Polynomial::new([0., -0.5, 0., 3. / 16., 0., -1. / 32.0]).eval_at(epsilon),
            Polynomial::new([0., 0., -1. / 16., 0., 1. / 32., 0., -9. / 2048.]).eval_at(epsilon),
            Polynomial::new([0., 0., 0., -1. / 48., 0., 3. / 256.]).eval_at(epsilon),
            Polynomial::new([0., 0., 0., 0., -5. / 512., 0., 3. / 512.]).eval_at(epsilon),
            Polynomial::new([0., 0., 0., 0., 0., -7. / 1280.]).eval_at(epsilon),
            Polynomial::new([0., 0., 0., 0., 0., 0., -7. / 2048.]).eval_at(epsilon),
        ];

        a1 * (sigma + c1xs.into_iter().enumerate().map(|(ix, c1x)| c1x * (2. * sigma * (ix + 1) as f64).sin()).sum::<f64>())
    }

    /// From Karney - Algorithms for geodesics eqn 20:
    /// sigma = tau + sum(1, inf, Cp1l * sin(2l * tau)) where tau = s / (b * A1)
    fn sigma(epsilon: f64, tau: f64) -> f64 {
        let cp1xs = [
            Polynomial::new([0., 1. / 2., 0., -9. / 32., 0., 205. / 1536.]).eval_at(epsilon),
            Polynomial::new([0., 0., 5. / 16., 0., -37. / 96., 0., 1335. / 4096.]).eval_at(epsilon),
            Polynomial::new([0., 0., 0., 29. / 96., 0., -75. / 128.]).eval_at(epsilon),
            Polynomial::new([0., 0., 0., 0., 539. / 1536., -2391. / 2560.]).eval_at(epsilon),
            Polynomial::new([0., 0., 0., 0., 0., 3467. / 7680.]).eval_at(epsilon),
            Polynomial::new([0., 0., 0., 0., 0., 0., 38081. / 61440.]).eval_at(epsilon),
        ];

        tau + cp1xs.into_iter().enumerate().map(|(ix, cp1x)| cp1x * (2. * tau * (ix + 1) as f64).sin()).sum::<f64>()
    }

    /// From Karney - Algorithms for geodesics eqn 24:
    fn a3(n: f64, epsilon: f64) -> f64 {
        Polynomial::new([
            1.,                                                 // x^0
            -Polynomial::new([1. / 2., -1. / 2.]).eval_at(n),            // x^1
            -Polynomial::new([1. / 4., 1. / 8., -3. / 8.]).eval_at(n),    // x^2
            -Polynomial::new([1. / 16., 3. / 16., 1. / 16.]).eval_at(n), // x^3
            -Polynomial::new([3. / 64., 1. / 32.]).eval_at(n),            // x^4
            -3. / 128.,
        ]).eval_at(epsilon)
    }

    /// From Karney - Algorithms for geodesics eqn 23
    /// I3(sigma) = A3 * (sigma + sum(1, inf, C3l * sin(2l * sigma))
    fn i3(a3: f64, n: f64, epsilon: f64, sigma: f64) -> f64 {
        let c3xs = [
            Polynomial::new([
                0.,  // x^0
                Polynomial::new([1. / 4., -1. / 4.]).eval_at(n),             // x^1
                Polynomial::new([1. / 8., 0., -1. / 8.]).eval_at(n),             // x^2
                Polynomial::new([3. / 64., 3. / 64., -1. / 64.]).eval_at(n), // x^3
                Polynomial::new([5. / 128., 1. / 64.]).eval_at(n),           // x^4
                3. / 128.,                                         // x^5
            ]).eval_at(epsilon),
            Polynomial::new([
                0.,
                0.,
                Polynomial::new([1. / 16., -3. / 32., 1. / 32.]).eval_at(n),
                Polynomial::new([3. / 64., -1. / 32., -3. / 64.]).eval_at(n),
                Polynomial::new([3. / 128., 1. / 128.]).eval_at(n),
                5. / 256.,
            ]).eval_at(epsilon),
            Polynomial::new([
                0.,
                0.,
                0.,
                Polynomial::new([5. / 192., -3. / 64., 5. / 192.]).eval_at(n),
                Polynomial::new([3. / 128., -5. / 192.]).eval_at(n),
                7. / 512.,
            ]).eval_at(epsilon),
            Polynomial::new([
                0.,
                0.,
                0.,
                0.,
                Polynomial::new([7. / 512., -7. / 256.]).eval_at(n),
                7. / 512.,
            ]).eval_at(epsilon),
            Polynomial::new([0., 0., 0., 0., 0., 21. / 2560.]).eval_at(epsilon),
        ];

        a3 * (sigma + c3xs.into_iter().enumerate().map(|(ix, c3x)| c3x * (2. * sigma * (ix + 1) as f64).sin()).sum::<f64>())
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

    macro_rules! def {
        ( ($name:literal, a = $a:expr, b = $b:expr ) ) => {
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
        ( ($name:literal, a = $a:expr, invf = $invf:literal ) ) => {
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
    use crate::units::length::{M, Meters};
    use approx::assert_abs_diff_eq;
    use crate::units::angle::Degrees;

    use super::*;

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
            wgs84.prime_meridional_radius(0.0),
            wgs84.a() * (1. - wgs84.e_sq()),
            epsilon = 1e-10
        );
        assert_abs_diff_eq!(
            wgs84.prime_meridional_radius(std::f64::consts::FRAC_PI_2),
            wgs84.a_sq() / wgs84.b(),
            epsilon = 1e-8
        );
    }

    #[test]
    fn test_prime_vertical_radius() {
        let wgs84 = consts::WGS84;
        assert_abs_diff_eq!(wgs84.prime_vertical_radius(0.0), wgs84.a(), epsilon = 1e-10);
        assert_abs_diff_eq!(
            wgs84.prime_vertical_radius(std::f64::consts::FRAC_PI_2),
            wgs84.a_sq() / wgs84.b(),
            epsilon = 1e-8
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

    #[test]
    fn test_solve_direct() {
        let wgs84 = consts::WGS84;
        let ([lon, lat], alpha2) = wgs84.solve_direct(&[0.0, Degrees(40.).rad()], Degrees(30.), Meters(10_000_000.0));
        assert_abs_diff_eq!(Degrees::from_rad(lon), Degrees(137.84490004377), epsilon=1e-11);
        assert_abs_diff_eq!(Degrees::from_rad(lat), Degrees(41.79331020506), epsilon=1e-10);
        assert_abs_diff_eq!(Degrees::from_rad(alpha2), Degrees(149.09016931807), epsilon=1e-10);
    }
}
