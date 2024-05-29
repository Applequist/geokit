use crate::cs::geodetic::Lon;
use crate::geodesy::Ellipsoid;
use crate::math::complex::Complex;
use crate::math::polynomial::Polynomial;
use crate::units::angle::Radians;

/// Solve direct and inverse geodesic problems on an ellipsoid.
#[derive(Debug, Clone)]
pub struct GeodesicSolver<'e> {
    ellipsoid: &'e Ellipsoid,
    a3: Polynomial<6>,
    c3xs: [Polynomial<6>; 5],
}

impl<'e> GeodesicSolver<'e> {
    /// Create a new solver instance.
    pub fn new(ellipsoid: &'e Ellipsoid) -> Self {
        let n = ellipsoid.n();
        // From Karney - Algorithms for geodesics eqn 24:
        let a3 = Polynomial::new([
            1.,                                                               // x^0
            -Polynomial::new([1. / 2., -1. / 2.]).eval_at(n),            // x^1
            -Polynomial::new([1. / 4., 1. / 8., -3. / 8.]).eval_at(n),   // x^2
            -Polynomial::new([1. / 16., 3. / 16., 1. / 16.]).eval_at(n), // x^3
            -Polynomial::new([3. / 64., 1. / 32.]).eval_at(n),           // x^4
            -3. / 128.,                                                       // x^5
        ]);
        let c3xs = [
            // C31
            Polynomial::new([
                0.,                                                               // x^0
                Polynomial::new([1. / 4., -1. / 4.]).eval_at(n),             // x^1
                Polynomial::new([1. / 8., 0., -1. / 8.]).eval_at(n),         // x^2
                Polynomial::new([3. / 64., 3. / 64., -1. / 64.]).eval_at(n), // x^3
                Polynomial::new([5. / 128., 1. / 64.]).eval_at(n),           // x^4
                3. / 128.,                                                        // x^5
            ]),
            // C32
            Polynomial::new([
                0.,
                0.,
                Polynomial::new([1. / 16., -3. / 32., 1. / 32.]).eval_at(n),
                Polynomial::new([3. / 64., -1. / 32., -3. / 64.]).eval_at(n),
                Polynomial::new([3. / 128., 1. / 128.]).eval_at(n),
                5. / 256.,
            ]),
            // C33
            Polynomial::new([
                0.,
                0.,
                0.,
                Polynomial::new([5. / 192., -3. / 64., 5. / 192.]).eval_at(n),
                Polynomial::new([3. / 128., -5. / 192.]).eval_at(n),
                7. / 512.,
            ]),
            // C34
            Polynomial::new([
                0.,
                0.,
                0.,
                0.,
                Polynomial::new([7. / 512., -7. / 256.]).eval_at(n),
                7. / 512.,
            ]),
            // C35
            Polynomial::new([0., 0., 0., 0., 0., 21. / 2560.]),
        ];
        Self {
            ellipsoid,
            a3,
            c3xs,
        }
    }

    /// Compute the coordinates and forward azimuth of the point at `s12` meters away
    /// from `p1` following the geodesic in the azimuth `alpha1` at `p1`.
    ///
    /// The algorithm is taken from Karney - Algorithms for geodesics.
    ///
    /// # Parameters
    ///
    /// - ellipsoid: the ellipsoid on which the geodesic is solved
    /// - p1: the **normalized geodetic** coordinates of the starting point
    /// - alpha1: the azimuth **in radians** of the geodesic at `p1`,
    /// - s12: the distance **in meters** along the geodesic from `p1` of the returned point coordinates.
    pub fn solve_direct(&self, p1: &[f64], alpha1: f64, s12: f64) -> ([f64; 2], f64) {
        // Step 1: solve NEP_1 to give alpha0, sigma1 and omega1
        let lat1 = p1[1];
        let beta1 = self.ellipsoid.reduced_latitude(lat1);
        let (sin_beta1, cos_beta1) = beta1.sin_cos();
        let (sin_alpha1, cos_alpha1) = alpha1.sin_cos();
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
        let k = self.ellipsoid.e_prime() * cos_alpha0;
        // Karney - Algorithms for geodesics eqn 16
        let t = (1. + k * k).sqrt();
        let epsilon = (t - 1.) / (t + 1.);

        let a1 = Self::a1(epsilon);
        // geodesic arc length from E to p1 in meters
        let s1 = self.ellipsoid.b() * Self::i1(a1, epsilon, sigma1);

        // geodesic arc length from E to p2 in meters (p2 is what we are looking for)
        let s2 = s1 + s12;

        let tau2 = s2 / (self.ellipsoid.b() * a1);
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
        let a3 = self.a3(epsilon);
        // Karney - Algorithms for geodesics eqn 8
        let lambda1 = omega1 - self.ellipsoid.f() * sin_alpha0 * self.i3(a3, epsilon, sigma1);
        let lambda2 = omega2 - self.ellipsoid.f() * sin_alpha0 * self.i3(a3, epsilon, sigma2);

        let lambda12 = lambda2 - lambda1;

        // from: tan(beta) = (1 - f)*tan(lat)
        let lat2 = (beta2.tan() / (1. - self.ellipsoid.f())).atan();

        // wrap the longitude in [-pi, pi]
        ([Lon::new(Radians(p1[0] + lambda12)).rad(), lat2], Radians(alpha2).wrap(0.))

    }

    pub fn solve_inverse(&self, p1: &[f64], p2: &[f64]) -> (f64, f64, f64) {
        // TODO:
        (0., 0., 0.)
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

    fn a2(epsilon: f64) -> f64 {
        Polynomial::new([1., 0., 1. / 4., 0., 9. / 64., 0., 25. / 256.]).eval_at(epsilon)
    }

    fn i2(a2: f64, epsilon: f64, sigma: f64) -> f64 {
        let c2xs = [
            Polynomial::new([0., 1./ 2., 0., 1. / 16., 0., 1. / 32.]).eval_at(epsilon),
            Polynomial::new([0., 0., 3. / 16., 0., 1. / 32., 0., 25. / 2048.]).eval_at(epsilon),
            Polynomial::new([0., 0., 0., 5. / 48., 0., 5. / 256.]).eval_at(epsilon),
            Polynomial::new([0., 0., 0., 0., 35. / 512., 0., 7. / 512.]).eval_at(epsilon),
            Polynomial::new([0., 0., 0., 0., 0., 63. / 1280.]).eval_at(epsilon),
            Polynomial::new([0., 0., 0., 0., 0., 0., 77. / 2048.]).eval_at(epsilon),
        ];

        a2 * (sigma + c2xs.into_iter().enumerate().map(|(ix, c2x)| c2x * (2. * sigma * (ix + 1) as f64).sin()).sum::<f64>())
    }

    /// From Karney - Algorithms for geodesics eqn 24:
    fn a3(&self, epsilon: f64) -> f64 {
        self.a3.eval_at(epsilon)
    }

    /// From Karney - Algorithms for geodesics eqn 23
    /// I3(sigma) = A3 * (sigma + sum(1, inf, C3l * sin(2l * sigma))
    fn i3(&self, a3: f64, epsilon: f64, sigma: f64) -> f64 {
        let c3xs = self.c3xs.map(|p| p.eval_at(epsilon));
        a3 * (sigma + c3xs.into_iter().enumerate().map(|(ix, c3x)| c3x * (2. * sigma * (ix + 1) as f64).sin()).sum::<f64>())
    }
}

#[cfg(test)]
mod tests {
    use approx::assert_abs_diff_eq;
    use crate::geodesy::ellipsoid::consts;
    use crate::geodesy::geodesics::GeodesicSolver;
    use crate::units::angle::{Angle, Degrees};
    use crate::units::length::{Length, Meters};

    #[test]
    fn test_solve_direct() {
        // From Karney - Algorithms for geodesics
        let wgs84 = consts::WGS84;
        let solver = GeodesicSolver::new(&wgs84);
        let ([lon, lat], alpha2) = solver.solve_direct(&[0.0, Degrees(40.).rad()], Degrees(30.).rad(), 10_000_000.0);
        assert_abs_diff_eq!(Degrees::from_rad(lon), Degrees(137.84490004377), epsilon=1e-11);
        assert_abs_diff_eq!(Degrees::from_rad(lat), Degrees(41.79331020506), epsilon=1e-10);
        assert_abs_diff_eq!(Degrees::from_rad(alpha2), Degrees(149.09016931807), epsilon=1e-10);

        // From Rapp - Geometric Geodesy 1.71 Standard Test Lines
        // Geodesics are given as (ph1, ph2, L, s, alpha12, alpha2)
        let intl = consts::INTL;
        let solver = GeodesicSolver::new(&intl);
        struct Geodesic {
            lat1: Degrees,
            lat2: Degrees,
            delta_lon: Degrees,
            s: Meters,
            alpha12: Degrees,
            alpha2: Degrees,
        }

        let standard_lines = [
            // Line 1
            Geodesic {
                lat1: Degrees::dms(37., 19., 54.95367),
                lat2: Degrees::dms(26., 7., 42.83946),
                delta_lon: Degrees::dms(41., 28., 35.50729),
                s: Meters(4_085_966.7026),
                alpha12: Degrees::dms(95., 27., 59.630888),
                alpha2: Degrees::dms(118., 5., 58.961608),
            },
            // Line 2
            Geodesic {
                lat1: Degrees::dms(35., 16., 11.24862),
                lat2: Degrees::dms(67., 22., 14.77638),
                delta_lon: Degrees::dms(137., 47., 28.31435),
                s: Meters(8_084_823.8383),
                alpha12: Degrees::dms(15., 44., 23.748498),
                alpha2: Degrees::dms(144., 55., 39.921473),
            },
            // Line 3
            Geodesic {
                lat1: Degrees::dms(1., 0., 0.),
                lat2: Degrees::dms(-0., 59., 53.83076),
                delta_lon: Degrees::dms(179., 17., 48.02997),
                s: Meters(19_959_999.9998),
                alpha12: Degrees::dms(88., 59., 59.998970),
                alpha2: Degrees::dms(91., 0., 6.118357),
            },
            // Line 4
            Geodesic {
                lat1: Degrees::dms(1., 0., 0.),
                lat2: Degrees::dms(1., 1., 15.18952),
                delta_lon: Degrees::dms(179., 46., 17.84244),
                s: Meters(19_780_006.5588),
                alpha12: Degrees::dms(4., 59., 59.999953),
                alpha2: Degrees::dms(174., 59., 59.884804),
            },
            // Line 5
            Geodesic {
                lat1: Degrees::dms(41., 41., 45.88000),
                lat2: Degrees::dms(41., 41., 46.20000),
                delta_lon: Degrees::dms(0., 0., 0.56000),
                s: Meters(16.2839751),
                alpha12: Degrees::dms(52., 40., 39.390667),
                alpha2: Degrees::dms(52., 40., 39.763168),
            },
            // Line 6
            Geodesic {
                lat1: Degrees::dms(30., 0., 0.),
                lat2: Degrees::dms(37., 53., 32.46584),
                delta_lon: Degrees::dms(116., 19., 16.68843),
                s: Meters(10002499.9999),
                alpha12: Degrees::dms(45., 0., 0.000004),
                alpha2: Degrees::dms(129., 8., 12.326010),
            },
            // Line 7
            Geodesic {
                lat1: Degrees::dms(37., 0., 0.),
                lat2: Degrees::dms(28., 15., 36.69535),
                delta_lon: Degrees::dms(-2., 37., 39.52918),
                s: Meters(1_000_000.0),
                alpha12: Degrees::dms(195., 0., 0.),
                alpha2: Degrees::dms(193., 34., 43.74060),
            }
        ];
        for (ix, g) in standard_lines.iter().enumerate() {
            println!("Test Line {}", ix + 1);
            let ([lon, lat], alpha2) = solver.solve_direct(&[0.0, g.lat1.rad()], g.alpha12.rad(), g.s.m());
            assert_abs_diff_eq!(Degrees::from_rad(lon), g.delta_lon, epsilon = 2e-9);
            assert_abs_diff_eq!(Degrees::from_rad(lat), g.lat2, epsilon = 2e-9);
            assert_abs_diff_eq!(Degrees::from_rad(alpha2), g.alpha2, epsilon = 1e-9);
        }

        // From Rapp - Geometric Geodesy Table 1.3
        let antipodal_lines = [
            // Line A
            Geodesic {
                lat1: Degrees::dms(41., 41., 45.88),
                lat2: Degrees::dms(-41., 41., 46.20),
                delta_lon: Degrees::dms(179., 59., 59.99985),
                s: Meters(20004566.7228),
                alpha12: Degrees::dms(179., 58., 49.16255),
                alpha2: Degrees::dms(0., 1., 10.8376),
            },
            // Line B
            Geodesic {
                lat1: Degrees(0.),
                lat2: Degrees(0.),
                delta_lon: Degrees(180.),
                s: Meters(19996147.4168),
                alpha12: Degrees::dms(29., 59., 59.9999),
                alpha2: Degrees(150.),
            },
            // Line C
            Geodesic {
                lat1: Degrees(30.),
                lat2: Degrees(-30.),
                delta_lon: Degrees(180.),
                s: Meters(19994364.6069),
                alpha12: Degrees::dms(39., 24., 51.8058),
                alpha2: Degrees::dms(140., 35., 8.1942),
            },
            // Line D
            Geodesic {
                lat1: Degrees(60.),
                lat2: Degrees::dms(-59., 59., 0.),
                delta_lon: Degrees::dms(179., 58., 53.03674),
                s: Meters(20000433.9629),
                alpha12: Degrees::dms(29., 11., 51.0700),
                alpha2: Degrees::dms(150., 49., 6.8680),
            },
            // Line E
            Geodesic {
                lat1: Degrees(30.),
                lat2: Degrees::dms(-29., 50., 0.),
                delta_lon: Degrees::dms(179., 56., 41.64754),
                s: Meters(19983420.1536),
                alpha12: Degrees::dms(16., 2., 28.3389),
                alpha2: Degrees::dms(163., 59., 10.3369),
            },
            // Line F
            Geodesic {
                lat1: Degrees(30.),
                lat2: Degrees::dms(-29., 55., 0.),
                delta_lon: Degrees::dms(179., 58., 3.57082),
                s: Meters(19992241.7634),
                alpha12: Degrees::dms(18., 38., 12.5568),
                alpha2: Degrees::dms(161., 22., 45.4373),
            }
        ];

        for (ix, g) in antipodal_lines.iter().enumerate() {
            println!("Test Line {}", ix + 1);
            let ([lon, lat], alpha2) = solver.solve_direct(&[0.0, g.lat1.rad()], g.alpha12.rad(), g.s.m());
            assert_abs_diff_eq!(Degrees::from_rad(lon), g.delta_lon, epsilon = 4e-1); // from 2e-9 to 4e-1 because of Line B
            assert_abs_diff_eq!(Degrees::from_rad(lat), g.lat2, epsilon = 2e-9);
            assert_abs_diff_eq!(Degrees::from_rad(alpha2), g.alpha2, epsilon = 1e-7);
        }

    }
}