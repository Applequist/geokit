#![allow(non_snake_case)]

use std::mem::swap;

use crate::cs::azimuth::Azimuth;
use crate::cs::geodetic::{Lat, Lon};
use crate::geodesy::geodesics::{Geodesic, GeodesicSolver};
use crate::geodesy::Ellipsoid;
use crate::math::complex::Complex;
use crate::math::polynomial::Polynomial;
use crate::math::utils::clenshaw_sin_sum;
use crate::math::Float;
use crate::quantities::angle::Angle;
use crate::quantities::length::{Arc, Length};
use crate::units::angle::RAD;

/// Solve direct and inverse geodesic problems on an ellipsoid using algorithms
/// for Karney - Algorithms for Geodesics.
#[derive(Debug, Clone)]
pub struct KarneyGeodesicSolver<'e> {
    ellipsoid: &'e Ellipsoid,
    A3: Polynomial<6>,
    C3L: [Polynomial<6>; 6],
}

impl<'e> KarneyGeodesicSolver<'e> {
    // Used to compute A1 Eq(17)
    const PA1: Polynomial<7> = Polynomial::new([1., 0., 1. / 4., 0., 1. / 64., 0., 1. / 256.]);
    // Used to comput I1 Eq(17)
    const C1L: [Polynomial<7>; 7] = [
        Polynomial::new([0.; 7]), // c10
        Polynomial::new([0., -0.5, 0., 3. / 16., 0., -1. / 32.0, 0.]),
        Polynomial::new([0., 0., -1. / 16., 0., 1. / 32., 0., -9. / 2048.]),
        Polynomial::new([0., 0., 0., -1. / 48., 0., 3. / 256., 0.]),
        Polynomial::new([0., 0., 0., 0., -5. / 512., 0., 3. / 512.]),
        Polynomial::new([0., 0., 0., 0., 0., -7. / 1280., 0.0]),
        Polynomial::new([0., 0., 0., 0., 0., 0., -7. / 2048.]),
    ];
    // Used to compute sigma Eq (21)
    const CP1L: [Polynomial<7>; 6] = [
        Polynomial::new([0., 1. / 2., 0., -9. / 32., 0., 205. / 1536., 0.]),
        Polynomial::new([0., 0., 5. / 16., 0., -37. / 96., 0., 1335. / 4096.]),
        Polynomial::new([0., 0., 0., 29. / 96., 0., -75. / 128., 0.]),
        Polynomial::new([0., 0., 0., 0., 539. / 1536., 0., -2391. / 2560.]),
        Polynomial::new([0., 0., 0., 0., 0., 3467. / 7680., 0.]),
        Polynomial::new([0., 0., 0., 0., 0., 0., 38081. / 61440.]),
    ];

    /// Create a new solver instance.
    pub fn new(ellipsoid: &'e Ellipsoid) -> Self {
        let n = ellipsoid.n();
        // From Karney - Algorithms for geodesics eqn 24:
        let A3 = Polynomial::new([
            1.,                                                               // x^0
            -Polynomial::new([1. / 2., -1. / 2.]).fast_eval_at(n),            // x^1
            -Polynomial::new([1. / 4., 1. / 8., -3. / 8.]).fast_eval_at(n),   // x^2
            -Polynomial::new([1. / 16., 3. / 16., 1. / 16.]).fast_eval_at(n), // x^3
            -Polynomial::new([3. / 64., 1. / 32.]).fast_eval_at(n),           // x^4
            -3. / 128.,                                                       // x^5
        ]);
        let C3L = [
            // C30
            Polynomial::new([0.; 6]),
            // C31
            Polynomial::new([
                0.,                                                               // eps^0
                Polynomial::new([1. / 4., -1. / 4.]).fast_eval_at(n),             // eps^1
                Polynomial::new([1. / 8., 0., -1. / 8.]).fast_eval_at(n),         // eps^2
                Polynomial::new([3. / 64., 3. / 64., -1. / 64.]).fast_eval_at(n), // eps^3
                Polynomial::new([5. / 128., 1. / 64.]).fast_eval_at(n),           // eps^4
                3. / 128.,                                                        // eps^5
            ]),
            // C32
            Polynomial::new([
                0.,                                                               // eps^0
                0.,                                                               // eps^1
                Polynomial::new([1. / 16., -3. / 32., 1. / 32.]).fast_eval_at(n), // eps^2
                Polynomial::new([3. / 64., -1. / 32., -3. / 64.]).fast_eval_at(n),
                Polynomial::new([3. / 128., 1. / 128.]).fast_eval_at(n),
                5. / 256.,
            ]),
            // C33
            Polynomial::new([
                0.,
                0.,
                0.,
                Polynomial::new([5. / 192., -3. / 64., 5. / 192.]).fast_eval_at(n),
                Polynomial::new([3. / 128., -5. / 192.]).fast_eval_at(n),
                7. / 512.,
            ]),
            // C34
            Polynomial::new([
                0.,
                0.,
                0.,
                0.,
                Polynomial::new([7. / 512., -7. / 256.]).fast_eval_at(n),
                7. / 512.,
            ]),
            // C35
            Polynomial::new([0., 0., 0., 0., 0., 21. / 2560.]),
        ];
        Self { ellipsoid, A3, C3L }
    }

    pub fn direct(&self, p1: (Lon, Lat), alpha1: Azimuth, s12: Length) -> Geodesic {
        let (lon1, lat1) = p1;
        let beta1 = self.ellipsoid.reduced_latitude(lat1);

        // Solve triangle NEP1 given P1's reduced/parametric latitude and forward azimuth
        let nep1 = solve_triangle_from_p(beta1, alpha1);

        // Determine sigma2
        let (sin_alpha0, cos_alpha0) = nep1.alpha0.sin_cos();
        let k2 = self.ellipsoid.e_prime_sq() * cos_alpha0 * cos_alpha0;
        let t = (1. + k2).sqrt();
        let epsilon = (t - 1.) / (t + 1.);

        let A1: Float = Self::A1(epsilon);
        let I1_sigma1 = Self::I1(epsilon, A1, nep1.sigma);

        let s1: Length = self.ellipsoid.b() * I1_sigma1;
        let s2 = s1 + s12;

        let tau2 = Arc(s2) / (self.ellipsoid.b() * A1);
        let sigma2 = Self::Sigma(epsilon, tau2);

        // Solve triange NEP2
        let nep2 = solve_triangle_from_e(nep1.alpha0, sigma2);

        let A3 = self.A3(epsilon);
        let I3_sigma1 = self.I3(A3, epsilon, nep1.sigma);
        let I3_sigma2 = self.I3(A3, epsilon, nep2.sigma);

        let lambda1 = nep1.omega - self.ellipsoid.f() * sin_alpha0 * I3_sigma1 * RAD;
        let lambda2 = nep2.omega - self.ellipsoid.f() * sin_alpha0 * I3_sigma2 * RAD;
        let lambda12 = lambda2 - lambda1;

        // from: tan(beta) = (1 - f)*tan(lat)
        let lat2 = (nep2.beta.tan() / (1. - self.ellipsoid.f())).atan() * RAD;

        Geodesic {
            p1,
            alpha1,
            p2: (lon1 + lambda12, Lat::new(lat2)),
            alpha2: nep2.alpha,
            s: s12,
        }
    }

    pub fn inverse(&self, p1: (Lon, Lat), p2: (Lon, Lat)) -> Geodesic {
        let (lon1, mut lat1) = p1;
        let (lon2, mut lat2) = p2;

        // 1. swap endpoints and coordinates signs so:
        // a. 0 <= lambda12 <= pi
        // b. lat1 <= 0 and lat1 <= lat2 <= -lat1
        let mut lambda12 = lon2 - lon1;
        let mut lon_sign = lambda12.signum();
        lambda12 *= lon_sign;

        let swap_lat = if lat1.abs() < lat2.abs() { -1. } else { 1. };
        if swap_lat < 0. {
            lon_sign *= -1.;
            swap(&mut lat1, &mut lat2);
        }

        let lat_sign = -lat1.signum();
        lat1 *= lat_sign;
        lat2 *= lat_sign;

        assert!(
            Lon::ZERO <= lambda12 && lambda12 <= Lon::MAX,
            "lambda12 = {}",
            lambda12
        );
        assert!(Lat::MIN <= lat1 && lat1 <= Lat::ZERO, "lat1 = {}", lat1);
        assert!(lat1 <= lat2 && lat2 <= -lat1, "lat2 = {}", lat2);

        println!("lat1 = {lat1}");
        println!("lat2 = {lat2}");
        println!("lambda12 = {lambda12}");

        // 2.
        let beta1 = self.ellipsoid.reduced_latitude(lat1);
        let (sin_beta1, cos_beta1) = beta1.sin_cos();
        let beta2 = self.ellipsoid.reduced_latitude(lat2);
        let (sin_beta2, cos_beta2) = beta2.sin_cos();

        println!("beta1 = {}", beta1.deg());
        println!("beta2 = {}", beta2.deg());

        let alpha1: Azimuth;
        let alpha2: Azimuth;

        // 2. Special cases
        // a. Meridional: lambda12 = 0 or pi -> alpha1 = lambda12
        let is_meridional_special_case = lambda12 == Lon::ZERO || lambda12 == Lon::MAX;
        // b. Equatorial: lat1 = lat2 = 0 with lambda12 <= (1 -f) * pi -> alpha1 = pi / 2
        let is_equatorial_special_case =
            p1.1 == p2.1 && lambda12.angle() <= (1. - self.ellipsoid.f()) * Angle::PI;

        let sigma12: Angle;
        let mut s12: Length = Length::ZERO;
        if is_meridional_special_case {
            alpha1 = Azimuth::new(lambda12.angle());
            alpha2 = Azimuth::new(lambda12.angle());
        } else if is_equatorial_special_case {
            alpha1 = Azimuth::EAST;
            alpha2 = Azimuth::EAST;
        } else {
            let omega_bar =
                (1. - self.ellipsoid.e_sq() * (0.5 * (cos_beta1 + cos_beta2)).powi(2)).sqrt();
            let omega12 = lambda12.angle() / omega_bar;
            let (sin_omega12, cos_omega12) = omega12.sin_cos();

            // Eq (49)
            let z1 = Complex::new(
                cos_beta1 * sin_beta2 - sin_beta1 * cos_beta2 * cos_omega12,
                cos_beta2 * sin_omega12,
            );
            // Eq (50)
            let z2 = Complex::new(
                -sin_beta1 * cos_beta2 + cos_beta1 * sin_beta2 * cos_omega12,
                cos_beta1 * sin_omega12,
            );
            alpha1 = Azimuth::new(z1.arg());
            alpha2 = Azimuth::new(z2.arg());

            sigma12 = Complex::new(
                sin_beta1 * sin_beta2 + cos_beta1 * cos_beta2 * cos_omega12,
                z1.abs(),
            )
            .arg();

            s12 = omega_bar * (self.ellipsoid.a() * sigma12).length();

            println!("omega_bar = {omega_bar}");
            println!("omega12 = {}", omega12.deg());
            println!("sigma12 = {}", sigma12.deg());
            println!("alpha1 = {}", alpha1.angle().deg());
            println!("alpha2 = {}", alpha2.angle().deg());
            println!("s12 = {s12}");
        }

        Geodesic {
            p1,
            alpha1,
            p2,
            alpha2,
            s: s12,
        }
    }

    /// From Karney - Algorithms for geodesics eqn 17:
    /// A1 = (1 + 1/4 eps^2 + 1/64 eps^4 + 1/256 eps^6 + ...) / (1 - eps)
    fn A1(epsilon: Float) -> Float {
        Self::PA1.fast_eval_at(epsilon) / (1. - epsilon)
    }

    /// From Karney - Algorithms for geodesics eqn 15:
    /// I1(sigma) = A1 * (sigma + sum(1, inf, C1l * sin(2l * sigma)))
    fn I1(epsilon: Float, A1: Float, sigma: Angle) -> Float {
        let C1L = Self::C1L.map(|p| p.fast_eval_at(epsilon));
        A1 * (sigma.rad() + clenshaw_sin_sum(&C1L, 2. * sigma.rad()))
    }

    /// From Karney - Algorithms for geodesics eqn 20:
    /// sigma = tau + sum(1, inf, Cp1l * sin(2l * tau)) where tau = s / (b * A1)
    fn Sigma(epsilon: Float, tau: Angle) -> Angle {
        let cp1xs = Self::CP1L.map(|p| p.fast_eval_at(epsilon));

        tau + cp1xs
            .into_iter()
            .enumerate()
            .map(|(ix, cp1x)| cp1x * (2. * tau * (ix + 1) as Float).sin())
            .sum::<Float>()
            * RAD
    }

    /// From Karney - Algorithms for geodesics eqn 24:
    fn A3(&self, epsilon: Float) -> Float {
        self.A3.fast_eval_at(epsilon)
    }

    /// From Karney - Algorithms for geodesics eqn 23
    /// I3(sigma) = A3 * (sigma + sum(1, inf, C3l * sin(2l * sigma))
    fn I3(&self, A3: Float, epsilon: Float, sigma: Angle) -> Float {
        let C3L = self.C3L.map(|p| p.fast_eval_at(epsilon));
        A3 * (sigma.rad() + clenshaw_sin_sum(&C3L, 2. * sigma.rad()))
    }
}

/// A triangle NEP on the auxiliary sphere corresponding to the goedesic going
/// through a point P:
/// NE: meridian arc from the North pole (N) to intersection of the geodesic
/// going through P with the equator (E),
/// EP: great circle from E to P,
/// PN: meridian arc from P to N.
#[derive(Debug)]
struct Triangle {
    /// Azimuth of the geodesic arc EP at E
    alpha0: Azimuth,
    /// Azimuth of the geodesic arc EP at P
    alpha: Azimuth,
    /// reduced/parametric latitude of P
    beta: Angle,
    /// angular arc length of EP
    /// related to the length of EP on the ellipsoid `s` by: `s = b * I1(sigma)`
    sigma: Angle,
    /// angle from meridian plane containing NE and meridian plane containing PN
    /// related to the longitude difference on the ellipsoid by:
    /// `lambda = omega - f * sin(alpha0) * I3(sigma)`
    omega: Angle,
}

/// Given the reduced/parametric latitude of a Point P on the ellipsoid,
/// and the forward azimuth of a geodesic at P, solve the corresponding [Triangle]
/// on the auxiliary sphere.
fn solve_triangle_from_p(beta: Angle, alpha: Azimuth) -> Triangle {
    let (sin_alpha, cos_alpha) = alpha.sin_cos();
    let (sin_beta, cos_beta) = beta.sin_cos();

    // Eq (10)
    let alpha0 = Complex::new(
        Complex::new(cos_alpha, sin_alpha * sin_beta).abs(),
        sin_alpha * cos_beta,
    )
    .arg();

    // Eq (11)
    let sigma = Complex::new(cos_alpha * cos_beta, sin_beta).arg();
    let (sin_sigma, cos_sigma) = sigma.sin_cos();

    // Eq (12): omega1 is the ~longitude~ of p1 from the point of the geodesic on the equator
    let sin_alpha0 = alpha0.sin();
    let omega = Complex::new(cos_sigma, sin_alpha0 * sin_sigma).arg();

    Triangle {
        alpha0: Azimuth::new(alpha0),
        alpha,
        beta,
        sigma,
        omega,
    }
}

/// Given the forward azimuth of a geodesic crossing the equator with azimuth
/// `alpha0` and the angular arc length `sigma` of a great circle arc, solve the
/// cooresponding [Triangle] on the auxiliary sphere.
fn solve_triangle_from_e(alpha0: Azimuth, sigma: Angle) -> Triangle {
    let (sin_alpha0, cos_alpha0) = alpha0.sin_cos();
    let (sin_sigma, cos_sigma) = sigma.sin_cos();

    // Eq (14)
    let alpha = Complex::new(cos_alpha0 * cos_sigma, sin_alpha0).arg();

    // Eq (13)
    let beta = Complex::new(
        Complex::new(cos_alpha0 * cos_sigma, sin_alpha0).abs(),
        cos_alpha0 * sin_sigma,
    )
    .arg();

    // Eq (12)
    let omega = Complex::new(cos_sigma, sin_alpha0 * sin_sigma).arg();

    Triangle {
        alpha0,
        alpha: Azimuth::new(alpha),
        beta,
        sigma,
        omega,
    }
}

impl<'e> GeodesicSolver for KarneyGeodesicSolver<'e> {
    fn solve_direct(
        &self,
        p1: (Lon, Lat),
        alpha1: Azimuth,
        s12: Length,
    ) -> Result<Geodesic, &'static str> {
        Ok(self.direct(p1, alpha1, s12))
    }

    fn solve_inverse(&self, p1: (Lon, Lat), p2: (Lon, Lat)) -> Result<Geodesic, &'static str> {
        Ok(self.inverse(p1, p2))
    }
}

#[cfg(test)]
mod tests {
    use crate::cs::azimuth::Azimuth;
    use crate::cs::geodetic::{Lat, Lon};
    use crate::geodesy::ellipsoid::consts::{self, WGS84};
    use crate::geodesy::geodesics::karney::KarneyGeodesicSolver;
    use crate::geodesy::geodesics::tests::{
        antipodal_lines, check_direct, check_inverse, equatorial_lines, geographiclib_lines,
        meridional_lines, standard_lines, DirectError, InverseError, LineData,
    };
    use crate::geodesy::geodesics::{Geodesic, GeodesicSolver};
    use crate::math::{PI, PI_2};
    use crate::quantities::angle::Angle;
    use crate::quantities::length::Length;
    use crate::units::angle::{DEG, RAD};
    use crate::units::length::M;
    use approx::assert_abs_diff_eq;

    fn test_on(tset: LineData, err_direct: &DirectError, _err_inverse: &InverseError) {
        let solver = KarneyGeodesicSolver::new(&tset.ellipsoid);
        for tcase in tset.testcases.into_iter() {
            let direct = solver
                .solve_direct(tcase.p1, tcase.alpha1, tcase.s)
                .unwrap();
            check_direct(&direct, &tcase, err_direct);

            // TODO: implement inverse method.
            //let inverse = solver.solve_inverse(tcase.p1, tcase.p2).unwrap();
            //check_inverse(&inverse, &tcase, err_inverse);
        }
    }

    #[test]
    fn on_geographiclib_lines() {
        let tset = geographiclib_lines();
        test_on(tset, &DirectError::default(), &InverseError::default());
    }

    #[test]
    fn on_standard_lines() {
        let tset = standard_lines();
        test_on(tset, &DirectError::default(), &InverseError::default());
    }

    #[test]
    fn on_equatorial_lines() {
        let tset = equatorial_lines();
        test_on(tset, &DirectError::default(), &InverseError::default());
    }

    #[test]
    fn on_meridional_lines() {
        let tset = meridional_lines();
        test_on(tset, &DirectError::default(), &InverseError::default());
    }

    #[test]
    fn on_antipodal_lines() {
        let tset = antipodal_lines();
        test_on(tset, &DirectError::default(), &InverseError::default());
    }

    #[test]
    fn direct_karney_sample() {
        let wgs84 = consts::WGS84;
        dbg!(wgs84.a());
        dbg!(wgs84.f());
        dbg!(wgs84.b());
        dbg!(wgs84.n());
        dbg!(wgs84.e_sq());
        dbg!(wgs84.e_prime_sq());

        let solver = KarneyGeodesicSolver::new(&wgs84);
        let computed = solver
            .solve_direct(
                (Lon::new(0. * DEG), Lat::new(40. * DEG)),
                Azimuth::new(30. * DEG),
                Length::new(10_000_000., M),
            )
            .unwrap();
        check_direct(
            &computed,
            &Geodesic {
                p2: (
                    Lon::new(137.844_900_043_77 * DEG),
                    Lat::new(41.793_310_205_06 * DEG),
                ),
                alpha2: Azimuth::new(149.090_169_318_07 * DEG),
                ..computed
            },
            &DirectError::default(),
        );
        assert_abs_diff_eq!(
            computed.p2.0,
            Lon::new(137.844_900_043_77 * DEG),
            epsilon = Angle::tiny()
        );
        assert_abs_diff_eq!(
            computed.p2.1,
            Lat::new(41.793_310_205_06 * DEG),
            epsilon = Angle::tiny()
        );
        assert_abs_diff_eq!(
            computed.alpha2,
            Azimuth::new(149.090_169_318_07 * DEG),
            epsilon = Angle::tiny()
        );
    }

    #[test]
    fn inverse_karney_sample_short() {
        let wgs84 = consts::WGS84;
        let solver = KarneyGeodesicSolver::new(&wgs84);
        let expected = Geodesic {
            p1: (Lon::new(0. * DEG), Lat::new(-30.12345 * DEG)),
            alpha1: Azimuth::new(77.043_533_542_37 * DEG),
            p2: (Lon::new(0.000_05 * DEG), Lat::new(-30.12344 * DEG)),
            alpha2: Azimuth::new(77.043_508_449_13 * DEG),
            s: 4.9444_208 * M,
        };
        let computed = solver.solve_inverse(expected.p1, expected.p2).unwrap();
        check_inverse(&computed, &expected, &InverseError::default());
    }

    #[test]
    fn inverse_karney_sample_antipodal() {
        let wgs84 = consts::WGS84;
        let solver = KarneyGeodesicSolver::new(&wgs84);
        let expected = Geodesic {
            p1: (Lon::ZERO, Lat::new(-30. * DEG)),
            alpha1: Azimuth::new(161.890_524_736_33 * DEG),
            p2: (Lon::new(179.8 * DEG), Lat::new(29.9 * DEG)),
            alpha2: Azimuth::new(18.090_737_245_74 * DEG),
            s: 19_989_832.82761 * M,
        };
        let computed = solver.solve_inverse(expected.p1, expected.p2).unwrap();
        check_inverse(&computed, &expected, &InverseError::default());
    }
}
