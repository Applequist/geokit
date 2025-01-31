#![allow(non_snake_case)]

use crate::cs::azimuth::Azimuth;
use crate::cs::geodetic::{Lat, Lon};
use crate::geodesy::geodesics::{Geodesic, GeodesicSolver};
use crate::geodesy::Ellipsoid;
use crate::math::complex::Complex;
use crate::math::polynomial::Polynomial;
use crate::math::Float;
use crate::quantities::angle::Angle;
use crate::quantities::length::{Arc, Length};
use crate::units::angle::{DEG, RAD};

/// Solve direct and inverse geodesic problems on an ellipsoid using algorithms
/// for Karney - Algorithms for Geodesics.
#[derive(Debug, Clone)]
pub struct KarneyGeodesicSolver<'e> {
    ellipsoid: &'e Ellipsoid,
    A3: Polynomial<6>,
    C3L: [Polynomial<6>; 5],
}

impl<'e> KarneyGeodesicSolver<'e> {
    // Used to compute A1 Eq(17)
    const PA1: Polynomial<7> = Polynomial::new([1., 0., 1. / 4., 0., 1. / 64., 0., 1. / 256.]);
    // Used to comput I1 Eq(17)
    const C1L: [Polynomial<7>; 6] = [
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

    /// From Karney - Algorithms for geodesics eqn 17:
    /// A1 = (1 + 1/4 eps^2 + 1/64 eps^4 + 1/256 eps^6 + ...) / (1 - eps)
    fn A1(epsilon: Float) -> Float {
        Self::PA1.fast_eval_at(epsilon) / (1. - epsilon)
    }

    /// From Karney - Algorithms for geodesics eqn 15:
    /// I1(sigma) = A1 * (sigma + sum(1, inf, C1l * sin(2l * sigma)))
    fn I1(epsilon: Float, A1: Float, sigma: Angle) -> Float {
        let C1Ls = Self::C1L.map(|p| p.fast_eval_at(epsilon));
        A1 * (sigma.rad()
            + C1Ls
                .into_iter()
                .enumerate()
                .map(|(ix, c1x)| c1x * (2. * sigma * (ix + 1) as Float).sin())
                .sum::<Float>())
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
        A3 * (sigma.rad()
            + C3L
                .into_iter()
                .enumerate()
                .map(|(ix, c3x)| c3x * (2. * sigma * (ix + 1) as Float).sin())
                .sum::<Float>())
    }
}

#[derive(Debug)]
struct Triangle {
    alpha0: Azimuth, // azimuth at equator
    alpha: Azimuth,
    beta: Lat,
    sigma: Angle,
    omega: Angle,
}

fn solve_triangle_from_p(beta: Lat, alpha: Azimuth) -> Triangle {
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
        beta: Lat::new(beta),
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
        unimplemented!()
    }
}

#[cfg(test)]
mod tests {
    use crate::cs::azimuth::Azimuth;
    use crate::cs::geodetic::{Lat, Lon};
    use crate::geodesy::ellipsoid::consts;
    use crate::geodesy::geodesics::karney::KarneyGeodesicSolver;
    use crate::geodesy::geodesics::tests::{
        antipodal_lines, standard_lines, vincenty_direct_deltas, vincenty_inverse_deltas,
        vincenty_lines,
    };
    use crate::geodesy::geodesics::GeodesicSolver;
    use crate::math::{PI, PI_2};
    use crate::quantities::length::Length;
    use crate::units::angle::{DEG, RAD};
    use crate::units::length::M;
    use approx::assert_abs_diff_eq;

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
        assert_abs_diff_eq!(
            computed.p2.0,
            Lon::new(137.844_900_043_77 * DEG),
            epsilon = 1e-12
        );
        assert_abs_diff_eq!(
            computed.p2.1,
            Lat::new(41.793_310_205_06 * DEG),
            epsilon = 1e-12
        );
        assert_abs_diff_eq!(
            computed.alpha2,
            Azimuth::new(149.090_169_318_07 * DEG),
            epsilon = 1e-12
        );
    }

    #[test]
    fn direct_vincenty_lines() {
        let chars = ['a', 'b', 'c', 'd', 'e', 'f'];
        for (ix, (input, delta)) in vincenty_lines()
            .into_iter()
            .zip(vincenty_direct_deltas())
            .enumerate()
        {
            let solver = KarneyGeodesicSolver::new(&input.ellipsoid);
            let computed = solver
                .solve_direct(input.geodesic.p1, input.geodesic.alpha1, input.geodesic.s)
                .unwrap();
            println!("-----------------------------------------------------------------------------------");
            println!("Input: Vincenty Line ({})", chars[ix]);
            println!("{}", input.geodesic);
            println!("Computed: ");
            println!("{}", computed);
            println!();

            let diff_lon_dms = (computed.p2.0.angle() - input.geodesic.p2.0.angle()).to_dms();
            println!("error on lon (computed - input) = {}", diff_lon_dms);

            let diff_lat_dms = (computed.p2.1.angle() - input.geodesic.p2.1.angle()).to_dms();
            println!("error on lat (computed - input) = {}", diff_lat_dms);

            let diff_az_dms = (computed.alpha2.angle() - input.geodesic.alpha2.angle()).to_dms();
            println!("error on az (computed - input) = {}", diff_az_dms);

            assert_eq!(diff_lon_dms.deg(), 0.0);
            assert_eq!(diff_lon_dms.min(), 0.0);
            assert!(diff_lon_dms.sec() < delta.delta_delta_lon.abs());

            assert_eq!(diff_lat_dms.deg(), 0.0);
            assert_eq!(diff_lat_dms.min(), 0.0);
            assert!(diff_lat_dms.sec() < delta.delta_lat2.abs());

            assert_eq!(diff_az_dms.deg(), 0.0);
            assert_eq!(diff_az_dms.min(), 0.0);
            assert!(diff_az_dms.sec() < delta.delta_alpha2.abs());
        }
    }

    #[test]
    fn direct_standard_lines() {
        // Direct problem
        for (ix, input) in standard_lines().iter().enumerate() {
            let solver = KarneyGeodesicSolver::new(&input.ellipsoid);
            let computed = solver
                .solve_direct(input.geodesic.p1, input.geodesic.alpha1, input.geodesic.s)
                .unwrap();
            println!("-----------------------------------------------------------------------------------");
            println!("Input: Standard Line ({})", ix);
            println!("{}", input.geodesic);
            println!("Computed: ");
            println!("{}", computed);
            println!();

            let diff_lon_dms = (computed.p2.0.angle() - input.geodesic.p2.0.angle()).to_dms();
            println!("error on lon (computed - input) = {}", diff_lon_dms);

            let diff_lat_dms = (computed.p2.1.angle() - input.geodesic.p2.1.angle()).to_dms();
            println!("error on lat (computed - input) = {}", diff_lat_dms);

            let diff_az_dms = (computed.alpha2.angle() - input.geodesic.alpha2.angle()).to_dms();
            println!("error on az (computed - input) = {}", diff_az_dms);

            assert_abs_diff_eq!(computed.p2.0, input.geodesic.p2.0, epsilon = 3e-11);
            assert_abs_diff_eq!(computed.p2.1, input.geodesic.p2.1, epsilon = 3e-11);
            assert_abs_diff_eq!(computed.alpha2, input.geodesic.alpha2, epsilon = 1e-10);
        }
    }

    #[test]
    fn direct_antipodal_lines() {
        for (ix, input) in antipodal_lines().iter().enumerate() {
            let solver = KarneyGeodesicSolver::new(&input.ellipsoid);
            let computed = solver
                .solve_direct(input.geodesic.p1, input.geodesic.alpha1, input.geodesic.s)
                .unwrap();
            println!("-----------------------------------------------------------------------------------");
            println!("Input: Antipodal Line ({})", ix);
            println!("{}", input.geodesic);
            println!("Computed: ");
            println!("{}", computed);
            println!();

            let diff_lon_dms = (computed.p2.0.angle() - input.geodesic.p2.0.angle()).to_dms();
            println!("error on lon (computed - input) = {}", diff_lon_dms);

            let diff_lat_dms = (computed.p2.1.angle() - input.geodesic.p2.1.angle()).to_dms();
            println!("error on lat (computed - input) = {}", diff_lat_dms);

            let diff_az_dms = (computed.alpha2.angle() - input.geodesic.alpha2.angle()).to_dms();
            println!("error on az (computed - input) = {}", diff_az_dms);

            assert_abs_diff_eq!(computed.p2.0, input.geodesic.p2.0, epsilon = 1e-11);
            assert_abs_diff_eq!(computed.p2.1, input.geodesic.p2.1, epsilon = 1e-11);
            assert_abs_diff_eq!(computed.alpha2, input.geodesic.alpha2, epsilon = 1e-9);
        }
    }

    #[test]
    fn direct_equator_lines() {
        let wgs84 = consts::WGS84;
        let solver = KarneyGeodesicSolver::new(&wgs84);
        let computed = solver
            .solve_direct(
                (Lon::new(0.0 * RAD), Lat::new(0.0 * RAD)),
                Azimuth::new(90.0 * DEG),
                Length::new(20_000.0, M),
            )
            .unwrap();
        assert_abs_diff_eq!(
            computed.p2.0.rad(),
            (0.17966306 * DEG).rad(),
            epsilon = 1e-8
        );
        assert_abs_diff_eq!(computed.p2.1.rad(), 0.0, epsilon = 1e-8);
        assert_abs_diff_eq!(computed.alpha2.rad(), PI_2, epsilon = 1e-8);

        let computed = solver
            .solve_direct(
                (Lon::new(170.0 * DEG), Lat::new(0.0 * DEG)),
                Azimuth::new(90.0 * DEG),
                Length::new(2_000_000.0, M),
            )
            .unwrap();
        assert_abs_diff_eq!(
            computed.p2.0,
            Lon::new(-172.03369432 * DEG),
            epsilon = 1e-10
        );
        assert_abs_diff_eq!(computed.p2.1, Lat::new(0.0 * RAD), epsilon = 1e-10);
        assert_abs_diff_eq!(computed.alpha2, Azimuth::new(PI_2 * RAD), epsilon = 1e-8);
    }

    #[test]
    fn direct_meridian_lines() {
        let wgs84 = consts::WGS84;
        let solver = KarneyGeodesicSolver::new(&wgs84);

        let computed = solver
            .solve_direct(
                (Lon::new(0. * DEG), Lat::new(-10. * DEG)),
                Azimuth::new(0. * DEG),
                Length::new(2_000_000.0, M),
            )
            .unwrap();
        assert_abs_diff_eq!(computed.p2.0, Lon::new(0.0 * RAD), epsilon = 1e-10);
        assert_abs_diff_eq!(computed.p2.1, Lat::new(8.08583903 * DEG), epsilon = 1e-10);
        assert_abs_diff_eq!(computed.alpha2, Azimuth::new(0.0 * RAD), epsilon = 1e-8);

        let computed = solver
            .solve_direct(
                (Lon::new(0. * DEG), Lat::new(80. * DEG)),
                Azimuth::new(0. * DEG),
                Length::new(2_000_000.0, M),
            )
            .unwrap();
        assert_abs_diff_eq!(computed.p2.0, Lon::new(PI * RAD), epsilon = 1e-10);
        assert_abs_diff_eq!(computed.p2.1, Lat::new(82.09240627 * DEG), epsilon = 1e-10);
        assert_abs_diff_eq!(computed.alpha2, Azimuth::new(PI * RAD), epsilon = 1e-8);
    }

    #[test]
    fn inverse_karney_sample_short() {
        let wgs84 = consts::WGS84;
        let solver = KarneyGeodesicSolver::new(&wgs84);
        let computed = solver
            .solve_inverse(
                (Lon::new(0. * DEG), Lat::new(-30.12345 * DEG)),
                (Lon::new(0.000_05 * DEG), Lat::new(-30.12344 * DEG)),
            )
            .unwrap();
        assert_abs_diff_eq!(
            computed.alpha1,
            Azimuth::new(77.043_533_542_37 * DEG),
            epsilon = 1e-12
        );
        assert_abs_diff_eq!(
            computed.alpha2,
            Azimuth::new(77.043_508_449_13 * DEG),
            epsilon = 1e-12
        );
        assert_abs_diff_eq!(computed.s.m(), 4.944_208, epsilon = 1e-6);
    }

    #[test]
    fn inverse_karney_sample_antipodal() {
        let wgs84 = consts::WGS84;
        let solver = KarneyGeodesicSolver::new(&wgs84);
        let computed = solver
            .solve_inverse(
                (Lon::new(0.0 * DEG), Lat::new(-30. * DEG)),
                (Lon::new(179.8 * DEG), Lat::new(29.9 * DEG)),
            )
            .unwrap();
        assert_abs_diff_eq!(
            computed.alpha1,
            Azimuth::new(161.890_524_736_33 * DEG),
            epsilon = 1e-12
        );
        assert_abs_diff_eq!(
            computed.alpha2,
            Azimuth::new(18.090_737_245_74 * DEG),
            epsilon = 1e-12
        );
        assert_abs_diff_eq!(computed.s.m(), 19_989_832.82761, epsilon = 1e-6);
    }

    #[test]
    fn inverse_vincenty_lines() {
        // Vincenty lines
        let chars = ['a', 'b', 'c', 'd', 'e', 'f'];
        for (ix, (input, delta)) in vincenty_lines()
            .into_iter()
            .zip(vincenty_inverse_deltas())
            .enumerate()
        {
            let solver = KarneyGeodesicSolver::new(&input.ellipsoid);
            let computed = solver
                .solve_inverse(input.geodesic.p1, input.geodesic.p2)
                .unwrap();
            println!("-----------------------------------------------------------------------------------");
            println!("Input: Vincenty Line ({})", chars[ix]);
            println!("{}", input.geodesic);
            println!("Computed: ");
            println!("{}", computed);
            println!();

            let diff_az1_dms = (computed.alpha1.angle() - input.geodesic.alpha1.angle()).to_dms();
            println!("error on az1 (computed - input) = {}", diff_az1_dms);

            let diff_az2_dms = (computed.alpha2.angle() - input.geodesic.alpha2.angle()).to_dms();
            println!("error on az2 (computed - input) = {}", diff_az2_dms);

            let diff_s_m = computed.s - input.geodesic.s;
            println!("error on s (computed - input) = {}", diff_s_m);

            assert_eq!(diff_az1_dms.deg(), 0.0);
            assert_eq!(diff_az1_dms.min(), 0.0);
            assert!(diff_az1_dms.sec() < delta.delta_alpha1.abs());

            assert_eq!(diff_az2_dms.deg(), 0.0);
            assert_eq!(diff_az2_dms.min(), 0.0);
            assert!(diff_az2_dms.sec() < delta.delta_alpha2.abs());

            assert!(diff_s_m.m() < delta.delta_s.abs());
        }
    }

    #[test]
    fn inverse_standard_lines() {
        // Standard lines
        for (ix, input) in standard_lines().iter().enumerate() {
            let solver = KarneyGeodesicSolver::new(&input.ellipsoid);
            let computed = solver
                .solve_inverse(input.geodesic.p1, input.geodesic.p2)
                .unwrap();
            println!("-----------------------------------------------------------------------------------");
            println!("Input: Standard Line ({})", ix);
            println!("{}", input.geodesic);
            println!("Computed: ");
            println!("{}", computed);
            println!();

            let diff_az1_dms = (computed.alpha1.angle() - input.geodesic.alpha1.angle()).to_dms();
            println!("error on az1 (computed - input) = {}", diff_az1_dms);

            let diff_az2_dms = (computed.alpha2.angle() - input.geodesic.alpha2.angle()).to_dms();
            println!("error on az2 (computed - input) = {}", diff_az2_dms);

            let diff_s_m = computed.s - input.geodesic.s;
            println!("error on s (computed - input) = {}", diff_s_m);

            assert_abs_diff_eq!(computed.alpha1, input.geodesic.alpha1, epsilon = 1e-10);
            assert_abs_diff_eq!(computed.alpha2, input.geodesic.alpha2, epsilon = 1e-10);
            assert_abs_diff_eq!(computed.s, input.geodesic.s, epsilon = 1e-3);
        }
    }

    #[test]
    fn inverse_antipodal_lines() {
        // Anti-podal lines
        for (ix, input) in antipodal_lines().iter().enumerate() {
            let solver = KarneyGeodesicSolver::new(&input.ellipsoid);
            let computed = solver
                .solve_inverse(input.geodesic.p1, input.geodesic.p2)
                .unwrap();
            println!("-----------------------------------------------------------------------------------");
            println!("Input: Antipodal Line ({})", ix);
            println!("{}", input.geodesic);
            println!("Computed: ");
            println!("{}", computed);
            println!();

            let diff_az1_dms = (computed.alpha1.angle() - input.geodesic.alpha1.angle()).to_dms();
            println!("error on az1 (computed - input) = {}", diff_az1_dms);

            let diff_az2_dms = (computed.alpha2.angle() - input.geodesic.alpha2.angle()).to_dms();
            println!("error on az2 (computed - input) = {}", diff_az2_dms);

            let diff_s_m = computed.s - input.geodesic.s;
            println!("error on s (computed - input) = {}", diff_s_m);

            assert_abs_diff_eq!(computed.alpha1, input.geodesic.alpha1, epsilon = 1e-10);
            assert_abs_diff_eq!(computed.alpha2, input.geodesic.alpha2, epsilon = 1e-10);
            assert_abs_diff_eq!(computed.s, input.geodesic.s, epsilon = 1e-3);
        }
    }

    #[test]
    fn inverse_equator_lines() {
        let wgs84 = consts::WGS84;
        let solver = KarneyGeodesicSolver::new(&wgs84);
        let computed = solver
            .solve_inverse(
                (Lon::new(-10. * DEG), Lat::new(0. * DEG)),
                (Lon::new(10. * DEG), Lat::new(0. * DEG)),
            )
            .unwrap();
        assert_abs_diff_eq!(computed.alpha1.rad(), PI_2, epsilon = 1e-10);
        assert_abs_diff_eq!(computed.alpha2.rad(), PI_2, epsilon = 1e-10);
        assert_abs_diff_eq!(computed.s, Length::new(2_226_389.816, M), epsilon = 1e-3);

        let computed = solver
            .solve_inverse(
                (Lon::new(10. * DEG), Lat::new(0. * DEG)),
                (Lon::new(-10. * DEG), Lat::new(0. * DEG)),
            )
            .unwrap();
        assert_abs_diff_eq!(computed.alpha1.rad(), -PI_2, epsilon = 1e-10);
        assert_abs_diff_eq!(computed.alpha2.rad(), -PI_2, epsilon = 1e-10);
        assert_abs_diff_eq!(computed.s, Length::new(2_226_389.816, M), epsilon = 1e-3);

        let computed = solver
            .solve_inverse(
                (Lon::new(170. * DEG), Lat::new(0. * DEG)),
                (Lon::new(-170. * DEG), Lat::new(0. * DEG)),
            )
            .unwrap();
        assert_abs_diff_eq!(computed.alpha1.rad(), PI_2, epsilon = 1e-10);
        assert_abs_diff_eq!(computed.alpha2.rad(), PI_2, epsilon = 1e-10);
        assert_abs_diff_eq!(computed.s, Length::new(2_226_389.816, M), epsilon = 1e-3);
    }

    #[test]
    fn inverse_meridian_lines() {
        let wgs84 = consts::WGS84;
        let solver = KarneyGeodesicSolver::new(&wgs84);
        let computed = solver
            .solve_inverse(
                (Lon::new(0. * DEG), Lat::new(-10. * DEG)),
                (Lon::new(0. * DEG), Lat::new(10. * DEG)),
            )
            .unwrap();
        assert_abs_diff_eq!(computed.alpha1.rad(), 0., epsilon = 1e-10);
        assert_abs_diff_eq!(computed.alpha2.rad(), 0., epsilon = 1e-10);
        assert_abs_diff_eq!(computed.s, Length::new(2_211_709.666, M), epsilon = 1e-3);

        let computed = solver
            .solve_inverse(
                (Lon::new(0. * DEG), Lat::new(10. * DEG)),
                (Lon::new(0. * DEG), Lat::new(-10. * DEG)),
            )
            .unwrap();
        assert_abs_diff_eq!(computed.alpha1.rad(), PI, epsilon = 1e-10);
        assert_abs_diff_eq!(computed.alpha2.rad(), PI, epsilon = 1e-10);
        assert_abs_diff_eq!(computed.s, Length::new(2_211_709.666, M), epsilon = 1e-3);

        let computed = solver
            .solve_inverse(
                (Lon::new(0. * DEG), Lat::new(80. * DEG)),
                (Lon::new(180. * DEG), Lat::new(80. * DEG)),
            )
            .unwrap();
        assert_abs_diff_eq!(computed.alpha1.rad(), 0., epsilon = 1e-10);
        assert_abs_diff_eq!(computed.alpha2.rad(), PI, epsilon = 1e-10);
        assert_abs_diff_eq!(computed.s, Length::new(2_233_651.715, M), epsilon = 1e-3);
    }
}
