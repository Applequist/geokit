use crate::cs::geodetic::{Lat, Lon};
use crate::cs::Azimuth;
use crate::geodesy::geodesics::{Geodesic, GeodesicSolver};
use crate::geodesy::Ellipsoid;
use crate::math::complex::Complex;
use crate::math::polynomial::Polynomial;
use crate::quantity::angle::units::RAD;
use std::f64::consts::{FRAC_PI_2, PI};

/// Solve direct and inverse geodesic problems on an ellipsoid using algorithms
/// for Karney - Algorithms for Geodesics.
#[derive(Debug, Clone)]
pub struct KarneyGeodesicSolver<'e> {
    ellipsoid: &'e Ellipsoid,
    a3: Polynomial<6>,
    c3xs: [Polynomial<6>; 5],
}

impl<'e> KarneyGeodesicSolver<'e> {
    /// Create a new solver instance.
    pub fn new(ellipsoid: &'e Ellipsoid) -> Self {
        let n = ellipsoid.n();
        // From Karney - Algorithms for geodesics eqn 24:
        let a3 = Polynomial::new([
            1.,                                                               // x^0
            -Polynomial::new([1. / 2., -1. / 2.]).fast_eval_at(n),            // x^1
            -Polynomial::new([1. / 4., 1. / 8., -3. / 8.]).fast_eval_at(n),   // x^2
            -Polynomial::new([1. / 16., 3. / 16., 1. / 16.]).fast_eval_at(n), // x^3
            -Polynomial::new([3. / 64., 1. / 32.]).fast_eval_at(n),           // x^4
            -3. / 128.,                                                       // x^5
        ]);
        let c3xs = [
            // C31
            Polynomial::new([
                0.,                                                               // x^0
                Polynomial::new([1. / 4., -1. / 4.]).fast_eval_at(n),             // x^1
                Polynomial::new([1. / 8., 0., -1. / 8.]).fast_eval_at(n),         // x^2
                Polynomial::new([3. / 64., 3. / 64., -1. / 64.]).fast_eval_at(n), // x^3
                Polynomial::new([5. / 128., 1. / 64.]).fast_eval_at(n),           // x^4
                3. / 128.,                                                        // x^5
            ]),
            // C32
            Polynomial::new([
                0.,
                0.,
                Polynomial::new([1. / 16., -3. / 32., 1. / 32.]).fast_eval_at(n),
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
        Self {
            ellipsoid,
            a3,
            c3xs,
        }
    }

    pub fn direct(&self, p1: (Lon, Lat), alpha1: Azimuth, s12: f64) -> Geodesic {
        let (lon1, lat1) = p1;
        let beta1 = self.ellipsoid.reduced_latitude(lat1.rad());
        let (sin_beta1, cos_beta1) = beta1.sin_cos();
        let (sin_alpha1, cos_alpha1) = alpha1.rad().sin_cos();

        // Eq (10): alpha0 is the azimuth of the geodesic when crossing the equator.
        let alpha0 = Complex::new(
            Complex::new(cos_alpha1, sin_alpha1 * sin_beta1).abs(),
            sin_alpha1 * cos_beta1,
        )
        .arg();
        let (sin_alpha0, cos_alpha0) = alpha0.sin_cos();

        // Eq (11): sigma1 is the spherical arc from the geodesic on the equator ot p1.
        let sigma1 = Complex::new(cos_alpha1 * cos_beta1, sin_beta1).arg();
        let (sin_sigma1, cos_sigma1) = sigma1.sin_cos();
        // Eq (12): omega1 is the ~longitude~ of p1 from the point of the geodesic on the equator
        let omega1 = Complex::new(cos_sigma1, sin_alpha0 * sin_sigma1).arg();

        // Eq (9)
        let k = self.ellipsoid.e_prime() * cos_alpha0;
        // Eq (16)
        let t = (1. + k * k).sqrt();
        let epsilon = (t - 1.) / (t + 1.);

        // Eq (17)
        let a1 = Self::a1(epsilon);
        // Eq (15) geodesic arc length from point on the geodesic at equator to p1 in meters
        let s1 = self.ellipsoid.b() * Self::i1(a1, epsilon, sigma1);

        // geodesic arc length from point on the geodesic at equator to p2 in meters (p2 is what we are looking for)
        let s2 = s1 + s12;

        let tau2 = s2 / (self.ellipsoid.b() * a1);
        // Eq (11) arc length on the auxiliary sphere from E to p2 in radians
        let sigma2 = Self::sigma(epsilon, tau2);
        let (sin_sigma2, cos_sigma2) = sigma2.sin_cos();

        // Eq (14)
        let alpha2 = Complex::new(cos_alpha0 * cos_sigma2, sin_alpha0).arg();
        // Eq (13)
        let beta2 = Complex::new(
            Complex::new(cos_alpha0 * cos_sigma2, sin_alpha0).abs(),
            cos_alpha0 * sin_sigma2,
        )
        .arg();

        // Eq (12)
        let omega2 = Complex::new(cos_sigma2, sin_alpha0 * sin_sigma2).arg();

        let a3 = self.a3(epsilon);
        // Eq (8)
        let lambda1 = omega1 - self.ellipsoid.f() * sin_alpha0 * self.i3(a3, epsilon, sigma1);
        let lambda2 = omega2 - self.ellipsoid.f() * sin_alpha0 * self.i3(a3, epsilon, sigma2);

        let lambda12 = lambda2 - lambda1;

        // from: tan(beta) = (1 - f)*tan(lat)
        let lat2 = (beta2.tan() / (1. - self.ellipsoid.f())).atan();

        Geodesic {
            p1,
            alpha1,
            p2: (lon1 + lambda12 * RAD, Lat::new(lat2 * RAD)),
            alpha2: Azimuth::new(alpha2 * RAD),
            s: s12,
        }
    }

    pub fn inverse(&self, p1: (Lon, Lat), p2: (Lon, Lat)) -> Geodesic {
        let f = self.ellipsoid.f();
        let a = self.ellipsoid.a();

        // TODO: swap  endpoints and coordinates signs so that lat1 <= 0,
        //       lat1 <= lat2 <= -lat1 and 0 <= lat2 - lat2 <= pi

        let beta1 = self.ellipsoid.reduced_latitude(p1.1.rad());
        println!("beta1 = {}", beta1.to_degrees());
        let beta2 = self.ellipsoid.reduced_latitude(p2.1.rad());
        let (sin_beta1, cos_beta1) = beta1.sin_cos();
        let (sin_beta2, cos_beta2) = beta2.sin_cos();

        let lambda12_0 = (p2.0 - p1.0).length().rad();
        println!("lambda12_0 = {}", lambda12_0.to_degrees());

        // Get initial value of alpha1
        let delta = f * a * PI * cos_beta1.powi(2);
        let x = (lambda12_0 - PI) * a * cos_beta1 / delta;
        println!("x = {}", x);
        let y = (beta2 + beta1) * a / delta;
        println!("y = {}", y);

        let y_sq = y * y;
        let mu = *Polynomial::new([-y_sq, -2. * y_sq, (1. + x * x + y_sq), 2., 1.])
            .solve::<f64>()
            .first()
            .expect("One positive root");
        let mut alpha1 = if y != 0.0 {
            Complex::new(y / mu, -x / (1. + mu)).arg()
        } else {
            // TODO
            0.0f64
        };
        println!("alpha1 = {}", alpha1);
        let mut iter_count = 0;
        loop {
            iter_count += 1;
            let (sin_alpha1, cos_alpha1) = alpha1.sin_cos();

            // Solve for alpha0
            // Eq (5)
            let sin_alpha0 = sin_alpha1 * cos_beta1;
            let cos_alpha0 = (1. - sin_alpha0.powi(2)).sqrt(); // cos_alpha0 > 0
            println!("alpha0 = {}", sin_alpha0.atan2(cos_alpha0).to_degrees());

            // Solve for alpha2
            // Eq (5)
            let sin_alpha2 = sin_alpha0 / cos_beta2;
            // Eq (45)
            let cos_alpha2 = (cos_alpha1.powi(2) * cos_beta1.powi(2)
                + (cos_beta2.powi(2) - cos_beta1.powi(2)))
            .sqrt()
                / cos_beta2;
            let alpha2 = sin_alpha2.atan2(cos_alpha2);
            println!("alpha2 = {}", alpha2.to_degrees());

            // Compute sigma1 and omega1
            // Eq (11)
            let sigma1 = Complex::new(cos_alpha1 * cos_beta1, sin_beta1).arg();
            println!("sigma1 = {}", sigma1.to_degrees());
            let (sin_sigma1, cos_sigma1) = sigma1.sin_cos();
            // Eq (12)
            let omega1 = Complex::new(cos_sigma1, sin_alpha0 * sin_sigma1).arg();
            println!("omega1 = {}", omega1.to_degrees());

            // Compute sigma2 and omega2
            // Eq (11)
            let sigma2 = Complex::new(cos_alpha2 * cos_beta2, sin_beta2).arg();
            println!("sigma2 = {}", sigma2.to_degrees());
            let (sin_sigma2, cos_sigma2) = sigma2.sin_cos();
            // Eq (12)
            let omega2 = Complex::new(cos_sigma2, sin_alpha0 * sin_sigma2).arg();
            println!("omega2 = {}", omega2.to_degrees());

            // Determine lambda12
            // Eq (9)
            let k = self.ellipsoid.e_prime() * cos_alpha0;
            let k_sq = k * k;
            println!("k_sq = {}", k_sq);
            // Eq (16)
            let t = (1. + k_sq).sqrt();
            let epsilon = (t - 1.) / (t + 1.);
            println!("epsilon = {}", epsilon);

            let a3 = self.a3(epsilon);
            // Eq (8)
            let lambda1 = omega1 - self.ellipsoid.f() * sin_alpha0 * self.i3(a3, epsilon, sigma1);
            println!("lambda1 = {}", lambda1.to_degrees());
            let lambda2 = omega2 - self.ellipsoid.f() * sin_alpha0 * self.i3(a3, epsilon, sigma2);
            println!("lambda2 = {}", lambda2.to_degrees());

            let lambda12 = lambda2 - lambda1;
            println!("lambda12 = {}", lambda12);

            let d_lambda12 = lambda12 - lambda12_0;
            println!("d_lambda12 = {}", d_lambda12);
            let d_lambda12_d_alpha1 = if beta2.abs() == beta1.abs() && alpha1 == FRAC_PI_2 {
                // Eq (47): special case
                let sign = (beta2 / beta1).signum();
                -(1. - self.ellipsoid.e_sq() * cos_beta1.powi(2)).sqrt() / sin_beta1
                    * (1. + sign * cos_alpha1)
            } else {
                let a1 = Self::a1(epsilon);
                let a2 = Self::a2(epsilon);
                // Eq (40)
                let j = |sigma| Self::i1(a1, epsilon, sigma) - Self::i2(a2, epsilon, sigma);

                let j_sigma1 = j(sigma1);
                println!("J(sigma1) = {}", j_sigma1);
                let j_sigma2 = j(sigma2);
                println!("J(sigma1) = {}", j_sigma1);
                let m12 = self.ellipsoid.b()
                    * ((1. + k_sq * sin_sigma2.powi(2)).sqrt() * cos_sigma1 * sin_sigma2
                        - (1. + k_sq * sin_sigma1.powi(2)).sqrt() * sin_sigma1 * cos_sigma2
                        - cos_sigma1 * cos_sigma2 * (j_sigma2 - j_sigma1));
                println!("m12 = {}", m12);

                m12 / (self.ellipsoid.a() * cos_alpha2 * cos_beta2)
            };
            println!("d_lambda12_d_alpha1 = {}", d_lambda12_d_alpha1);

            let d_alpha1 = -d_lambda12 / d_lambda12_d_alpha1;
            println!("d_alpha1 = {}", d_alpha1);
            alpha1 += d_alpha1;
            println!("alpha1 = {}", alpha1.to_degrees());

            if iter_count == 2 {
                break;
            }
        }

        Geodesic {
            p1,
            p2,
            ..Default::default()
        }
    }

    /// From Karney - Algorithms for geodesics eqn 17:
    /// A1 = (1 + 1/4 eps^2 + 1/64 eps^4 + 1/256 eps^6 + ...) / (1 - eps)
    fn a1(epsilon: f64) -> f64 {
        Polynomial::new([1., 0., 1. / 4., 0., 1. / 64., 0., 1. / 256.]).fast_eval_at(epsilon)
            / (1. - epsilon)
    }

    /// From Karney - Algorithms for geodesics eqn 15:
    /// I1(sigma) = A1 * (sigma + sum(1, inf, C1l * sin(2l * sigma)))
    fn i1(a1: f64, epsilon: f64, sigma: f64) -> f64 {
        let c1xs = [
            Polynomial::new([0., -0.5, 0., 3. / 16., 0., -1. / 32.0]).fast_eval_at(epsilon),
            Polynomial::new([0., 0., -1. / 16., 0., 1. / 32., 0., -9. / 2048.])
                .fast_eval_at(epsilon),
            Polynomial::new([0., 0., 0., -1. / 48., 0., 3. / 256.]).fast_eval_at(epsilon),
            Polynomial::new([0., 0., 0., 0., -5. / 512., 0., 3. / 512.]).fast_eval_at(epsilon),
            Polynomial::new([0., 0., 0., 0., 0., -7. / 1280.]).fast_eval_at(epsilon),
            Polynomial::new([0., 0., 0., 0., 0., 0., -7. / 2048.]).fast_eval_at(epsilon),
        ];

        a1 * (sigma
            + c1xs
                .into_iter()
                .enumerate()
                .map(|(ix, c1x)| c1x * (2. * sigma * (ix + 1) as f64).sin())
                .sum::<f64>())
    }

    /// From Karney - Algorithms for geodesics eqn 20:
    /// sigma = tau + sum(1, inf, Cp1l * sin(2l * tau)) where tau = s / (b * A1)
    fn sigma(epsilon: f64, tau: f64) -> f64 {
        let cp1xs = [
            Polynomial::new([0., 1. / 2., 0., -9. / 32., 0., 205. / 1536.]).fast_eval_at(epsilon),
            Polynomial::new([0., 0., 5. / 16., 0., -37. / 96., 0., 1335. / 4096.])
                .fast_eval_at(epsilon),
            Polynomial::new([0., 0., 0., 29. / 96., 0., -75. / 128.]).fast_eval_at(epsilon),
            Polynomial::new([0., 0., 0., 0., 539. / 1536., -2391. / 2560.]).fast_eval_at(epsilon),
            Polynomial::new([0., 0., 0., 0., 0., 3467. / 7680.]).fast_eval_at(epsilon),
            Polynomial::new([0., 0., 0., 0., 0., 0., 38081. / 61440.]).fast_eval_at(epsilon),
        ];

        tau + cp1xs
            .into_iter()
            .enumerate()
            .map(|(ix, cp1x)| cp1x * (2. * tau * (ix + 1) as f64).sin())
            .sum::<f64>()
    }

    fn a2(epsilon: f64) -> f64 {
        Polynomial::new([1., 0., 1. / 4., 0., 9. / 64., 0., 25. / 256.]).fast_eval_at(epsilon)
    }

    fn i2(a2: f64, epsilon: f64, sigma: f64) -> f64 {
        let c2xs = [
            Polynomial::new([0., 1. / 2., 0., 1. / 16., 0., 1. / 32.]).fast_eval_at(epsilon),
            Polynomial::new([0., 0., 3. / 16., 0., 1. / 32., 0., 25. / 2048.])
                .fast_eval_at(epsilon),
            Polynomial::new([0., 0., 0., 5. / 48., 0., 5. / 256.]).fast_eval_at(epsilon),
            Polynomial::new([0., 0., 0., 0., 35. / 512., 0., 7. / 512.]).fast_eval_at(epsilon),
            Polynomial::new([0., 0., 0., 0., 0., 63. / 1280.]).fast_eval_at(epsilon),
            Polynomial::new([0., 0., 0., 0., 0., 0., 77. / 2048.]).fast_eval_at(epsilon),
        ];

        a2 * (sigma
            + c2xs
                .into_iter()
                .enumerate()
                .map(|(ix, c2x)| c2x * (2. * sigma * (ix + 1) as f64).sin())
                .sum::<f64>())
    }

    /// From Karney - Algorithms for geodesics eqn 24:
    fn a3(&self, epsilon: f64) -> f64 {
        self.a3.fast_eval_at(epsilon)
    }

    /// From Karney - Algorithms for geodesics eqn 23
    /// I3(sigma) = A3 * (sigma + sum(1, inf, C3l * sin(2l * sigma))
    fn i3(&self, a3: f64, epsilon: f64, sigma: f64) -> f64 {
        let c3xs = self.c3xs.map(|p| p.fast_eval_at(epsilon));
        a3 * (sigma
            + c3xs
                .into_iter()
                .enumerate()
                .map(|(ix, c3x)| c3x * (2. * sigma * (ix + 1) as f64).sin())
                .sum::<f64>())
    }
}

impl<'e> GeodesicSolver for KarneyGeodesicSolver<'e> {
    fn solve_direct(
        &self,
        p1: (Lon, Lat),
        alpha1: Azimuth,
        s12: f64,
    ) -> Result<Geodesic, &'static str> {
        Ok(self.direct(p1, alpha1, s12))
    }

    fn solve_inverse(&self, p1: (Lon, Lat), p2: (Lon, Lat)) -> Result<Geodesic, &'static str> {
        Ok(self.inverse(p1, p2))
    }
}

#[cfg(test)]
mod tests {
    use crate::cs::geodetic::{Lat, Lon};
    use crate::cs::Azimuth;
    use crate::geodesy::ellipsoid::consts;
    use crate::geodesy::geodesics::karney::KarneyGeodesicSolver;
    use crate::geodesy::geodesics::tests::{
        antipodal_lines, standard_lines, vincenty_direct_deltas, vincenty_inverse_deltas,
        vincenty_lines,
    };
    use crate::geodesy::geodesics::GeodesicSolver;
    use crate::quantity::angle::units::{DEG, RAD};
    use approx::assert_abs_diff_eq;
    use std::f64::consts::{FRAC_PI_2, PI};

    #[test]
    fn direct_karney_sample() {
        let wgs84 = consts::WGS84;
        let solver = KarneyGeodesicSolver::new(&wgs84);
        let computed = solver
            .solve_direct(
                (Lon::new(0. * DEG), Lat::new(40. * DEG)),
                Azimuth::new(30. * DEG),
                10_000_000.,
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

            assert_abs_diff_eq!(computed.p2.0, input.geodesic.p2.0, epsilon = 1e-10);
            assert_abs_diff_eq!(computed.p2.1, input.geodesic.p2.1, epsilon = 1e-10);
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

            assert_abs_diff_eq!(computed.p2.0, input.geodesic.p2.0, epsilon = 1e-10);
            assert_abs_diff_eq!(computed.p2.1, input.geodesic.p2.1, epsilon = 1e-10);
            assert_abs_diff_eq!(computed.alpha2, input.geodesic.alpha2, epsilon = 1e-10);
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
                20_000.0,
            )
            .unwrap();
        assert_abs_diff_eq!(
            computed.p2.0.rad(),
            (0.17966306 * DEG).rad(),
            epsilon = 1e-8
        );
        assert_abs_diff_eq!(computed.p2.1.rad(), 0.0, epsilon = 1e-8);
        assert_abs_diff_eq!(computed.alpha2.rad(), FRAC_PI_2, epsilon = 1e-8);

        let computed = solver
            .solve_direct(
                (Lon::new(170.0 * DEG), Lat::new(0.0 * DEG)),
                Azimuth::new(90.0 * DEG),
                2_000_000.0,
            )
            .unwrap();
        assert_abs_diff_eq!(
            computed.p2.0,
            Lon::new(-172.03369432 * DEG),
            epsilon = 1e-10
        );
        assert_abs_diff_eq!(computed.p2.1, Lat::new(0.0 * RAD), epsilon = 1e-10);
        assert_abs_diff_eq!(
            computed.alpha2,
            Azimuth::new(FRAC_PI_2 * RAD),
            epsilon = 1e-8
        );
    }

    #[test]
    fn direct_meridian_lines() {
        let wgs84 = consts::WGS84;
        let solver = KarneyGeodesicSolver::new(&wgs84);

        let computed = solver
            .solve_direct(
                (Lon::new(0. * DEG), Lat::new(-10. * DEG)),
                Azimuth::new(0. * DEG),
                2_000_000.0,
            )
            .unwrap();
        assert_abs_diff_eq!(computed.p2.0, Lon::new(0.0 * RAD), epsilon = 1e-10);
        assert_abs_diff_eq!(computed.p2.1, Lat::new(8.08583903 * DEG), epsilon = 1e-10);
        assert_abs_diff_eq!(computed.alpha2, Azimuth::new(0.0 * RAD), epsilon = 1e-8);

        let computed = solver
            .solve_direct(
                (Lon::new(0. * DEG), Lat::new(80. * DEG)),
                Azimuth::new(0. * DEG),
                2_000_000.0,
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
        assert_abs_diff_eq!(computed.s, 4.944_208, epsilon = 1e-6);
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
        assert_abs_diff_eq!(computed.s, 19_989_832.82761, epsilon = 1e-6);
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

            assert!(diff_s_m < delta.delta_s.abs());
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
        assert_abs_diff_eq!(computed.alpha1.rad(), FRAC_PI_2, epsilon = 1e-10);
        assert_abs_diff_eq!(computed.alpha2.rad(), FRAC_PI_2, epsilon = 1e-10);
        assert_abs_diff_eq!(computed.s, 2_226_389.816, epsilon = 1e-3);

        let computed = solver
            .solve_inverse(
                (Lon::new(10. * DEG), Lat::new(0. * DEG)),
                (Lon::new(-10. * DEG), Lat::new(0. * DEG)),
            )
            .unwrap();
        assert_abs_diff_eq!(computed.alpha1.rad(), -FRAC_PI_2, epsilon = 1e-10);
        assert_abs_diff_eq!(computed.alpha2.rad(), -FRAC_PI_2, epsilon = 1e-10);
        assert_abs_diff_eq!(computed.s, 2_226_389.816, epsilon = 1e-3);

        let computed = solver
            .solve_inverse(
                (Lon::new(170. * DEG), Lat::new(0. * DEG)),
                (Lon::new(-170. * DEG), Lat::new(0. * DEG)),
            )
            .unwrap();
        assert_abs_diff_eq!(computed.alpha1.rad(), FRAC_PI_2, epsilon = 1e-10);
        assert_abs_diff_eq!(computed.alpha2.rad(), FRAC_PI_2, epsilon = 1e-10);
        assert_abs_diff_eq!(computed.s, 2_226_389.816, epsilon = 1e-3);
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
        assert_abs_diff_eq!(computed.s, 2_211_709.666, epsilon = 1e-3);

        let computed = solver
            .solve_inverse(
                (Lon::new(0. * DEG), Lat::new(10. * DEG)),
                (Lon::new(0. * DEG), Lat::new(-10. * DEG)),
            )
            .unwrap();
        assert_abs_diff_eq!(computed.alpha1.rad(), PI, epsilon = 1e-10);
        assert_abs_diff_eq!(computed.alpha2.rad(), PI, epsilon = 1e-10);
        assert_abs_diff_eq!(computed.s, 2_211_709.666, epsilon = 1e-3);

        let computed = solver
            .solve_inverse(
                (Lon::new(0. * DEG), Lat::new(80. * DEG)),
                (Lon::new(180. * DEG), Lat::new(80. * DEG)),
            )
            .unwrap();
        assert_abs_diff_eq!(computed.alpha1.rad(), 0., epsilon = 1e-10);
        assert_abs_diff_eq!(computed.alpha2.rad(), PI, epsilon = 1e-10);
        assert_abs_diff_eq!(computed.s, 2_233_651.715, epsilon = 1e-3);
    }
}
