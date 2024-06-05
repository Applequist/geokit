use crate::cs::geodetic::{Lat, Lon};
use crate::cs::Azimuth;
use crate::geodesy::geodesics::{Geodesic, GeodesicSolver};
use crate::geodesy::Ellipsoid;
use crate::math::complex::Complex;
use crate::math::polynomial::Polynomial;

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
            1.,                                                          // x^0
            -Polynomial::new([1. / 2., -1. / 2.]).eval_at(n),            // x^1
            -Polynomial::new([1. / 4., 1. / 8., -3. / 8.]).eval_at(n),   // x^2
            -Polynomial::new([1. / 16., 3. / 16., 1. / 16.]).eval_at(n), // x^3
            -Polynomial::new([3. / 64., 1. / 32.]).eval_at(n),           // x^4
            -3. / 128.,                                                  // x^5
        ]);
        let c3xs = [
            // C31
            Polynomial::new([
                0.,                                                          // x^0
                Polynomial::new([1. / 4., -1. / 4.]).eval_at(n),             // x^1
                Polynomial::new([1. / 8., 0., -1. / 8.]).eval_at(n),         // x^2
                Polynomial::new([3. / 64., 3. / 64., -1. / 64.]).eval_at(n), // x^3
                Polynomial::new([5. / 128., 1. / 64.]).eval_at(n),           // x^4
                3. / 128.,                                                   // x^5
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

    pub fn direct(&self, p1: (Lon, Lat), alpha1: Azimuth, s12: f64) -> Geodesic {
        // Step 1: solve NEP_1 to give alpha0, sigma1 and omega1
        let (lon1, lat1) = p1;
        let beta1 = self.ellipsoid.reduced_latitude(lat1.rad());
        let (sin_beta1, cos_beta1) = beta1.sin_cos();
        let (sin_alpha1, cos_alpha1) = alpha1.rad().sin_cos();
        // alpha0 is the azimuth at the equator of the spherical triangle NEP_1, N: north pole,
        // E: intersection of the geodesic and the equator
        // Karney - Algorithms for geodesics eqn 10
        let alpha0 = Complex::new(
            Complex::new(cos_alpha1, sin_alpha1 * sin_beta1).abs(),
            sin_alpha1 * cos_beta1,
        )
            .theta();
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
        let beta2 = Complex::new(
            Complex::new(cos_alpha0 * cos_sigma2, sin_alpha0).abs(),
            cos_alpha0 * sin_sigma2,
        )
            .theta();
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

        Geodesic {
            p1,
            alpha1,
            p2: (lon1 + lambda12, Lat::new(lat2)),
            alpha2: Azimuth::new(alpha2),
            s: s12,
        }
    }

    pub fn inverse(&self, p1: (Lon, Lat), p2: (Lon, Lat)) -> Geodesic {
        // TODO:
        Geodesic {
            p1,
            p2,
            ..Default::default()
        }
    }

    /// From Karney - Algorithms for geodesics eqn 17:
    /// A1 = (1 + 1/4 eps^2 + 1/64 eps^4 + 1/256 eps^6 + ...) / (1 - eps)
    fn a1(epsilon: f64) -> f64 {
        Polynomial::new([1., 0., 1. / 4., 0., 1. / 64., 0., 1. / 256.]).eval_at(epsilon)
            / (1. - epsilon)
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
            Polynomial::new([0., 1. / 2., 0., -9. / 32., 0., 205. / 1536.]).eval_at(epsilon),
            Polynomial::new([0., 0., 5. / 16., 0., -37. / 96., 0., 1335. / 4096.]).eval_at(epsilon),
            Polynomial::new([0., 0., 0., 29. / 96., 0., -75. / 128.]).eval_at(epsilon),
            Polynomial::new([0., 0., 0., 0., 539. / 1536., -2391. / 2560.]).eval_at(epsilon),
            Polynomial::new([0., 0., 0., 0., 0., 3467. / 7680.]).eval_at(epsilon),
            Polynomial::new([0., 0., 0., 0., 0., 0., 38081. / 61440.]).eval_at(epsilon),
        ];

        tau + cp1xs
            .into_iter()
            .enumerate()
            .map(|(ix, cp1x)| cp1x * (2. * tau * (ix + 1) as f64).sin())
            .sum::<f64>()
    }

    fn a2(epsilon: f64) -> f64 {
        Polynomial::new([1., 0., 1. / 4., 0., 9. / 64., 0., 25. / 256.]).eval_at(epsilon)
    }

    fn i2(a2: f64, epsilon: f64, sigma: f64) -> f64 {
        let c2xs = [
            Polynomial::new([0., 1. / 2., 0., 1. / 16., 0., 1. / 32.]).eval_at(epsilon),
            Polynomial::new([0., 0., 3. / 16., 0., 1. / 32., 0., 25. / 2048.]).eval_at(epsilon),
            Polynomial::new([0., 0., 0., 5. / 48., 0., 5. / 256.]).eval_at(epsilon),
            Polynomial::new([0., 0., 0., 0., 35. / 512., 0., 7. / 512.]).eval_at(epsilon),
            Polynomial::new([0., 0., 0., 0., 0., 63. / 1280.]).eval_at(epsilon),
            Polynomial::new([0., 0., 0., 0., 0., 0., 77. / 2048.]).eval_at(epsilon),
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
        self.a3.eval_at(epsilon)
    }

    /// From Karney - Algorithms for geodesics eqn 23
    /// I3(sigma) = A3 * (sigma + sum(1, inf, C3l * sin(2l * sigma))
    fn i3(&self, a3: f64, epsilon: f64, sigma: f64) -> f64 {
        let c3xs = self.c3xs.map(|p| p.eval_at(epsilon));
        a3 * (sigma
            + c3xs
                .into_iter()
                .enumerate()
                .map(|(ix, c3x)| c3x * (2. * sigma * (ix + 1) as f64).sin())
                .sum::<f64>())
    }
}

impl<'e> GeodesicSolver for KarneyGeodesicSolver<'e> {
    fn solve_direct(&self, p1: (Lon, Lat), alpha1: Azimuth, s12: f64) -> Geodesic {
        self.direct(p1, alpha1, s12)
    }

    fn solve_inverse(&self, p1: (Lon, Lat), p2: (Lon, Lat)) -> Geodesic {
        self.inverse(p1, p2)
    }
}

#[cfg(test)]
mod tests {
    use approx::assert_abs_diff_eq;
    use crate::geodesy::geodesics::GeodesicSolver;
    use crate::geodesy::geodesics::karney::KarneyGeodesicSolver;
    use crate::geodesy::geodesics::tests::{antipodal_lines, standard_lines, vincenty_direct_deltas, vincenty_inverse_deltas, vincenty_lines};
    use crate::quantity::angle::DMS;

    #[test]
    fn direct_vincenty_lines() {
        let chars = ['a', 'b', 'c', 'd', 'e', 'f'];
        for (ix, (input, delta)) in vincenty_lines()
            .into_iter()
            .zip(vincenty_direct_deltas())
            .enumerate()
        {
            let solver = KarneyGeodesicSolver::new(&input.ellipsoid);
            let computed =
                solver.solve_direct(input.geodesic.p1, input.geodesic.alpha1, input.geodesic.s);
            println!("-----------------------------------------------------------------------------------");
            println!("Input: Vincenty Line ({})", chars[ix]);
            println!("{}", input.geodesic);
            println!("Computed: ");
            println!("{}", computed);
            println!();

            let diff_lon_dms = DMS::from_rad(computed.p2.0.rad() - input.geodesic.p2.0.rad());
            println!("error on lon (computed - input) = {}", diff_lon_dms);

            let diff_lat_dms = DMS::from_rad(computed.p2.1.rad() - input.geodesic.p2.1.rad());
            println!("error on lat (computed - input) = {}", diff_lat_dms);

            let diff_az_dms = DMS::from_rad(computed.alpha2.rad() - input.geodesic.alpha2.rad());
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
            let computed = solver.solve_direct(input.geodesic.p1, input.geodesic.alpha1, input.geodesic.s);
            println!("-----------------------------------------------------------------------------------");
            println!("Input: Standard Line ({})", ix);
            println!("{}", input.geodesic);
            println!("Computed: ");
            println!("{}", computed);
            println!();

            let diff_lon_dms = DMS::from_rad(computed.p2.0.rad() - input.geodesic.p2.0.rad());
            println!("error on lon (computed - input) = {}", diff_lon_dms);

            let diff_lat_dms = DMS::from_rad(computed.p2.1.rad() - input.geodesic.p2.1.rad());
            println!("error on lat (computed - input) = {}", diff_lat_dms);

            let diff_az_dms = DMS::from_rad(computed.alpha2.rad() - input.geodesic.alpha2.rad());
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
            let computed = solver.solve_direct(input.geodesic.p1, input.geodesic.alpha1, input.geodesic.s);
            println!("-----------------------------------------------------------------------------------");
            println!("Input: Antipodal Line ({})", ix);
            println!("{}", input.geodesic);
            println!("Computed: ");
            println!("{}", computed);
            println!();

            let diff_lon_dms = DMS::from_rad(computed.p2.0.rad() - input.geodesic.p2.0.rad());
            println!("error on lon (computed - input) = {}", diff_lon_dms);

            let diff_lat_dms = DMS::from_rad(computed.p2.1.rad() - input.geodesic.p2.1.rad());
            println!("error on lat (computed - input) = {}", diff_lat_dms);

            let diff_az_dms = DMS::from_rad(computed.alpha2.rad() - input.geodesic.alpha2.rad());
            println!("error on az (computed - input) = {}", diff_az_dms);

            assert_abs_diff_eq!(computed.p2.0, input.geodesic.p2.0, epsilon = 1e-10);
            assert_abs_diff_eq!(computed.p2.1, input.geodesic.p2.1, epsilon = 1e-10);
            assert_abs_diff_eq!(computed.alpha2, input.geodesic.alpha2, epsilon = 1e-10);
        }
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
            let computed = solver.solve_inverse(input.geodesic.p1, input.geodesic.p2);
            println!("-----------------------------------------------------------------------------------");
            println!("Input: Vincenty Line ({})", chars[ix]);
            println!("{}", input.geodesic);
            println!("Computed: ");
            println!("{}", computed);
            println!();

            let diff_az1_dms = DMS::from_rad(computed.alpha1.rad() - input.geodesic.alpha1.rad());
            println!("error on az1 (computed - input) = {}", diff_az1_dms);

            let diff_az2_dms = DMS::from_rad(computed.alpha2.rad() - input.geodesic.alpha2.rad());
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
            let computed = solver.solve_inverse(input.geodesic.p1, input.geodesic.p2);
            println!("-----------------------------------------------------------------------------------");
            println!("Input: Standard Line ({})", ix);
            println!("{}", input.geodesic);
            println!("Computed: ");
            println!("{}", computed);
            println!();

            let diff_az1_dms = DMS::from_rad(computed.alpha1.rad() - input.geodesic.alpha1.rad());
            println!("error on az1 (computed - input) = {}", diff_az1_dms);

            let diff_az2_dms = DMS::from_rad(computed.alpha2.rad() - input.geodesic.alpha2.rad());
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
            let computed = solver.solve_inverse(input.geodesic.p1, input.geodesic.p2);
            println!("-----------------------------------------------------------------------------------");
            println!("Input: Antipodal Line ({})", ix);
            println!("{}", input.geodesic);
            println!("Computed: ");
            println!("{}", computed);
            println!();

            let diff_az1_dms = DMS::from_rad(computed.alpha1.rad() - input.geodesic.alpha1.rad());
            println!("error on az1 (computed - input) = {}", diff_az1_dms);

            let diff_az2_dms = DMS::from_rad(computed.alpha2.rad() - input.geodesic.alpha2.rad());
            println!("error on az2 (computed - input) = {}", diff_az2_dms);

            let diff_s_m = computed.s - input.geodesic.s;
            println!("error on s (computed - input) = {}", diff_s_m);

            assert_abs_diff_eq!(computed.alpha1, input.geodesic.alpha1, epsilon = 1e-10);
            assert_abs_diff_eq!(computed.alpha2, input.geodesic.alpha2, epsilon = 1e-10);
            assert_abs_diff_eq!(computed.s, input.geodesic.s, epsilon = 1e-3);
        }
    }
}
