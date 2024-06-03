use crate::cs::Azimuth;
use crate::cs::geodetic::{Lat, Lon};
use crate::geodesy::Ellipsoid;
use crate::geodesy::geodesics::{Geodesic, GeodesicSolver};
use crate::math::polynomial::Polynomial;

pub struct RappIterativeGeodisicSolver<'e> {
    ellipsoid: &'e Ellipsoid,
    axs: [Polynomial<4>; 4],
}

impl<'e> RappIterativeGeodisicSolver<'e> {
    pub fn new(ellipsoid: &'e Ellipsoid) -> Self {
        // From eq 1.56 to eval at cos_alpha_sq
        let f = ellipsoid.f();
        let axs = [
            // A0
            Polynomial::new([
                1.,
                -Polynomial::new([0., 1. / 4., 1. / 4., 1. / 4.]).eval_at(f),
                Polynomial::new([0., 0., 3. / 16., 27. / 64.]).eval_at(f),
                -Polynomial::new([0., 0., 0., 25. / 128.]).eval_at(f),
            ]),
            // A2
            Polynomial::new([
                0.,
                Polynomial::new([0., 1. / 4., 1. / 4., 1. / 4.]).eval_at(f),
                -Polynomial::new([0., 0., 1. / 4., 9. / 16.]).eval_at(f),
                Polynomial::new([0., 0., 0., 75. / 256.]).eval_at(f),
            ]),
            // A4
            Polynomial::new([
                0.,
                0.,
                Polynomial::new([0., 0., 1. / 32., 9. / 128.]).eval_at(f),
                -Polynomial::new([0., 0., 0., 15. / 256.]).eval_at(f),
            ]),
            // A6
            Polynomial::new([
                0.,
                0.,
                0.,
                Polynomial::new([0., 0., 0., 5. / 768.]).eval_at(f),
            ])
        ];
        Self {
            ellipsoid,
            axs,
        }
    }
}

impl<'e> GeodesicSolver for RappIterativeGeodisicSolver<'e> {
    fn solve_direct(&self, p1: (Lon, Lat), alpha1: Azimuth, s12: f64) -> Geodesic {
        todo!()
    }

    fn solve_inverse(&self, p1: (Lon, Lat), p2: (Lon, Lat)) -> Geodesic {
        // Compute the difference between lambda (delta lon on the sphere) and L (delta lon on
        // the ellipsoid
        fn step(f: f64, axs: &[Polynomial<4>], (sin_beta1, cos_beta1): (f64, f64), (sin_beta2, cos_beta2): (f64, f64), lambda: f64) -> (f64, f64, f64, [f64; 4], [f64; 4]) {
            let (sin_lambda, cos_lambda) = lambda.sin_cos();

            // Compute sigma
            let cos_sigma = sin_beta1 * sin_beta2 + cos_beta1 * cos_beta2 * cos_lambda; // Eq 1.57
            let sin_sigma = ((sin_lambda * cos_beta2).powi(2) + (sin_beta2 * cos_beta1 - sin_beta1 * cos_beta2 * cos_lambda).powi(2)).sqrt(); // Eq 1.58
            let sigma = sin_sigma.atan2(cos_sigma);
            let sin_p_sigmas = [
                sin_sigma,
                (2. * sigma).sin(),
                (3. * sigma).sin(),
                (4. * sigma).sin(),
            ];

            // Compute alpha
            let sin_alpha = sin_lambda * cos_beta1 * cos_beta2 / sin_sigma;
            let alpha = sin_alpha.asin();
            let cos_alpha_sq = alpha.cos().powi(2);

            // Compute cos 2p sigma_m
            let cos_2sigma_m = cos_sigma - 2. * sin_beta1 * sin_beta2 / cos_alpha_sq;
            let cos_4sigma_m = 2. * cos_2sigma_m.powi(2) - 1.;
            let cos_6sigma_m = cos_2sigma_m.powi(3) - 3. * cos_2sigma_m * (1. - cos_2sigma_m.powi(2));
            let cos_8sigma_m = cos_2sigma_m.powi(4) - 6. * cos_2sigma_m.powi(2) * (1. - cos_2sigma_m.powi(2)) + (1. - cos_2sigma_m.powi(2)).powi(2);
            let cos_2p_sigma_ms = [
                cos_2sigma_m,
                cos_4sigma_m,
                cos_6sigma_m,
                cos_8sigma_m,
            ];

            let a0 = axs[0].eval_at(cos_alpha_sq);
            let a2 = axs[1].eval_at(cos_alpha_sq);
            let a4 = axs[2].eval_at(cos_alpha_sq);
            let a6 = axs[3].eval_at(cos_alpha_sq);
            (
                f * sin_alpha * (
                    a0 * sigma +
                        a2 * sin_p_sigmas[0] * cos_2sigma_m +
                        a4 * sin_p_sigmas[1] * cos_4sigma_m +
                        a6 * sin_p_sigmas[2] * cos_6sigma_m),
                alpha,
                sigma,
                sin_p_sigmas,
                cos_2p_sigma_ms
            )
        }

        fn distance(u_sq: f64, b: f64, sigma: f64, sin_p_sigmas: [f64; 4], cos_2p_sigma_ms: [f64; 4]) -> f64 {
            let b0 = Polynomial::new([1., 1. / 4., -3. / 64., 5. / 256., -175. / 16384.]).eval_at(u_sq);
            let b2 = Polynomial::new([0., -1. / 4., 1. / 16., -15. / 512., 35. / 2048.]).eval_at(u_sq);
            let b4 = Polynomial::new([0., 0., -1. / 128., 3. / 512., -35. / 8192.]).eval_at(u_sq);
            let b6 = Polynomial::new([0., 0., 0., -1. / 1536., 5. / 6144.]).eval_at(u_sq);
            let b8 = Polynomial::new([0., 0., 0., 0., -5. / 65536.]).eval_at(u_sq);

            b * (
                b0 * sigma +
                    b2 * sin_p_sigmas[0] * cos_2p_sigma_ms[0] +
                    b4 * sin_p_sigmas[1] * cos_2p_sigma_ms[1] +
                    b6 * sin_p_sigmas[2] * cos_2p_sigma_ms[2] +
                    b8 * sin_p_sigmas[3] * cos_2p_sigma_ms[3])
        }

        let (lon1, lat1) = p1;
        let (lon2, lat2) = p2;
        let beta1 = self.ellipsoid.reduced_latitude(lat1.rad());
        let (sin_beta1, cos_beta1) = beta1.sin_cos();
        let beta2 = self.ellipsoid.reduced_latitude(lat2.rad());
        let (sin_beta2, cos_beta2) = beta2.sin_cos();
        let L = (lon2 - lon1).length().unwrap_or(0.);

        // first approximation of lambda = L = lon2 - lon1
        let mut lambda = L;
        let mut s = 0.;
        let mut az1 = 0.;
        let mut az2 = 0.;
        let mut delta_lambda = 0.;
        let mut iter_count = 0;
        loop {
            let (delta_lambda_new, alpha, sigma, sin_p_sigmas, cos_2p_sigma_ms) = step(self.ellipsoid.f(), &self.axs, (sin_beta1, cos_beta1), (sin_beta2, cos_beta2), lambda);
            let sum = delta_lambda + delta_lambda_new;
            iter_count += 1;

            if (delta_lambda - delta_lambda_new).abs() < 0.5e-14 || sum.abs() < 3e-12 {
                s = distance(self.ellipsoid.e_prime_sq() * alpha.cos().powi(2), self.ellipsoid.b(), sigma, sin_p_sigmas, cos_2p_sigma_ms);
                let sin_alpha1 = alpha.sin() / beta1.cos();
                let sin_alpha2 = alpha.sin() / beta2.cos();
                let tan_alpha12 = if s < 1000.0 {
                    lambda.sin() * cos_beta2 / ((beta2 - beta1).sin() + 2. * sin_beta1 * cos_beta2 * (lambda / 2.).sin().powi(2))
                } else {
                    lambda.sin() * cos_beta2 / (sin_beta2 * cos_beta1 - lambda.cos() * sin_beta1 * cos_beta2)
                };
                let tan_alpha21 = if s < 1000.0 {
                    lambda.sin() * cos_beta1 / ((beta2 - beta1).sin() - 2. * cos_beta1 * sin_beta2 * (lambda / 2.).sin().powi(2))
                } else {
                    lambda.sin() * cos_beta1 / (sin_beta2 * cos_beta1 * lambda.cos() - sin_beta1 * cos_beta2)
                };
                az1 = sin_alpha1.atan2(sin_alpha1 / tan_alpha12);
                az2 = sin_alpha2.atan2(sin_alpha2 / tan_alpha21);
                break;
            }
            delta_lambda = delta_lambda_new;
            lambda = L + delta_lambda;
        }

        println!("Iteration count = {}", iter_count);
        Geodesic {
            p1: (lon1, lat1),
            alpha1: Azimuth::new(az1),
            p2: (lon2, lat2),
            alpha2: Azimuth::new(az2),
            s,
        }
    }

}

#[cfg(test)]
mod tests {
    use approx::assert_abs_diff_eq;
    use crate::cs::geodetic::Lon;
    use crate::geodesy::ellipsoid::consts;
    use crate::geodesy::geodesics::GeodesicSolver;
    use crate::geodesy::geodesics::rapp::RappIterativeGeodisicSolver;
    use crate::geodesy::geodesics::tests::{antipodal_lines, standard_lines};

    #[test]
    fn solve_direct() {
        let intl = consts::INTL;
        let solver = RappIterativeGeodisicSolver::new(&intl);

        // Direct problem
        for (ix, g) in standard_lines().iter().enumerate() {
            println!("Standard Test Line {}", ix + 1);
            let res = solver.solve_direct(g.p1, g.alpha1, g.s);
            let (lon, lat) = res.p2;
            let (expected_lon, expected_lat) = g.p2;
            assert_abs_diff_eq!(lon, expected_lon, epsilon = 2e-9);
            assert_abs_diff_eq!(lat, expected_lat, epsilon = 2e-9);
            assert_abs_diff_eq!(res.alpha2, g.alpha2, epsilon = 1e-7);
        }

        for (ix, g) in antipodal_lines().iter().enumerate() {
            println!("Anti-podal Test Line {}", ix + 1);
            let res = solver.solve_direct(g.p1, g.alpha1, g.s);
            let (lon, lat) = res.p2;
            let (expected_lon, expected_lat) = g.p2;
            assert_abs_diff_eq!(lon, expected_lon, epsilon = 2e-9);
            assert_abs_diff_eq!(lat, expected_lat, epsilon = 2e-9);
            assert_abs_diff_eq!(res.alpha2, g.alpha2, epsilon = 1e-7);
        }
    }
    #[test]
    fn solve_inverse() {
        let intl = consts::INTL;
        let solver = RappIterativeGeodisicSolver::new(&intl);

        for (ix, g) in standard_lines().iter().enumerate() {
            println!("Standard Test Line {}", ix + 1);
            let res = solver.solve_inverse(g.p1, g.p2);
            assert_abs_diff_eq!(res.alpha1, g.alpha1, epsilon = 2e-10);
            assert_abs_diff_eq!(res.alpha2, g.alpha2, epsilon = 1e-10);
            assert_abs_diff_eq!(res.s, g.s, epsilon = 1e-3);
        }

        for (ix, g) in antipodal_lines().iter().enumerate() {
            println!("Anti-podal Test Line {}", ix + 1);
            let res = solver.solve_inverse(g.p1, g.p2);
            assert_abs_diff_eq!(res.alpha1, g.alpha1, epsilon = 2e-10);
            assert_abs_diff_eq!(res.alpha2, g.alpha2, epsilon = 1e-10);
            assert_abs_diff_eq!(res.s, g.s, epsilon = 1e-3);
        }

    }
}