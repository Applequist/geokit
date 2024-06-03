use crate::cs::geodetic::{Lat, Lon};
use crate::cs::Azimuth;
use crate::geodesy::geodesics::{Geodesic, GeodesicSolver};
use crate::geodesy::Ellipsoid;
use crate::math::polynomial::Polynomial;
use crate::math::utils::iter_fn;

pub struct VincentyGeodesicSolver<'e> {
    ellipsoid: &'e Ellipsoid,
}

impl<'e> VincentyGeodesicSolver<'e> {
    pub fn new(ellipsoid: &'e Ellipsoid) -> Self {
        Self { ellipsoid }
    }
}

impl<'e> GeodesicSolver for VincentyGeodesicSolver<'e> {
    fn solve_direct(&self, p1: (Lon, Lat), alpha1: Azimuth, s12: f64) -> Geodesic {
        let (lon1, lat1)  = p1;
        let f = self.ellipsoid.f();
        let beta1 = self.ellipsoid.reduced_latitude(lat1.rad());
        let (sin_alpha1, cos_alpha1) = alpha1.rad().sin_cos();

        // Eq (1): indeterminate for equatorial lines !!!
        let tan_sigma1 = beta1.tan() / cos_alpha1;
        let sigma1 = tan_sigma1.atan();
        // Eq (2)
        let sin_alpha = beta1.cos() / sin_alpha1;
        let cos_alpha_sq = 1. - sin_alpha.powi(2);
        let u_sq = cos_alpha_sq * self.ellipsoid.e_prime_sq();

        // Eq (3) & (4)
        let A = Polynomial::new([
            1.,
            4096. / 16384.,
            -768. / 16384.,
            320. / 16384.,
            -175. / 16384.,
        ])
        .eval_at(u_sq);
        let B = Polynomial::new([0., 256. / 1024., -128. / 1024., 74. / 1024., -47. / 1024.])
            .eval_at(u_sq);

        let delta_sigma = move |sigma: [f64; 1]| -> [f64; 1] {
            let (sin_sigma, cos_sigma) = sigma[0].sin_cos();
            // Eq (5)
            let two_sigma_m = 2. * sigma1 + sigma[0];

            let cos_2_sigma_m = two_sigma_m.cos();

            // Eq (6)
            [B * sin_sigma
                * (cos_2_sigma_m
                    + 0.25
                        * B
                        * (cos_sigma * (-1. + 2. * cos_2_sigma_m.powi(2))
                            - (1. / 6.)
                                * B
                                * cos_2_sigma_m
                                * (-3. + 4. * sin_sigma.powi(2))
                                * (-3. + 4. * cos_2_sigma_m.powi(2))))]
        };

        let [sigma] = iter_fn([s12 / (self.ellipsoid.b() * A)], &delta_sigma, [0.5e-14]);

        let (sin_beta1, cos_beta1) = beta1.sin_cos();
        let (sin_sigma, cos_sigma) = sigma.sin_cos();
        let lat2 = (sin_beta1 * cos_sigma + cos_beta1 * sin_sigma * cos_alpha1).atan2(
            (1. - f)
                * (sin_alpha.powi(2)
                    + (sin_beta1 * sin_sigma - cos_beta1 * cos_sigma * cos_alpha1).powi(2))
                .sqrt(),
        );
        let lambda = (sin_sigma * sin_alpha1)
            .atan2(cos_beta1 * cos_sigma - sin_beta1 * sin_sigma * cos_alpha1);

        let C = Polynomial::new([
            0.,
            Polynomial::new([0., 1. / 4., 1. / 4.]).eval_at(f),
            Polynomial::new([0., 0., -3. / 16.]).eval_at(f),
        ])
        .eval_at(cos_alpha_sq);

        let cos_2_sigma_m = (2. * sigma1 + sigma).cos();
        let L = lambda
            + (1. - C)
                * f
                * sin_alpha
                * (sigma
                    + C * sin_sigma
                        * (cos_2_sigma_m + C * cos_sigma * (-1. + 2. * cos_2_sigma_m.powi(2))));

        let alpha2 = sin_alpha.atan2(-sin_beta1 * sin_sigma + cos_beta1 * cos_sigma * cos_alpha1);

        Geodesic {
            p1: (lon1, lat1), alpha1,
            p2: (lon1 + L, Lat::new(lat2)), alpha2: Azimuth::new(alpha2),
            s: s12,
        }
    }

    fn solve_inverse(&self, p1: (Lon, Lat), p2: (Lon, Lat)) -> Geodesic {
        todo!()
    }
}

#[cfg(test)]
mod tests {
    use crate::cs::geodetic::{Lat, Lon};
    use crate::cs::Azimuth;
    use crate::geodesy::ellipsoid::consts;
    use crate::geodesy::geodesics::vincenty::VincentyGeodesicSolver;
    use crate::geodesy::geodesics::{Geodesic, GeodesicSolver};
    use crate::geodesy::Ellipsoid;
    use crate::operation::Inv;
    use crate::quantity::angle::dms;
    use crate::quantity::angle::units::DEG;
    use approx::assert_abs_diff_eq;

    struct DirectDeltas {
        delta_lat2: f64,
        delta_delta_lon: f64,
        delta_alpha2: f64,
    }

    struct InverseDeltas {
        delta_alpha1: f64,
        delta_alpha2: f64,
        delta_s: f64,
    }

    pub struct LineData {
        ellipsoid: Ellipsoid,
        geodesic: Geodesic,
        direct_deltas: DirectDeltas,
        inverse_deltas: InverseDeltas,
    }

    /// From Vincenty - Direct and inverse solutions of geodesics on the ellipsoid with
    ///                 application of nested equations
    /// Table II - Results of solutions
    fn lines() -> Vec<LineData> {
        vec![
            // Line (a)
            LineData {
                ellipsoid: consts::BESSEL,
                geodesic: Geodesic {
                    p1: (Lon::new(0.), Lat::new(dms(55., 45., 0.))),
                    alpha1: Azimuth::new(dms(96., 36., 8.79960)),
                    p2: (Lon::new(dms(108., 13., 0.)), Lat::new(dms(-33., 26., 0.))),
                    alpha2: Azimuth::new(dms(137., 52., 22.01454)),
                    s: 14_110_526.170,
                },
                direct_deltas: DirectDeltas {
                    delta_lat2: dms(0., 0., -1.2e-5),
                    delta_delta_lon: dms(0., 0., 0.7e-5),
                    delta_alpha2: dms(0., 0., -1.2e-5),
                },
                inverse_deltas: InverseDeltas {
                    delta_alpha1: dms(0., 0., -0.4e-5),
                    delta_alpha2: dms(0., 0., -0.5e-5),
                    delta_s: -0.4e-3,
                },
            },
            // Line (b)
            LineData {
                ellipsoid: consts::INTL,
                geodesic: Geodesic {
                    p1: (Lon::new(0.), Lat::new(dms(37., 19., 54.95367))),
                    alpha1: Azimuth::new(dms(95., 27., 59.63089)),
                    p2: (Lon::new(dms(41., 28., 35.50729)), Lat::new(dms(26., 7., 42.83946))),
                    alpha2: Azimuth::new(dms(118., 5., 58.96161)),
                    s: 4_085_966.703,
                },
                direct_deltas: DirectDeltas {
                    delta_lat2: dms(0., 0., -0.7e-5),
                    delta_delta_lon: dms(0., 0., 1.2e-5),
                    delta_alpha2: dms(0., 0., 0.5e-5),
                },
                inverse_deltas: InverseDeltas {
                    delta_alpha1: dms(0., 0., -0.2e-5),
                    delta_alpha2: dms(0., 0., -0.2e-5),
                    delta_s: -0.4e-3,
                },
            },
            // Line (c)
            LineData {
                ellipsoid: consts::INTL,
                geodesic: Geodesic {
                    p1: (Lon::new(0.0), Lat::new(dms(35., 16., 11.24862))),
                    alpha1: Azimuth::new(dms(15., 44., 23.74850)),
                    p2: (Lon::new(dms(137., 47., 28.31435)), Lat::new(dms(67., 22., 14.77638))),
                    alpha2: Azimuth::new(dms(144., 55., 39.92147)),
                    s: 8_084_823.839,
                },
                direct_deltas: DirectDeltas {
                    delta_lat2: dms(0., 0., -2e-5),
                    delta_delta_lon: dms(0., 0., 2.9e-5),
                    delta_alpha2: dms(0., 0., 3e-5),
                },
                inverse_deltas: InverseDeltas {
                    delta_alpha1: dms(0., 0., -0.2e-5),
                    delta_alpha2: dms(0., 0., 0.3e-5),
                    delta_s: -0.7e-3,
                },
            },
        ]
    }

    #[test]
    fn solve_direct() {
        for line_data in lines() {
            let solver = VincentyGeodesicSolver::new(&line_data.ellipsoid);
            let result = solver.solve_direct(
                line_data.geodesic.p1,
                line_data.geodesic.alpha2,
                line_data.geodesic.s,
            );
            assert_abs_diff_eq!(
                result.p2.1,
                line_data.geodesic.p2.1,
                epsilon = line_data.direct_deltas.delta_lat2.abs()
            );
        }
    }
}
