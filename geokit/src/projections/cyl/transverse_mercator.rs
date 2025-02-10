#![allow(non_snake_case)]

use crate::{
    cs::{
        cartesian::ENH,
        geodetic::{Lat, Lon, LLH},
    },
    geodesy::Ellipsoid,
    math::{polynomial::Polynomial, Float},
    projections::{Projection, ProjectionError},
    quantities::{length::Length},
    units::{
        angle::{DEG, RAD},
        length::M,
    },
};

/// The [TransverseMercator] projection (EPSG:9807).
///
/// This implementation is based on the JHS formulae as recommended in
/// EPSG Guidance 7 part 2 of November 2019.
///
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct TransverseMercator {
    e: Float,
    lon0: Lon,
    k0: Float,
    false_easting: Length,
    false_northing: Length,
    B: Length,
    h: [Float; 4],
    h_p: [Float; 4],
    M0: Length,
}

impl TransverseMercator {
    pub(crate) fn new_utm_north(ellipsoid: &Ellipsoid, zone: UTMZone) -> Self {
        Self::new(
            ellipsoid,
            utm_lon0(zone),
            UTM_LAT0,
            UTM_K0,
            UTM_FE,
            UTM_FN_NORTH,
        )
    }

    pub(crate) fn new_utm_south(ellipsoid: &Ellipsoid, zone: u8) -> Self {
        Self::new(
            ellipsoid,
            utm_lon0(zone),
            UTM_LAT0,
            UTM_K0,
            UTM_FE,
            UTM_FN_SOUTH,
        )
    }

    pub(crate) fn new(
        ellipsoid: &Ellipsoid,
        lon0: Lon,
        lat0: Lat,
        k0: Float,
        false_easting: Length,
        false_northing: Length,
    ) -> Self {
        let e = ellipsoid.e();
        let n = ellipsoid.n();
        let B = ellipsoid.a() / (1.0 + n)
            * Polynomial::new([1.0, 0.0, 0.25, 0.0, 1.0 / 64.0]).fast_eval_at(n);
        let h = [
            Polynomial::new([0.0, 0.5, -2.0 / 3.0, 5.0 / 16.0, 41.0 / 180.0]).fast_eval_at(n),
            Polynomial::new([0.0, 0.0, 13.0 / 48.0, -3.0 / 5.0, 557.0 / 1440.0]).fast_eval_at(n),
            Polynomial::new([0.0, 0.0, 0.0, 61.0 / 240.0, -103.0 / 140.0]).fast_eval_at(n),
            Polynomial::new([0.0, 0.0, 0.0, 0.0, 49561.0 / 161280.0]).fast_eval_at(n),
        ];

        let h_p = [
            Polynomial::new([0.0, 0.5, -2.0 / 3.0, 37.0 / 96.0, -1.0 / 360.0]).fast_eval_at(n),
            Polynomial::new([0.0, 0.0, 1.0 / 48.0, 1.0 / 15.0, -437.0 / 1440.0]).fast_eval_at(n),
            Polynomial::new([0.0, 0.0, 0.0, 17.0 / 480.0, -37.0 / 840.0]).fast_eval_at(n),
            Polynomial::new([0.0, 0.0, 0.0, 0.0, 4397.0 / 161280.0]).fast_eval_at(n),
        ];

        let M0 = M0(lat0, B, e, h);

        Self {
            e,
            lon0,
            k0,
            false_easting,
            false_northing,
            B,
            h,
            h_p,
            M0,
        }
    }
}

impl Projection for TransverseMercator {
    fn proj(&self, input: LLH) -> Result<ENH, ProjectionError> {
        let Q = Q(self.e, input.lat);
        let beta = beta(Q);
        let eta_0 = {
            let lon = input.lon;
            let lon0 = self.lon0;
            (beta.cos() * (lon - lon0).sin()).atanh()
        };
        let xi_0 = (beta.sin() * eta_0.cosh()).asin();

        let xis = [
            self.h[0] * (2.0 * xi_0).sin() * (2.0 * eta_0).cosh(),
            self.h[1] * (4.0 * xi_0).sin() * (4.0 * eta_0).cosh(),
            self.h[2] * (6.0 * xi_0).sin() * (6.0 * eta_0).cosh(),
            self.h[3] * (8.0 * xi_0).sin() * (8.0 * eta_0).cosh(),
        ];
        let xi = xi_0 + xis.into_iter().sum::<Float>();

        let etas = [
            self.h[0] * (2.0 * xi_0).cos() * (2.0 * eta_0).sinh(),
            self.h[1] * (4.0 * xi_0).cos() * (4.0 * eta_0).sinh(),
            self.h[2] * (6.0 * xi_0).cos() * (6.0 * eta_0).sinh(),
            self.h[3] * (8.0 * xi_0).cos() * (8.0 * eta_0).sinh(),
        ];
        let eta = eta_0 + etas.into_iter().sum::<Float>();

        Ok(ENH {
            easting: self.false_easting + self.k0 * self.B * eta,
            northing: self.false_northing + self.k0 * (self.B * xi - self.M0),
            height: input.height,
        })
    }

    fn unproj(&self, input: ENH) -> Result<LLH, ProjectionError> {
        let eta_p = (input.easting - self.false_easting) / (self.B * self.k0);
        let xi_p = (input.northing - self.false_northing + self.k0 * self.M0) / (self.B * self.k0);

        let xi_ps = [
            self.h_p[0] * (2.0 * xi_p).sin() * (2.0 * eta_p).cosh(),
            self.h_p[1] * (4.0 * xi_p).sin() * (4.0 * eta_p).cosh(),
            self.h_p[2] * (6.0 * xi_p).sin() * (6.0 * eta_p).cosh(),
            self.h_p[3] * (8.0 * xi_p).sin() * (8.0 * eta_p).cosh(),
        ];
        let xi_p0 = xi_p - xi_ps.into_iter().sum::<Float>();

        let eta_ps = [
            self.h_p[0] * (2.0 * xi_p).cos() * (2.0 * eta_p).sinh(),
            self.h_p[1] * (4.0 * xi_p).cos() * (4.0 * eta_p).sinh(),
            self.h_p[2] * (6.0 * xi_p).cos() * (6.0 * eta_p).sinh(),
            self.h_p[3] * (8.0 * xi_p).cos() * (8.0 * eta_p).sinh(),
        ];
        let eta_p0 = eta_p - eta_ps.into_iter().sum::<Float>();

        let beta_p = (xi_p0.sin() / eta_p0.cosh()).asin();
        let Q_p = beta_p.tan().asinh();

        fn iter_fn(e: Float, Q_p: Float, Q_pp: Float) -> Float {
            Q_p + (e * (e * Q_pp.tanh()).atanh())
        }

        let mut Q_pp = iter_fn(self.e, Q_p, Q_p);
        loop {
            let old_Q_pp = Q_pp;
            Q_pp = iter_fn(self.e, Q_p, Q_pp);
            if (Q_pp - old_Q_pp).abs() < 1e-12 {
                break;
            }
        }

        Ok(LLH {
            lon: self.lon0 + (eta_p0.tanh() / beta_p.cos()).asin() * RAD,
            lat: Lat::new(Q_pp.sinh().atan() * RAD),
            height: input.height,
        })
    }
}

//
// UTM zone calculations
//

type UTMZone = u8;

const UTM_ZONE_WIDTH_DEG: Float = 6.0;
const UTM_ZONE_COUNT: Float = 360.0 / UTM_ZONE_WIDTH_DEG;
const UTM_ZONE_1_WEST_LIMIT_DEG: Float = -177.0;
const UTM_LAT0: Lat = Lat::ZERO;
const UTM_K0: Float = 0.9996;
const UTM_FE: Length = Length::new(500_000.0, M);
const UTM_FN_NORTH: Length = Length::ZERO;
const UTM_FN_SOUTH: Length = Length::new(10_000_000.0, M);

/// Compute the zone number from a longitude.
fn utm_zone(lon: Lon) -> UTMZone {
    let lon_deg = lon.angle().val(DEG);
    let delta_lon = ((lon_deg - UTM_ZONE_1_WEST_LIMIT_DEG) / UTM_ZONE_WIDTH_DEG).floor();
    let z = delta_lon - (delta_lon / UTM_ZONE_COUNT).floor() * UTM_ZONE_COUNT;
    z as u8
}

/// Compute the central longitude of a zone by
fn utm_lon0(zone: UTMZone) -> Lon {
    let lon0_deg =
        UTM_ZONE_1_WEST_LIMIT_DEG + (zone as Float) * UTM_ZONE_WIDTH_DEG - UTM_ZONE_WIDTH_DEG / 2.0;
    Lon::new(lon0_deg * DEG)
}

//
// Helper functions for the projection
//

fn Q(e: Float, lat: Lat) -> Float {
    lat.tan().asinh() - e * (e * lat.sin()).atanh()
}

fn beta(Q: Float) -> Float {
    Q.sinh().atan()
}

fn M0(lat0: Lat, B: Length, e: Float, h: [Float; 4]) -> Length {
    if lat0 == Lat::ZERO {
        0.0 * M
    } else if lat0.abs() == Lat::MAX {
        B * lat0.angle().rad()
    } else {
        let Q0 = Q(e, lat0);
        let beta0 = beta(Q0);
        let xi_o0 = beta0; // From note: simplified beta0.sin().asin() to beta0;
        let xi_os = [
            xi_o0,
            h[0] * (2.0 * xi_o0).sin(),
            h[1] * (4.0 * xi_o0).sin(),
            h[2] * (6.0 * xi_o0).sin(),
            h[3] * (8.0 * xi_o0).sin(),
        ];
        let xi_o = xi_os.into_iter().sum::<Float>();
        B * xi_o
    }
}

#[cfg(test)]
mod test {
    use super::TransverseMercator;
    use crate::{
        cs::{
            cartesian::ENH,
            geodetic::{GeodeticErrors, Lat, Lon, LLH},
        },
        geodesy::ellipsoid,
        projections::Projection,
        units::{angle::DEG, length::M},
    };
    use approx::assert_abs_diff_eq;

    #[test]
    fn transverse_mercator() {
        // From EPSG Guidance Note 7 part 2
        let proj = TransverseMercator::new(
            &ellipsoid::consts::AIRY,
            Lon::new(-2.0 * DEG), // 2 W
            Lat::new(49.0 * DEG), // 49 N
            0.9996012717,
            400_000.0 * M,
            -100_000.0 * M,
        );

        let enh = proj
            .proj(LLH {
                lon: Lon::dms(0., 30., 0.),  //  0d 30' E
                lat: Lat::dms(50., 30., 0.), // 50d 30' N
                height: 0.0 * M,
            })
            .unwrap();
        assert_abs_diff_eq!(enh.easting.m(), 577_274.99, epsilon = 1e-2);
        assert_abs_diff_eq!(enh.northing.m(), 69_740.5, epsilon = 1e-2);

        let llh = proj
            .unproj(ENH {
                easting: 577_274.99 * M,
                northing: 69_740.50 * M,
                height: 0.0 * M,
            })
            .unwrap();
        assert_abs_diff_eq!(llh.lon, Lon::new(0.5 * DEG), epsilon = 1e-7 * DEG);
        assert_abs_diff_eq!(llh.lat, Lat::new(50.5 * DEG), epsilon = 7e-8 * DEG);
    }

    #[test]
    fn transverse_mercator_bound() {
        let proj = TransverseMercator::new(
            &ellipsoid::consts::WGS84,
            Lon::ZERO,
            Lat::ZERO,
            0.9996,
            0.0 * M,
            0.0 * M,
        );

        let north_pole_enh = proj
            .proj(LLH {
                lon: Lon::ZERO,
                lat: Lat::MAX,
                height: 0.0 * M,
            })
            .unwrap();
        println!("north_pole = {north_pole_enh:?}");
        let east_enh = proj
            .proj(LLH {
                lon: Lon::new(90. * DEG),
                lat: Lat::new(10.0 * DEG),
                height: 0.0 * M,
            })
            .unwrap();
        println!("90E 10N = {east_enh:?}");
        let west_enh = proj
            .proj(LLH {
                lon: Lon::new(-90. * DEG),
                lat: Lat::new(10. * DEG),
                height: 0.0 * M,
            })
            .unwrap();
        println!("90W 10N = {west_enh:?}");

        let south_pole_enh = proj
            .proj(LLH {
                lon: Lon::ZERO,
                lat: Lat::new(-90.0 * DEG),
                height: 0.0 * M,
            })
            .unwrap();
        println!("south_pole = {south_pole_enh:?}");
    }

    #[test]
    fn transverse_mercator_roundtrip() {
        let proj = TransverseMercator::new(
            &ellipsoid::consts::WGS84,
            Lon::ZERO,
            Lat::ZERO,
            0.9996,
            0.0 * M,
            0.0 * M,
        );

        let north_pole = LLH {
            lon: Lon::new(10. * DEG),
            lat: Lat::new(90.0 * DEG),
            height: 0.0 * M,
        };
        let north_pole_enh = proj.proj(north_pole).unwrap();
        println!("north_pole (enh) = {north_pole_enh:?}");

        let north_pole_llh = proj.unproj(north_pole_enh).unwrap();
        println!("north_pole (llh)= {north_pole_llh:?}");

        north_pole_llh.approx_eq(&north_pole, GeodeticErrors::default());

        let south_pole = LLH {
            lon: Lon::new(-10. * DEG),
            lat: Lat::new(-90.0 * DEG),
            height: 0.0 * M,
        };
        let south_pole_enh = proj.proj(south_pole).unwrap();
        println!("south_pole (enh) = {south_pole_enh:?}");

        let south_pole_llh = proj.unproj(south_pole_enh).unwrap();
        println!("south_pole (llh)= {south_pole_llh:?}");

        south_pole_llh.approx_eq(&south_pole, GeodeticErrors::default());
    }
}
