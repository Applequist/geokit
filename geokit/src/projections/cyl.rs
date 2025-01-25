use crate::cs::cartesian::ENH;
use crate::cs::geodetic::{Lat, Lon, LLH};
use crate::geodesy::Ellipsoid;
use crate::math::polynomial::Polynomial;
use crate::math::{Float, PI_2, PI_4};
use crate::quantities::angle::Angle;
use crate::quantities::length::{Arc, Length};
use crate::units::angle::{DEG, RAD};
use crate::units::length::M;

use super::ProjectionError;

/// The [Mercator] map projection.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct Mercator {
    a: Length,
    e: Float,
    e2: Float,
    lon0: Lon,
    lat0: Lat,
    k0: Float,
    false_easting: Length,
    false_northing: Length,
}

impl Mercator {
    /// Creates a new [Mercator] projection instance known as *Mercator (1SP)*. The projection is
    /// defined with the equator as the single standard parallel and the scale at the equator also
    /// defined. False grid coordinates are applied at the *natural origin* of the projection, the
    /// intersection of the equator and the *longitude of origin* `lon0`.
    /// Know as **EPSG:9804**.
    pub(crate) fn new_1_sp(
        ellipsoid: &Ellipsoid,
        lon0: Lon,
        k0: Float,
        false_easting: Length,
        false_northing: Length,
    ) -> Self {
        Self {
            a: ellipsoid.a(),
            e: ellipsoid.e(),
            e2: ellipsoid.e_sq(),
            lon0,
            lat0: Lat::ZERO,
            k0,
            false_easting,
            false_northing,
        }
    }

    /// Creates new [Mercator] projection instance defined through the latitude of 2 parallel
    /// equidistant either side of the equator upon which the grid scale is true. False grid
    /// coordinates are applied at the *natural origin* of the projection, the intersection of the
    /// equator and the *longitude of origin* `lon0`.
    /// Known as EPSG:9805
    pub(crate) fn new_2_sp(
        ellipsoid: &Ellipsoid,
        lon0: Lon,
        lat0: Lat,
        false_easting: Length,
        false_northing: Length,
    ) -> Self {
        let (sin_lat1, cos_lat1) = lat0.angle().abs().sin_cos();
        Self {
            a: ellipsoid.a(),
            e: ellipsoid.e(),
            e2: ellipsoid.e_sq(),
            lon0,
            lat0: lat0.abs(),
            k0: cos_lat1 / (1.0 - ellipsoid.e_sq() * sin_lat1 * sin_lat1).sqrt(),
            false_easting,
            false_northing,
        }
    }

    fn proj(&self, input: LLH) -> Result<ENH, ProjectionError> {
        let sin_lat = input.lat.sin();
        let e_sin_lat = self.e * sin_lat;
        let lat1 = input.lat / 2.0 + PI_4 * RAD;
        let r = ((1.0 - e_sin_lat) / (1.0 + e_sin_lat)).powf(self.e / 2.0);

        Ok(ENH {
            easting: self.false_easting + self.a * self.k0 * (input.lon - self.lon0.angle()).rad(),
            northing: self.false_northing + self.a * self.k0 * (lat1.tan() * r).ln(),
            height: input.height,
        })
    }

    fn unproj(&self, input: ENH) -> Result<LLH, ProjectionError> {
        let t = ((self.false_northing - input.northing) / (self.a * self.k0)).exp();
        let xi = PI_2 - 2.0 * t.atan();
        let sin_2xi = (2.0 * xi).sin();
        let sin_4xi = (4.0 * xi).sin();
        let sin_6xi = (6.0 * xi).sin();
        let sin_8xi = (8.0 * xi).sin();

        let e2 = self.e * self.e;

        let x = xi
            + Polynomial::new([0.0, 0.5, 5.0 / 24.0, 1.0 / 12.0, 13.0 / 360.0]).eval_at(e2)
                * sin_2xi
            + Polynomial::new([0.0, 0.0, 7.0 / 48.0, 29.0 / 240.0, 811.0 / 11520.0]).eval_at(e2)
                * sin_4xi
            + Polynomial::new([0.0, 0.0, 0.0, 7.0 / 120.0, 81.0 / 1120.0]).eval_at(e2) * sin_6xi
            + Polynomial::new([0.0, 0.0, 0.0, 0.0, 4279.0 / 161280.0]).eval_at(e2) * sin_8xi;

        Ok(LLH {
            lon: self.lon0 + Arc(input.easting - self.false_easting) / (self.a * self.k0),
            lat: Lat::new(x * RAD),
            height: input.height,
        })
    }
}

/// The [TransverseMercator] projection.
/// Known as EPSG:9807.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct TransverseMercator {
    e: Float,
    lon0: Lon,
    k0: Float,
    false_easting: Length,
    false_northing: Length,
    upper_b: Length,
    h: [Float; 4],
    h_p: [Float; 4],
    m0: Length,
}

impl TransverseMercator {
    const UTM_ZONE_WIDTH_DEG: Float = 6.0;
    const UTM_ZONE_COUNT: Float = 360.0 / Self::UTM_ZONE_WIDTH_DEG;
    const UTM_ZONE_1_WEST_LIMIT_DEG: Float = -177.0;
    const UTM_LAT0: Lat = Lat::ZERO;
    const UTM_K0: Float = 0.9996;
    const UTM_FE: Length = Length::new(500_000.0, M);
    const UTM_FN_NORTH: Length = Length::ZERO;
    const UTM_FN_SOUTH: Length = Length::new(10_000_000.0, M);

    fn utm_zone(lon: Lon) -> u8 {
        let lon_deg = lon.val(DEG);
        let delta_lon =
            ((lon_deg - Self::UTM_ZONE_1_WEST_LIMIT_DEG) / Self::UTM_ZONE_WIDTH_DEG).floor();
        let z = delta_lon - (delta_lon / Self::UTM_ZONE_COUNT).floor() * Self::UTM_ZONE_COUNT;
        z as u8
    }

    fn utm_lon0(zone: u8) -> Lon {
        let lon0_deg = Self::UTM_ZONE_1_WEST_LIMIT_DEG + (zone as Float) * Self::UTM_ZONE_WIDTH_DEG
            - Self::UTM_ZONE_WIDTH_DEG / 2.0;
        Lon::new(lon0_deg * DEG)
    }

    pub(crate) fn new_utm_north(ellipsoid: &Ellipsoid, zone: u8) -> Self {
        Self::new(
            ellipsoid,
            Self::utm_lon0(zone),
            Self::UTM_LAT0,
            Self::UTM_K0,
            Self::UTM_FE,
            Self::UTM_FN_NORTH,
        )
    }

    pub(crate) fn new_utm_south(ellipsoid: &Ellipsoid, zone: u8) -> Self {
        Self::new(
            ellipsoid,
            Self::utm_lon0(zone),
            Self::UTM_LAT0,
            Self::UTM_K0,
            Self::UTM_FE,
            Self::UTM_FN_SOUTH,
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
        let f = 1.0 / ellipsoid.invf();
        let n = f / (2.0 - f);
        let upper_b = ellipsoid.a() / (1.0 + n)
            * Polynomial::new([1.0, 0.0, 0.25, 0.0, 1.0 / 64.0]).eval_at(n);
        let h = [
            Polynomial::new([0.0, 0.5, -2.0 / 3.0, 5.0 / 16.0, 41.0 / 180.0]).eval_at(n),
            Polynomial::new([0.0, 0.0, 13.0 / 48.0, -3.0 / 5.0, 557.0 / 1440.0]).eval_at(n),
            Polynomial::new([0.0, 0.0, 0.0, 61.0 / 240.0, -103.0 / 140.0]).eval_at(n),
            Polynomial::new([0.0, 0.0, 0.0, 0.0, 49561.0 / 161280.0]).eval_at(n),
        ];

        let h_p = [
            Polynomial::new([0.0, 0.5, -2.0 / 3.0, 37.0 / 96.0, -1.0 / 360.0]).eval_at(n),
            Polynomial::new([0.0, 0.0, 1.0 / 48.0, 1.0 / 15.0, -437.0 / 1440.0]).eval_at(n),
            Polynomial::new([0.0, 0.0, 0.0, 17.0 / 480.0, -37.0 / 840.0]).eval_at(n),
            Polynomial::new([0.0, 0.0, 0.0, 0.0, 4397.0 / 161280.0]).eval_at(n),
        ];

        let m0 = Self::m0(lat0, upper_b, e, h);

        Self {
            e,
            lon0,
            k0,
            false_easting,
            false_northing,
            upper_b,
            h,
            h_p,
            m0,
        }
    }

    fn q(e: Float, lat: Lat) -> Float {
        lat.tan().asinh() - e * (e * lat.sin()).atanh()
    }

    fn beta(q: Float) -> Float {
        q.sinh().atan()
    }

    fn eta_0(beta: Float, lon: Lon, lon0: Lon) -> Float {
        (beta.cos() * (lon - lon0.angle()).sin()).atanh()
    }

    fn m0(lat0: Lat, upper_b: Length, e: Float, h: [Float; 4]) -> Length {
        if lat0 == Lat::ZERO {
            0.0 * M
        } else if lat0.abs() == Lat::MAX {
            upper_b * lat0.rad()
        } else {
            let q0 = Self::q(e, lat0);
            let beta0 = Self::beta(q0);
            let xi_o0 = beta0; // From note: simplified beta0.sin().asin() to beta0;
            let xi_os = [
                xi_o0,
                h[0] * (2.0 * xi_o0).sin(),
                h[1] * (4.0 * xi_o0).sin(),
                h[2] * (6.0 * xi_o0).sin(),
                h[3] * (8.0 * xi_o0).sin(),
            ];
            let xi_o = xi_os.into_iter().sum::<Float>();
            upper_b * xi_o
        }
    }

    fn proj(&self, input: LLH) -> Result<ENH, ProjectionError> {
        let q = Self::q(self.e, input.lat);
        let beta = Self::beta(q);
        let eta_0 = Self::eta_0(beta, input.lon, self.lon0);
        let xi_0 = (beta.sin() * eta_0.cosh()).asin();
        let xis = [
            xi_0,
            self.h[0] * (2.0 * xi_0).sin() * (2.0 * eta_0).cosh(),
            self.h[1] * (4.0 * xi_0).sin() * (4.0 * eta_0).cosh(),
            self.h[2] * (6.0 * xi_0).sin() * (6.0 * eta_0).cosh(),
            self.h[3] * (8.0 * xi_0).sin() * (8.0 * eta_0).cosh(),
        ];
        let xi = xis.into_iter().sum::<Float>();
        let etas = [
            self.h[0] * (2.0 * xi_0).cos() * (2.0 * eta_0).sinh(),
            self.h[1] * (4.0 * xi_0).cos() * (4.0 * eta_0).sinh(),
            self.h[2] * (6.0 * xi_0).cos() * (6.0 * eta_0).sinh(),
            self.h[3] * (8.0 * xi_0).cos() * (8.0 * eta_0).sinh(),
        ];
        let eta = eta_0 + etas.into_iter().sum::<Float>();

        Ok(ENH {
            easting: self.false_easting + self.k0 * self.upper_b * eta,
            northing: self.false_northing + self.k0 * (self.upper_b * xi - self.m0),
            height: input.height,
        })
    }

    fn unproj(&self, input: ENH) -> Result<LLH, ProjectionError> {
        let eta_p = (input.easting - self.false_easting) / (self.upper_b * self.k0);
        let xi_p =
            (input.northing - self.false_northing + self.k0 * self.m0) / (self.upper_b * self.k0);

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
        let q_p = beta_p.tan().asinh();

        fn iter_fn(e: Float, q_p: Float, q_pp: Float) -> Float {
            q_p + (e * (e * q_pp.tanh()).atanh())
        }

        let mut q_pp = iter_fn(self.e, q_p, q_p);
        loop {
            let old_q_pp = q_pp;
            q_pp = iter_fn(self.e, q_p, q_pp);
            if (q_pp - old_q_pp).abs() < 1e-10 {
                break;
            }
        }

        Ok(LLH {
            lon: self.lon0 + (eta_p0.tanh() / beta_p.cos()).asin() * RAD,
            lat: Lat::new(q_pp.sinh().atan() * RAD),
            height: input.height,
        })
    }
}

/// [WebMercator](epsg:1024) also known as 'Pseudo-Mercator' is a projection method
/// used by some popular web mapping and visualisation applications.
/// Strictly speaking the name is misleading as it is **NOT** a Mercator projection.
#[derive(Debug, Clone, PartialEq)]
pub struct WebMercator {
    a: Length,
    /// longitude of *natural origin* in radians.
    lon0: Lon,
    /// latitude of *natural origin* in radians.
    lat0: Lat,
    /// False easting in meters.
    false_easting: Length,
    /// False norhting in meters.
    false_northing: Length,
}

impl WebMercator {
    /// Creates a new [WebMercator] projection instance.
    pub(crate) fn new(
        ellipsoid: &Ellipsoid,
        lon0: Lon,
        lat0: Lat,
        false_easting: Length,
        false_northing: Length,
    ) -> Self {
        Self {
            a: ellipsoid.a(),
            lon0,
            lat0,
            false_easting,
            false_northing,
        }
    }

    fn proj(&self, input: LLH) -> Result<ENH, ProjectionError> {
        Ok(ENH {
            easting: self.false_easting + (self.a * (input.lon - self.lon0).angle()).length(),
            northing: self.false_northing + self.a * (input.lat / 2.0 + PI_4 * RAD).tan().ln(),
            height: input.height,
        })
    }

    fn unproj(&self, input: ENH) -> Result<LLH, ProjectionError> {
        let d = (self.false_northing - input.northing) / self.a;

        Ok(LLH {
            lon: self.lon0 + Arc(input.easting - self.false_easting) / self.a,
            lat: Lat::new(Angle::PI_2 - 2.0 * d.exp().atan() * RAD),
            height: input.height,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::Mercator;
    use super::TransverseMercator;
    use super::WebMercator;
    use crate::cs::cartesian::ENH;
    use crate::cs::geodetic::Lat;
    use crate::cs::geodetic::Lon;
    use crate::cs::geodetic::LLH;
    use crate::geodesy::ellipsoid;
    use crate::quantities::angle::Angle;
    use crate::quantities::length::Length;
    use crate::units::angle::DEG;
    use crate::units::angle::RAD;
    use crate::units::length::M;
    use approx::assert_abs_diff_eq;

    #[test]
    fn web_mercator() {
        let proj = WebMercator::new(
            &ellipsoid::consts::WGS84,
            Lon::ZERO,
            Lat::ZERO,
            0.0 * M,
            0.0 * M,
        );
        let enh = proj
            .proj(LLH {
                lon: Lon::new(-1.751147016 * RAD),
                lat: Lat::new(0.42554246 * RAD),
                height: Length::ZERO,
            })
            .unwrap();
        assert_abs_diff_eq!(enh.easting.m(), -11_169_055.58, epsilon = 1e-2);
        assert_abs_diff_eq!(enh.northing.m(), 2_800_000.0, epsilon = 1e-2);

        let llh = proj
            .unproj(ENH {
                easting: -11_169_055.58 * M,
                northing: 2_810_000.0 * M,
                height: Length::ZERO,
            })
            .unwrap();
        assert_abs_diff_eq!(llh.lon.rad(), -1.751147016, epsilon = 1e-9);
        assert_abs_diff_eq!(llh.lat.rad(), 0.426970023, epsilon = 1e-9);
    }

    #[test]
    fn mercator_1_sp() {
        // From EPSG Guidance Note 7 Part 2 3.5.1 Mercator
        let proj = Mercator::new_1_sp(
            &ellipsoid::consts::BESSEL,
            Lon::new(110.0 * DEG),
            0.997,
            3_900_000.0 * M,
            900_000.0 * M,
        );
        let enh = proj
            .proj(LLH {
                lon: Lon::new(120.0 * DEG),
                lat: Lat::new(-3.0 * DEG),
                height: 0.0 * M,
            })
            .unwrap();
        assert_abs_diff_eq!(enh.easting.m(), 5009726.58, epsilon = 2e-2);
        assert_abs_diff_eq!(enh.northing.m(), 569150.82, epsilon = 2e-2);

        let llh = proj
            .unproj(ENH {
                easting: 5009726.58 * M,
                northing: 569150.82 * M,
                height: 0.0 * M,
            })
            .unwrap();
        assert_abs_diff_eq!(llh.lon.val(DEG), 120.0, epsilon = 3e-8);
        assert_abs_diff_eq!(llh.lat.val(DEG), -3.0, epsilon = 3e-8);
    }

    #[test]
    fn mercator_2_sp() {
        // From EPSG Guidance Note 7 Part 2 3.5.1 Mercator
        let proj = Mercator::new_2_sp(
            &ellipsoid::consts::KRASS,
            Lon::new(51.0 * DEG),
            Lat::new(42.0 * DEG),
            0.0 * M,
            0.0 * M,
        );
        let enh = proj
            .proj(LLH {
                lon: Lon::new(53.0 * DEG),
                lat: Lat::new(53.0 * DEG),
                height: 0.0 * M,
            })
            .unwrap();
        assert_abs_diff_eq!(enh.easting.m(), 165704.29, epsilon = 4e-3);
        assert_abs_diff_eq!(enh.northing.m(), 5171848.07, epsilon = 3e-3);

        let llh = proj
            .unproj(ENH {
                easting: 165704.29 * M,
                northing: 5171848.07 * M,
                height: 0.0 * M,
            })
            .unwrap();
        assert_abs_diff_eq!(llh.lon.val(DEG), 53.0, epsilon = 5e-8);
        assert_abs_diff_eq!(llh.lat.val(DEG), 53.0, epsilon = 4e-8);
    }

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
                lon: Lon::new(0.5 * DEG),
                lat: Lat::new(50.5 * DEG),
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
        assert_abs_diff_eq!(llh.lon.val(DEG), 0.5, epsilon = 1e-7);
        assert_abs_diff_eq!(llh.lat.val(DEG), 50.5, epsilon = 7e-8);
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

        assert_abs_diff_eq!(&north_pole_llh, &north_pole, epsilon = 1e-9);

        let south_pole = LLH {
            lon: Lon::new(-10. * DEG),
            lat: Lat::new(-90.0 * DEG),
            height: 0.0 * M,
        };
        let south_pole_enh = proj.proj(south_pole).unwrap();
        println!("south_pole (enh) = {south_pole_enh:?}");

        let south_pole_llh = proj.unproj(south_pole_enh).unwrap();
        println!("south_pole (llh)= {south_pole_llh:?}");

        assert_abs_diff_eq!(&south_pole_llh, &south_pole, epsilon = 1e-9);
    }
}
