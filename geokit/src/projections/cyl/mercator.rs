use crate::cs::cartesian::ENH;
use crate::cs::geodetic::{Lat, Lon, LLH};
use crate::geodesy::Ellipsoid;
use crate::math::polynomial::Polynomial;
use crate::math::{Float, PI_2, PI_4};
use crate::projections::{Projection, ProjectionError};
use crate::quantities::length::{Arc, Length};
use crate::units::angle::RAD;

/// The [Mercator] map projection is a conformal projection of the ellipsoid on
/// a cylindrical surface whose axis is aligned with ellipsoid pole axis.
///
/// In the '1SP' case, the cylinder is tangent to the ellipsoid on the equator.
/// In the '2SP' case, the cylinder *cut* through the ellipsoid at the 2 standard
/// parallels.
///
/// Meridians are transformed into equally straight spaced lines.
/// Parallels are transformed into unequally spaced straight lines, closest to each other
/// near the equator and cut the transformed meridians at right angles.
///
/// # Sources
///
/// The formulae for this implementation are taken from 'EPSG guidance Note number 7, part 2 - November 2005'
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
    c_sin_nxi: [Float; 4],
}

impl Mercator {
    /// Creates a new [Mercator] projection instance known as *Mercator (1SP)* (EPSG:9804).
    /// The projection is defined with the equator as the single standard parallel and
    /// the scale at the equator `k0` is also defined.
    /// False grid coordinates are applied at the *natural origin* of the projection, the
    /// intersection of the equator and the *longitude of origin* `lon0`.
    ///
    /// # Arguments
    ///
    /// - 'ellipsoid': the ellipsoid used to derive some projection parameters, eg 'a', 'e'...
    /// - 'lon0': the longitude of *natural origin* (latitude of *natural origin* is 0).
    /// - 'k0': scale factor at the natural origin (eg on the equator)
    /// - 'false_easting':
    /// - 'false_northing'.
    pub(crate) fn new_1_sp(
        ellipsoid: &Ellipsoid,
        lon0: Lon,
        k0: Float,
        false_easting: Length,
        false_northing: Length,
    ) -> Self {
        Self::new(
            ellipsoid,
            lon0,
            Lat::ZERO,
            k0,
            false_easting,
            false_northing,
        )
    }

    /// Creates new [Mercator] projection instance known as *Mercator (2SP)* (EPSG:9805).
    /// The projecteion is defined through the latitude of 2 parallels equidistant on
    /// either side of the equator upon which the grid scale is true.
    /// False grid coordinates are applied at the *natural origin* of the projection,
    /// the intersection of the equator and the *longitude of origin* `lon0`.
    ///
    /// # Arguments
    ///
    /// - 'ellipsoid': the [Ellipsoid] used to derive some projection parameters, eg 'a', 'e'...
    /// - 'lon0': the longitude of *natural origin* (latitude of *natural origin* is 0).
    /// - 'lat0': the latitude of one standard parallel (the other being '-lat0').
    /// - 'k0': scale factor on the standard parallels
    /// - 'false_easting':
    /// - 'false_northing'.
    pub(crate) fn new_2_sp(
        ellipsoid: &Ellipsoid,
        lon0: Lon,
        lat0: Lat,
        false_easting: Length,
        false_northing: Length,
    ) -> Self {
        let (sin_lat1, cos_lat1) = lat0.abs().sin_cos();
        let e2 = ellipsoid.e_sq();
        let k0 = cos_lat1 / (1.0 - e2 * sin_lat1 * sin_lat1).sqrt();
        Self::new(ellipsoid, lon0, lat0, k0, false_easting, false_northing)
    }

    fn new(
        ellipsoid: &Ellipsoid,
        lon0: Lon,
        lat0: Lat,
        k0: Float,
        false_easting: Length,
        false_northing: Length,
    ) -> Self {
        let e2 = ellipsoid.e_sq();
        let c_sin_nxi = [
            Polynomial::new([0.0, 0.5, 5.0 / 24.0, 1.0 / 12.0, 13.0 / 360.0]).eval_at(e2),
            Polynomial::new([0.0, 0.0, 7.0 / 48.0, 29.0 / 240.0, 811.0 / 11520.0]).eval_at(e2),
            Polynomial::new([0.0, 0.0, 0.0, 7.0 / 120.0, 81.0 / 1120.0]).eval_at(e2),
            Polynomial::new([0.0, 0.0, 0.0, 0.0, 4279.0 / 161280.0]).eval_at(e2),
        ];
        Self {
            a: ellipsoid.a(),
            e: ellipsoid.e(),
            e2,
            lon0,
            lat0: lat0.abs(),
            k0,
            false_easting,
            false_northing,
            c_sin_nxi,
        }
    }
}

impl Projection for Mercator {
    fn proj(&self, input: LLH) -> Result<ENH, ProjectionError> {
        let sin_lat = input.lat.sin();
        let e_sin_lat = self.e * sin_lat;
        let lat1 = input.lat / 2.0 + PI_4 * RAD;
        let r = ((1.0 - e_sin_lat) / (1.0 + e_sin_lat)).powf(self.e / 2.0);

        Ok(ENH {
            easting: self.false_easting
                + ((self.a * self.k0) * (input.lon - self.lon0).angle()).length(),
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

        let x = xi
            + self.c_sin_nxi[0] * sin_2xi
            + self.c_sin_nxi[1] * sin_4xi
            + self.c_sin_nxi[2] * sin_6xi
            + self.c_sin_nxi[3] * sin_8xi;

        Ok(LLH {
            lon: self.lon0 + Arc(input.easting - self.false_easting) / (self.a * self.k0),
            lat: Lat::new(x * RAD),
            height: input.height,
        })
    }
}

#[cfg(test)]
mod test {
    use approx::assert_abs_diff_eq;

    use super::Mercator;
    use crate::{
        cs::{
            cartesian::ENH,
            geodetic::{Lat, Lon, LLH},
        },
        geodesy::ellipsoid,
        projections::Projection,
        units::{angle::DEG, length::M},
    };

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
        assert_abs_diff_eq!(enh.easting.m(), 165_704.29, epsilon = 4e-3);
        assert_abs_diff_eq!(enh.northing.m(), 5_171_848.07, epsilon = 3e-3);

        let llh = proj
            .unproj(ENH {
                easting: 165_704.29 * M,
                northing: 5_171_848.07 * M,
                height: 0.0 * M,
            })
            .unwrap();
        assert_abs_diff_eq!(llh.lon.val(DEG), 53.0, epsilon = 5e-8);
        assert_abs_diff_eq!(llh.lat.val(DEG), 53.0, epsilon = 4e-8);
    }
}
