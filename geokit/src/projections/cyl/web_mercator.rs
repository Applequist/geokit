use crate::{
    cs::{
        cartesian::ENH,
        geodetic::{Lat, Lon, LLH},
    },
    geodesy::Ellipsoid,
    math::PI_4,
    projections::{Projection, ProjectionError},
    quantities::{
        angle::Angle,
        length::{Arc, Length},
    },
    units::angle::RAD,
};

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
}

impl Projection for WebMercator {
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
    use super::Projection;
    use super::WebMercator;
    use crate::cs::cartesian::ENH;
    use crate::cs::geodetic::Lat;
    use crate::cs::geodetic::Lon;
    use crate::cs::geodetic::LLH;
    use crate::geodesy::ellipsoid;
    use crate::quantities::length::Length;
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
}
