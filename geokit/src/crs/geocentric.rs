use super::Crs;
use crate::cs::Coord;
use crate::cs::cartesian::geocentric::GeocentricAxes;
use crate::geodesy::GeodeticDatum;
use crate::quantities::length::Length;
use smol_str::SmolStr;

/// A [GeocentricCrs] is a **3D cartesian coordinates reference system** in which
/// coordinates are given by distance **in meters** along the following axes:
/// - X: axis from the center of the datum's ellipsoid in the equatorial and prime meridian plane,
/// - Y: axis from the center of the datum's ellipsoid in the equatorial plane and 90 degrees
/// meridian plane (east)
/// - Z: axis from the center of the datum's ellipsoid through the North Pole.
#[derive(Debug, Clone, PartialEq)]
pub struct GeocentricCrs {
    pub id: SmolStr,
    pub datum: GeodeticDatum,
    pub axes: GeocentricAxes,
}

impl Crs for GeocentricCrs {
    #[inline]
    fn id(&self) -> &str {
        &self.id
    }

    #[inline]
    fn dim(&self) -> usize {
        self.axes.dim()
    }

    /// Returns the cartesian distance between the 2 points.
    fn dist(&self, a: &Coord, b: &Coord) -> Result<Length, &'static str> {
        assert_eq!(a.len(), self.dim());
        assert_eq!(b.len(), self.dim());
        let na = self.axes.normalize(a);
        let nb = self.axes.normalize(b);
        Ok(na.dist_to(nb))
    }
}

#[cfg(test)]
mod tests {
    use super::GeocentricCrs;
    use crate::{
        crs::Crs,
        cs::cartesian::CartesianTolerance,
        geodesy::{GeodeticDatum, ellipsoid::consts::WGS84, prime_meridian::consts::GREENWICH},
        units::length::M,
    };
    use approx::assert_abs_diff_eq;
    use std::any::{Any, TypeId};

    #[test]
    fn type_id() {
        let crs: &dyn Any = &GeocentricCrs {
            id: "WGS84".into(),
            datum: GeodeticDatum::new("WGS84", WGS84, GREENWICH),
            axes: crate::cs::cartesian::geocentric::GeocentricAxes::XYZ,
        };
        assert_eq!(crs.type_id(), TypeId::of::<GeocentricCrs>());
        assert!(crs.is::<GeocentricCrs>());
    }

    #[test]
    fn approx_eq() {
        let crs = GeocentricCrs {
            id: "WGS84".into(),
            datum: GeodeticDatum::new("WGS84", WGS84, GREENWICH),
            axes: crate::cs::cartesian::geocentric::GeocentricAxes::XYZ,
        };

        let tol = crs.axes.denormalize_tol(CartesianTolerance::small());

        assert!(crs.approx_eq(&[1., 1., 1. + 1e-3], &[1., 1., 1.], &tol));
        assert!(!crs.approx_eq(&[1., 1., 1. + 1e-1], &[1., 1., 1.], &tol));
    }

    #[test]
    fn dist() {
        let crs = GeocentricCrs {
            id: "WGS84".into(),
            datum: GeodeticDatum::new("WGS84", WGS84, GREENWICH),
            axes: crate::cs::cartesian::geocentric::GeocentricAxes::XYZ,
        };
        let d = crs.dist(&[0., 0., 0.], &[1., 1., 1.]);
        assert!(d.is_ok());
        assert_abs_diff_eq!(d.unwrap(), 3.0f64.sqrt() * M);
    }
}
