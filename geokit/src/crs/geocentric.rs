use super::Crs;
use crate::cs::cartesian::geocentric::{GeocentricAxes, XYZ};
use crate::cs::cartesian::CartesianErrors;
use crate::geodesy::GeodeticDatum;
use crate::math::fp::Float;
use crate::quantities::length::Length;
use crate::transformations::{
    ToXYZTransformation, ToXYZTransformationProvider, TransformationError,
};
use approx::AbsDiffEq;
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
    type Tolerance = CartesianErrors;

    fn id(&self) -> &str {
        &self.id
    }

    fn dim(&self) -> usize {
        self.axes.dim()
    }

    fn approx_eq(&self, a: &[Float], b: &[Float], err: Self::Tolerance) -> bool {
        let err_m = err.length().m();
        a[0].abs_diff_eq(&b[0], err_m)
            && a[1].abs_diff_eq(&b[1], err_m)
            && a[2].abs_diff_eq(&b[2], err_m)
    }

    /// Returns the cartesian distance between the 2 points.
    fn dist(&self, a: &[Float], b: &[Float]) -> Length {
        let na = self.axes.normalize(a);
        let nb = self.axes.normalize(b);
        na.dist_to(nb)
    }
}

impl ToXYZTransformationProvider for GeocentricCrs {
    fn to_xyz_transformation<'a>(&self) -> Box<dyn ToXYZTransformation + 'a> {
        Box::new(self.clone())
    }
}

impl ToXYZTransformation for GeocentricCrs {
    fn to_xyz(&self, coords: &[Float]) -> Result<XYZ, TransformationError> {
        Ok(self.axes.normalize(coords))
    }

    fn from_xyz(&self, xyz: XYZ, coords: &mut [Float]) -> Result<(), TransformationError> {
        Ok(self.axes.denormalize(xyz, coords))
    }
}

#[cfg(test)]
mod tests {
    use super::GeocentricCrs;
    use crate::geodesy::{
        ellipsoid::consts::WGS84, prime_meridian::consts::GREENWICH, GeodeticDatum,
    };
    use std::any::{Any, TypeId};

    #[test]
    fn test_any() {
        let crs: &dyn Any = &GeocentricCrs {
            id: "WGS84".into(),
            datum: GeodeticDatum::new("WGS84", WGS84, GREENWICH),
            axes: crate::cs::cartesian::geocentric::GeocentricAxes::XYZ,
        };
        assert_eq!(crs.type_id(), TypeId::of::<GeocentricCrs>());
        assert!(crs.is::<GeocentricCrs>());
    }
}
