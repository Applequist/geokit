use super::Crs;
use crate::cs::cartesian::{GeocentricAxes, XYZ};
use crate::geodesy::GeodeticDatum;
use crate::math::fp::Float;
use crate::transformations::{
    ToXYZTransformation, ToXYZTransformationProvider, TransformationError,
};
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
    fn id(&self) -> &str {
        &self.id
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
