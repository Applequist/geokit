//! This module defines types to perform coordinates transformation through normalized geocentric
//! coordinates.

use super::{CrsTransformation, TransformationError};
use crate::{
    crs::{geocentric::GeocentricCrs, geographic::GeographicCrs, projected::ProjectedCrs},
    cs::cartesian::{geocentric::XYZ, projected::ProjectedAxes},
    geodesy::GeodeticDatum,
    math::fp::Float,
    projections::{Projection, ProjectionError},
    quantities::{angle::Angle, length::Length, scale::PPM},
    units::length::M,
};
use nalgebra::{Matrix3, Vector3};

/// Coordinates transformation between two CRS via [normalized geocentric coordinates](XYZ).
///
/// A [CrsXYZTransformation] from CRS A to CRS B happens in 4 steps:
/// - convert coordinates from CRS A to [XYZ] coordinates in CRS A's datum using
///   a [ToXYZ] transformation provided by CRS A.
/// - transform the [XYZ] coordinates in CRS A's datum to a [XYZ] coordinates in a reference datum
///   using a [XYZTransformation].
/// - transform the [XYZ] coordinates in the reference datum to [XYZ] coordinates in CRS B's datum
///   using another [XYZTransformation].
/// - finally convert the [XYZ] coordinates in CRS B's datum into coordinates in CRS B using the
///   [ToXYZ] transformation provided by CRS B.
pub struct CrsXYZTransformation<'a> {
    src_to_xyz: Box<dyn ToXYZ + 'a>,
    src_xyz_to_ref_xyz: Box<dyn XYZTransformation + 'a>,
    dst_xyz_to_ref_xyz: Box<dyn XYZTransformation + 'a>,
    dst_to_xyz: Box<dyn ToXYZ + 'a>,
}

impl<'a> CrsXYZTransformation<'a> {
    pub fn new<S: ToXYZProvider, D: ToXYZProvider>(
        src: S,
        dst: D,
        src_to_ref: impl XYZTransformation + 'a,
        dst_to_ref: impl XYZTransformation + 'a,
    ) -> Self {
        CrsXYZTransformation {
            src_to_xyz: src.to_xyz_transformation(),
            src_xyz_to_ref_xyz: Box::new(src_to_ref),
            dst_xyz_to_ref_xyz: Box::new(dst_to_ref),
            dst_to_xyz: dst.to_xyz_transformation(),
        }
    }
}

impl<'a> CrsTransformation for CrsXYZTransformation<'a> {
    fn fwd(&self, src: &[Float], dst: &mut [Float]) -> Result<(), TransformationError> {
        let src_xyz = self.src_to_xyz.to_xyz(src)?;
        let ref_xyz = self.src_xyz_to_ref_xyz.to_ref(src_xyz);
        let dst_xyz = self.dst_xyz_to_ref_xyz.from_ref(ref_xyz);
        self.dst_to_xyz.from_xyz(dst_xyz, dst)
    }

    fn bwd(&self, dst: &[Float], src: &mut [Float]) -> Result<(), TransformationError> {
        let dst_xyz = self.dst_to_xyz.to_xyz(dst)?;
        let ref_xyz = self.dst_xyz_to_ref_xyz.to_ref(dst_xyz);
        let src_xyz = self.src_xyz_to_ref_xyz.from_ref(ref_xyz);
        self.src_to_xyz.from_xyz(src_xyz, src)
    }
}

/// Trait implemented by transformation from Crs coordinates to
/// [normalized geocentric coordinates](crate::cs::cartesian::geocentric::XYZ).
pub trait ToXYZ {
    fn to_xyz(&self, coords: &[Float]) -> Result<XYZ, TransformationError>;
    fn from_xyz(&self, xyz: XYZ, coords: &mut [Float]) -> Result<(), TransformationError>;
}

/// Trait implemented by [Crs](crate::crs::Crs) to provide a [ToXYZ] transformation.
pub trait ToXYZProvider {
    // NOTE: as of v1.84, this doesn't work for our use case:
    // fn to_xyz_transformation(&self) -> impl ToXYZTransformation;
    // `use<...>` precise capturing syntax is currently not allowed in return-position
    // `impl Trait` in traits
    // See issue #130044 <https://github.com/rust-lang/rust/issues/130044>
    fn to_xyz_transformation<'a>(&self) -> Box<dyn ToXYZ + 'a>;
}

impl ToXYZProvider for GeocentricCrs {
    fn to_xyz_transformation<'a>(&self) -> Box<dyn ToXYZ + 'a> {
        Box::new(self.clone())
    }
}

impl ToXYZ for GeocentricCrs {
    fn to_xyz(&self, coords: &[Float]) -> Result<XYZ, TransformationError> {
        Ok(self.axes.normalize(coords))
    }

    fn from_xyz(&self, xyz: XYZ, coords: &mut [Float]) -> Result<(), TransformationError> {
        Ok(self.axes.denormalize(xyz, coords))
    }
}

impl ToXYZProvider for GeographicCrs {
    fn to_xyz_transformation<'a>(&self) -> Box<dyn ToXYZ + 'a> {
        Box::new(self.clone())
    }
}

impl ToXYZ for GeographicCrs {
    fn to_xyz(&self, coords: &[Float]) -> Result<XYZ, TransformationError> {
        let llh = self.axes.normalize(coords);
        Ok(self.datum.llh_to_xyz(llh))
    }

    fn from_xyz(&self, xyz: XYZ, coords: &mut [Float]) -> Result<(), TransformationError> {
        let llh = self.datum.xyz_to_llh(xyz);
        Ok(self.axes.denormalize(llh, coords))
    }
}

impl ToXYZProvider for ProjectedCrs {
    fn to_xyz_transformation<'a>(&self) -> Box<dyn ToXYZ + 'a> {
        Box::new(ProjectedToXYZ {
            datum: self.datum.clone(),
            axes: self.axes,
            projection: self.projection.applied_to(self.datum.ellipsoid()),
        })
    }
}

struct ProjectedToXYZ<'a> {
    pub datum: GeodeticDatum,
    pub axes: ProjectedAxes,
    pub projection: Box<dyn Projection + 'a>,
}

impl From<ProjectionError> for TransformationError {
    fn from(_err: ProjectionError) -> Self {
        // TODO: do proper error conversion
        Self::OutOfBounds
    }
}

impl<'a> ToXYZ for ProjectedToXYZ<'a> {
    fn to_xyz(&self, coords: &[Float]) -> Result<XYZ, TransformationError> {
        let enh = self.axes.normalize(coords);
        let llh = self.projection.unproj(enh)?;
        Ok(self.datum.llh_to_xyz(llh))
    }

    fn from_xyz(&self, xyz: XYZ, coords: &mut [Float]) -> Result<(), TransformationError> {
        let llh = self.datum.xyz_to_llh(xyz);
        let enh = self.projection.proj(llh)?;
        Ok(self.axes.denormalize(enh, coords))
    }
}

/// A trait for transformations between ***normalized geocentric coordinates***.
pub trait XYZTransformation {
    fn to_ref(&self, src: XYZ) -> XYZ;
    fn from_ref(&self, dst: XYZ) -> XYZ;
}

/// The identity [XYZTransformation]
pub struct XYZIdentity;

impl XYZTransformation for XYZIdentity {
    fn to_ref(&self, src: XYZ) -> XYZ {
        src
    }

    fn from_ref(&self, dst: XYZ) -> XYZ {
        dst
    }
}

/// The [GeocentricTranslation] transforms **normalized geocentric coordinates** between two
/// geocentric CRS whose [GeodeticDatum] are related by a simple translation of the origin,
/// such that `xyz_dst = xyz_src + t`.
///
/// This epsg:1031.
#[derive(Debug, Clone, Copy)]
pub(crate) struct GeocentricTranslation {
    tx: Length,
    ty: Length,
    tz: Length,
}

impl GeocentricTranslation {
    /// Create a new instance with the given translation in metres.
    ///
    /// # Arguments
    ///
    /// - tx :the X coordinates of the origin of source CRS in the target CRS
    /// - ty: the Y coordinates of the origin of source CRS in the target CRS
    /// - tz: the Z coordinates of the origin of source CRS in the target CRS
    pub(crate) fn new(tx: Length, ty: Length, tz: Length) -> Self {
        GeocentricTranslation { tx, ty, tz }
    }
}

impl XYZTransformation for GeocentricTranslation {
    fn to_ref(&self, src: XYZ) -> XYZ {
        XYZ {
            x: src.x + self.tx,
            y: src.y + self.ty,
            z: src.z + self.tz,
        }
    }

    fn from_ref(&self, dst: XYZ) -> XYZ {
        XYZ {
            x: dst.x - self.tx,
            y: dst.y - self.ty,
            z: dst.z - self.tz,
        }
    }
}

#[derive(Debug, Default, Clone, Copy, PartialEq)]
pub enum RotationConvention {
    /// Rotation convention where positive angle rotations are clockwise when seen from
    /// the origin of the CS in the positive direction of the axis of rotation,
    /// or counter-clockwise when looking at the origin along the axis of rotation.
    #[default]
    PositionVector,
    /// Rotation convention opposite of the [PositionVector] convention.
    CoordinateFrame,
}

/// [Helmert7Params] transforms **normalized geocentric coordinates** between two geocentric CRS
/// whose [GeodeticDatum] are related by rotation around the origin, translation and scale change,
/// using the Helmert 7-parameters coordinates transformation
///
/// This is epsg:1033 (with PositionVector convention) or
/// epsg:1032 (with the CoordinateFrame convention).
#[derive(Debug, Clone, Copy)]
pub(crate) struct Helmert7Params {
    rot: Matrix3<Float>,
    inv_rot: Matrix3<Float>,
    t: Vector3<Float>,
    ppm: PPM,
}

impl Helmert7Params {
    /// Create a new [Helmert7Params] transformation.
    pub(crate) fn new(
        conv: RotationConvention,
        rotation: [Angle; 3],
        translation: [Length; 3],
        ppm: PPM,
    ) -> Self {
        let [rx, ry, rz] = rotation.map(|a: Angle| a.rad());
        let [tx, ty, tz] = translation.map(|l: Length| l.m());
        let rot = match conv {
            RotationConvention::PositionVector => {
                Matrix3::new(0.0, -rz, ry, rz, 0.0, -rx, -ry, rx, 0.0)
            }
            RotationConvention::CoordinateFrame => {
                Matrix3::new(0.0, rz, -ry, -rz, 0.0, rx, ry, -rx, 0.0)
            }
        };
        let t = Vector3::new(tx, ty, tz);
        Self {
            rot: Matrix3::identity() + rot,
            inv_rot: Matrix3::identity() - rot,
            t,
            ppm,
        }
    }
}

impl XYZTransformation for Helmert7Params {
    fn to_ref(&self, src: XYZ) -> XYZ {
        let s = Vector3::new(src.x.m(), src.y.m(), src.z.m());
        let x = self.rot * self.ppm.factor() * s + self.t;

        XYZ {
            x: x[0] * M,
            y: x[1] * M,
            z: x[2] * M,
        }
    }

    fn from_ref(&self, dst: XYZ) -> XYZ {
        let t = Vector3::new(dst.x.m(), dst.y.m(), dst.z.m());
        let s = self.inv_rot * self.ppm.inv_factor() * t - self.t;

        XYZ {
            x: s[0] * M,
            y: s[1] * M,
            z: s[2] * M,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::GeocentricTranslation;
    use super::XYZTransformation;
    use super::{Helmert7Params, RotationConvention};
    use crate::cs::cartesian::geocentric::XYZ;
    use crate::quantities::scale::PPM;
    use crate::units::angle::RAD;
    use crate::units::length::M;
    use approx::assert_relative_eq;

    #[test]
    fn geocentric_translation_fwd() {
        let t = GeocentricTranslation::new(84.87 * M, 96.49 * M, 116.95 * M);
        let src = XYZ {
            x: 3_771_793.97 * M,
            y: 140_253.34 * M,
            z: 5_124_304.35 * M,
        };
        let dst = t.to_ref(src);

        let expected_dst = XYZ {
            x: 3_771_878.84 * M,
            y: 140_349.83 * M,
            z: 5_124_421.30 * M,
        };
        assert_relative_eq!(dst.x.m(), expected_dst.x.m(), epsilon = 1e-6);
        assert_relative_eq!(dst.y.m(), expected_dst.y.m(), epsilon = 1e-6);
        assert_relative_eq!(dst.z.m(), expected_dst.z.m(), epsilon = 1e-6);
    }

    #[test]
    fn geocentric_translation_bwd() {
        let t = GeocentricTranslation::new(84.87 * M, 96.49 * M, 116.95 * M);
        let dst = XYZ {
            x: 3_771_878.84 * M,
            y: 140_349.83 * M,
            z: 5_124_421.30 * M,
        };
        let src = t.from_ref(dst);

        let expected_src = XYZ {
            x: 3_771_793.97 * M,
            y: 140_253.34 * M,
            z: 5_124_304.35 * M,
        };
        assert_relative_eq!(src.x.m(), expected_src.x.m(), epsilon = 1e-6);
        assert_relative_eq!(src.y.m(), expected_src.y.m(), epsilon = 1e-6);
        assert_relative_eq!(src.z.m(), expected_src.z.m(), epsilon = 1e-6);
    }

    #[test]
    fn helmert_fwd() {
        let t = Helmert7Params::new(
            RotationConvention::PositionVector,
            [0.0 * RAD, 0.0 * RAD, 2.685868e-6 * RAD],
            [0.0 * M, 0.0 * M, 4.5 * M],
            PPM(0.219),
        );
        let src = XYZ {
            x: 3_657_660.66 * M,
            y: 255_768.55 * M,
            z: 5_201_382.11 * M,
        };
        let dst = t.to_ref(src);

        let expected_dst = XYZ {
            x: 3_657_660.78 * M,
            y: 255_778.43 * M,
            z: 5_201_387.75 * M,
        };
        assert_relative_eq!(dst.x.m(), expected_dst.x.m(), epsilon = 1e-2);
        assert_relative_eq!(dst.y.m(), expected_dst.y.m(), epsilon = 1e-2);
        assert_relative_eq!(dst.z.m(), expected_dst.z.m(), epsilon = 1e-2);
    }

    #[test]
    fn helmert_bwd() {
        let t = Helmert7Params::new(
            RotationConvention::PositionVector,
            [0.0 * RAD, 0.0 * RAD, 2.685868e-6 * RAD],
            [0.0 * M, 0.0 * M, 4.5 * M],
            PPM(0.219),
        );
        let dst = XYZ {
            x: 3_657_660.78 * M,
            y: 255_778.43 * M,
            z: 5_201_387.75 * M,
        };
        let src = t.from_ref(dst);

        let expected_src = XYZ {
            x: 3_657_660.66 * M,
            y: 255_768.55 * M,
            z: 5_201_382.11 * M,
        };
        assert_relative_eq!(src.x.m(), expected_src.x.m(), epsilon = 1e-2);
        assert_relative_eq!(src.y.m(), expected_src.y.m(), epsilon = 1e-2);
        assert_relative_eq!(src.z.m(), expected_src.z.m(), epsilon = 1e-2);
    }
}
