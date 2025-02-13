use super::XYZTransformation;
use crate::{
    cs::cartesian::XYZ,
    math::fp::Float,
    quantities::{angle::Angle, length::Length, scale::PPM},
    units::length::M,
};
use nalgebra::{Matrix3, Vector3};

pub struct XYZIdentity;

impl XYZTransformation for XYZIdentity {
    fn src_to_dst(&self, src: XYZ) -> XYZ {
        src
    }

    fn dst_to_src(&self, dst: XYZ) -> XYZ {
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
    fn src_to_dst(&self, src: XYZ) -> XYZ {
        XYZ {
            x: src.x + self.tx,
            y: src.y + self.ty,
            z: src.z + self.tz,
        }
    }

    fn dst_to_src(&self, dst: XYZ) -> XYZ {
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
    fn src_to_dst(&self, src: XYZ) -> XYZ {
        let s = Vector3::new(src.x.m(), src.y.m(), src.z.m());
        let x = self.rot * self.ppm.factor() * s + self.t;

        XYZ {
            x: x[0] * M,
            y: x[1] * M,
            z: x[2] * M,
        }
    }

    fn dst_to_src(&self, dst: XYZ) -> XYZ {
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
    use approx::assert_relative_eq;

    use crate::cs::cartesian::XYZ;
    use crate::quantities::scale::PPM;
    use crate::transformations::XYZTransformation;
    use crate::units::angle::RAD;
    use crate::units::length::M;

    use super::GeocentricTranslation;
    use super::{Helmert7Params, RotationConvention};

    #[test]
    fn geocentric_translation_fwd() {
        let t = GeocentricTranslation::new(84.87 * M, 96.49 * M, 116.95 * M);
        let src = XYZ {
            x: 3_771_793.97 * M,
            y: 140_253.34 * M,
            z: 5_124_304.35 * M,
        };
        let dst = t.src_to_dst(src);

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
        let src = t.dst_to_src(dst);

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
        let dst = t.src_to_dst(src);

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
        let src = t.dst_to_src(dst);

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
