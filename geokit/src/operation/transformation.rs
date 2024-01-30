use na::{Matrix3, Vector3};

use super::{DynOperation, Result};

/// The `GeocentricTranslation` transforms **normalized goecentric coordinates** between 2 [GeocentricCrs] whose
/// [GeodeticDatum] are related by a simple translation of the origin, such that
/// `xyz_dst = xyz_src + t`.
///
/// This epsg:1031.
#[derive(Debug, Clone, Copy)]
pub struct GeocentricTranslation {
    t: Vector3<f64>,
}

impl GeocentricTranslation {
    /// Create a new instance with the given translation in metres.
    ///
    /// # Arguments
    ///
    /// - tx :the X coordinates of the origin of source CRS in the target CRS
    /// - ty: the Y coordinates of the origin of source CRS in the target CRS
    /// - tz: the Z coordinates of the origin of source CRS in the target CRS
    pub fn new(tx: f64, ty: f64, tz: f64) -> Self {
        GeocentricTranslation {
            t: Vector3::new(tx, ty, tz),
        }
    }
}

impl DynOperation for GeocentricTranslation {
    fn fwd_in_dim(&self) -> usize {
        3
    }

    fn fwd_out_dim(&self) -> usize {
        3
    }

    fn fwd(&self, xyz_src: &[f64], xyz_dst: &mut [f64]) -> Result<()> {
        xyz_dst[0] = xyz_src[0] + self.t[0];
        xyz_dst[1] = xyz_src[1] + self.t[1];
        xyz_dst[2] = xyz_src[2] + self.t[2];
        Ok(())
    }

    fn bwd(&self, xyz_dst: &[f64], xyz_src: &mut [f64]) -> Result<()> {
        xyz_src[0] = xyz_dst[0] - self.t[0];
        xyz_src[1] = xyz_dst[1] - self.t[1];
        xyz_src[2] = xyz_dst[2] - self.t[2];
        Ok(())
    }
}

#[derive(Debug, Default, Clone, Copy)]
pub enum RotationConvention {
    /// Rotation convention where positive angle rotations are clockwise when seen from
    /// the origin of the CS in the positive direction of the axis of rotation,
    /// or counter-clockwise when looking at the origin along the axis of rotation.
    #[default]
    PositionVector,
    /// Rotation convention opposite of the [PositionVector] convention.
    CoordinateFrame,
}

/// [Helmert7Params] transforms **normalized geocentric coordinates** between 2 [GeocentricCrs]
/// whose [GeodeticDatum] are related by rotation around the origin, translation and scale change,
/// using the Helmert 7-parameters coordinates transformation
///
/// This is epsg:1033 (with PositionVector convention) or
/// epsg:1032 (with the CoordinateFrame convention).
#[derive(Debug, Clone, Copy)]
pub struct Helmert7Params {
    rot: Matrix3<f64>,
    inv_rot: Matrix3<f64>,
    t: Vector3<f64>,
    ds_ppm: f64,
}

impl Helmert7Params {
    /// Create a new [Helmert7Params] transformation.
    pub fn new(
        (conv, rx, ry, rz): (RotationConvention, f64, f64, f64),
        (tx, ty, tz): (f64, f64, f64),
        ds_ppm: f64,
    ) -> Self {
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
            ds_ppm,
        }
    }
}

impl DynOperation for Helmert7Params {
    fn fwd_in_dim(&self) -> usize {
        3
    }

    fn fwd_out_dim(&self) -> usize {
        3
    }

    fn fwd(&self, input: &[f64], output: &mut [f64]) -> Result<()> {
        let s = Vector3::new(input[0], input[1], input[2]);
        let x = self.rot * (1.0 + self.ds_ppm * 1e-6) * s + self.t;
        output.copy_from_slice(x.as_ref());
        Ok(())
    }

    fn bwd(&self, input: &[f64], output: &mut [f64]) -> Result<()> {
        let t = Vector3::new(input[0], input[1], input[2]);
        let s = self.inv_rot * (1.0 - self.ds_ppm * 1e-6) * t - self.t;
        output.copy_from_slice(s.as_ref());
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use approx::assert_relative_eq;

    use crate::operation::DynOperation;

    use super::{GeocentricTranslation, Helmert7Params, RotationConvention};

    #[test]
    fn geocentric_translation_fwd() {
        let t = GeocentricTranslation::new(84.87, 96.49, 116.95);
        let source_xyz = [3_771_793.97, 140_253.34, 5_124_304.35];
        let mut xyz = [0.; 3];
        let expected_xyz = [3_771_878.84, 140_349.83, 5_124_421.30];
        t.fwd(&source_xyz, &mut xyz).unwrap();
        assert_relative_eq!(xyz[0], expected_xyz[0], epsilon = 1e-6);
        assert_relative_eq!(xyz[1], expected_xyz[1], epsilon = 1e-6);
        assert_relative_eq!(xyz[2], expected_xyz[2], epsilon = 1e-6);
    }

    #[test]
    fn geocentric_translation_bwd() {
        let t = GeocentricTranslation::new(84.87, 96.49, 116.95);
        let xs = [3_771_793.97, 140_253.34, 5_124_304.35];
        let mut xyz = [0.; 3];
        let xt = [3_771_878.84, 140_349.83, 5_124_421.30];
        t.bwd(&xt, &mut xyz).unwrap();
        assert_relative_eq!(xyz[0], xs[0], epsilon = 1e-6);
        assert_relative_eq!(xyz[1], xs[1], epsilon = 1e-6);
        assert_relative_eq!(xyz[2], xs[2], epsilon = 1e-6);
    }

    #[test]
    fn helmert_fwd() {
        let t = Helmert7Params::new(
            (RotationConvention::PositionVector, 0.0, 0.0, 2.685868e-6),
            (0.0, 0.0, 4.5),
            0.219,
        );
        let xs = [3_657_660.66, 255_768.55, 5_201_382.11];
        let mut xyz = [0.; 3];
        let xt = [3_657_660.78, 255_778.43, 5_201_387.75];
        t.fwd(&xs, &mut xyz).unwrap();
        assert_relative_eq!(xyz[0], xt[0], epsilon = 1e-2);
        assert_relative_eq!(xyz[1], xt[1], epsilon = 1e-2);
        assert_relative_eq!(xyz[2], xt[2], epsilon = 1e-2);
    }

    #[test]
    fn helmert_bwd() {
        let t = Helmert7Params::new(
            (RotationConvention::PositionVector, 0.0, 0.0, 2.685868e-6),
            (0.0, 0.0, 4.5),
            0.219,
        );
        let xs = [3_657_660.66, 255_768.55, 5_201_382.11];
        let mut xyz = [0.; 3];
        let xt = [3_657_660.78, 255_778.43, 5_201_387.75];
        t.bwd(&xt, &mut xyz).unwrap();
        assert_relative_eq!(xyz[0], xs[0], epsilon = 1e-2);
        assert_relative_eq!(xyz[1], xs[1], epsilon = 1e-2);
        assert_relative_eq!(xyz[2], xs[2], epsilon = 1e-2);
    }
}
