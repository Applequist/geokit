use std::fmt::Debug;

use dyn_clone::DynClone;
use geomath::linalg::{matrix::Mat3, vector::Vec3};
use num::One;
use smallvec::SmallVec;
use thiserror::Error;

use crate::geodesy::{Ellipsoid, GeodeticDatum, PrimeMeridian};

/// Value returned when [`Transformation::apply()`] cannot be completed.
/// `index` indicates the index of **coordinates** (as opposed to the index of individual floats)
/// in the transformation input where the error occured.
/// `reason` is a human readable reason.
/// TODO: Do we want an enum for reason instead ?
#[derive(Error, Debug)]
#[error("Transformation failed at input index {}: {}", index, reason)]
pub struct TransformationError {
    index: usize,
    reason: String,
}

/// The return type for [`Transformation::apply()`] and [`Transformation::apply_seq()`].
pub type Result<T> = std::result::Result<T, TransformationError>;

/// Base trait for coordinate transformations.
pub trait Transformation: DynClone {
    fn in_dim(&self) -> usize;

    fn out_dim(&self) -> usize;

    /// Perform the forward transformation on the given input.
    fn fwd(&self, input: &[f64], output: &mut [f64]) -> Result<()>;

    /// Apply the forward transformation on coordinates in the `input` slice, storing the transformed
    /// coordinates into the `output` slice, and returning the number of successful applications when either:
    /// * the input slice has been fully processed
    /// * or the output slice has been filled,
    ///
    /// # Errors
    ///
    /// If an error occurs while transforming the sequence, this function returns an [`TransformationError`] containing:
    /// * the index of the transformation application, eg 0 on 1st application, 1 on 2nd...
    /// * the reason for the error
    ///
    /// If the index is greater that 0, the previous transformation results are avalaible in `output`.
    fn fwd_seq(&self, input: &[f64], output: &mut [f64]) -> Result<usize> {
        let input_chunks = input.chunks_exact(self.in_dim());
        let output_chunks = output.chunks_exact_mut(self.out_dim());
        let mut index: usize = 0;
        for (i, o) in input_chunks.zip(output_chunks) {
            match self.fwd(i, o) {
                Ok(_) => index += 1,
                Err(TransformationError { index: _, reason }) => {
                    return Err(TransformationError { index, reason })
                }
            };
        }
        Ok(index)
    }

    fn bwd(&self, input: &[f64], output: &mut [f64]) -> Result<()>;

    /// Apply the backward transformation on coordinates in the `input` slice, storing the transformed
    /// coordinates into the `output` slice, and returning the number of successful applications when either:
    /// * the input slice has been fully processed
    /// * or the output slice has been filled,
    ///
    /// # Errors
    ///
    /// If an error occurs while transforming the sequence, this function returns an [`TransformationError`] containing:
    /// * the index of the transformation application, eg 0 on 1st application, 1 on 2nd...
    /// * the reason for the error
    ///
    /// If the index is greater that 0, the previous transformation results are avalaible in `output`.
    fn bwd_seq(&self, input: &[f64], output: &mut [f64]) -> Result<usize> {
        let input_chunks = input.chunks_exact(self.out_dim());
        let output_chunks = output.chunks_exact_mut(self.in_dim());
        let mut index: usize = 0;
        for (i, o) in input_chunks.zip(output_chunks) {
            match self.bwd(i, o) {
                Ok(_) => index += 1,
                Err(TransformationError { index: _, reason }) => {
                    return Err(TransformationError { index, reason })
                }
            };
        }
        Ok(index)
    }

    /// Chain this transformation with the given one.
    /// It returns a new transformation that perform this transformation,
    /// followed by the given one.
    fn and_then<T>(self, t: T) -> Chain<Self, T>
    where
        T: Transformation,
        Self: Sized,
    {
        Chain {
            first: self,
            then: t,
        }
    }
}

dyn_clone::clone_trait_object!(Transformation);

/// A chained transformation.
pub struct Chain<A, B> {
    first: A,
    then: B,
}

impl<A, B> Debug for Chain<A, B>
where
    A: Debug,
    B: Debug,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Chain")
            .field("first", &self.first)
            .field("then", &self.then)
            .finish()
    }
}

impl<A, B> Clone for Chain<A, B>
where
    A: Clone,
    B: Clone,
{
    fn clone(&self) -> Self {
        Self {
            first: self.first.clone(),
            then: self.then.clone(),
        }
    }
}

impl<A, B> Transformation for Chain<A, B>
where
    A: Transformation + Clone,
    B: Transformation + Clone,
{
    fn in_dim(&self) -> usize {
        self.first.in_dim()
    }

    fn out_dim(&self) -> usize {
        self.then.out_dim()
    }

    fn fwd(&self, i: &[f64], o: &mut [f64]) -> Result<()> {
        let mut os: SmallVec<[f64; 3]> = SmallVec::from_elem(0.0, self.first.out_dim());
        self.first.fwd(i, &mut os)?;
        self.then.fwd(&os, o)
    }

    fn bwd(&self, i: &[f64], o: &mut [f64]) -> Result<()> {
        let mut os: SmallVec<[f64; 3]> = SmallVec::from_elem(0.0, self.first.out_dim());
        self.then.bwd(i, &mut os)?;
        self.first.bwd(&os, o)
    }
}

/// [Inv] inverses a transformation. Its [fwd] transformation performs the inner's [bwd]
/// tranformation and vice-versa.
pub struct Inv<T> {
    inner: T,
}

impl<T> Debug for Inv<T>
where
    T: Debug,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Inv").field("inner", &self.inner).finish()
    }
}

impl<T> Clone for Inv<T>
where
    T: Clone,
{
    fn clone(&self) -> Self {
        Self {
            inner: self.inner.clone(),
        }
    }
}

impl<I> Transformation for Inv<I>
where
    I: Transformation + Clone,
{
    fn in_dim(&self) -> usize {
        self.inner.out_dim()
    }

    fn out_dim(&self) -> usize {
        self.inner.in_dim()
    }

    fn fwd(&self, i: &[f64], o: &mut [f64]) -> Result<()> {
        self.inner.bwd(i, o)
    }

    fn bwd(&self, i: &[f64], o: &mut [f64]) -> Result<()> {
        self.inner.fwd(i, o)
    }
}

/// Return the inverse of the given transformation, eg a transformation whose [fwd] transformation
/// performs the [bwd] transformation of the given one.
pub fn inv<T: Transformation>(t: T) -> Inv<T> {
    Inv { inner: t }
}

// impl Transformation for Box<dyn Transformation> {
//     fn in_dim(&self) -> usize {
//         self.as_ref().in_dim()
//     }
//
//     fn out_dim(&self) -> usize {
//         self.as_ref().out_dim()
//     }
//
//     fn fwd(&self, i: &[f64], o: &mut [f64]) -> Result<()> {
//         self.as_ref().fwd(i, o)
//     }
//
//     fn bwd(&self, i: &[f64], o: &mut [f64]) -> Result<()> {
//         self.as_ref().bwd(i, o)
//     }
// }

#[derive(Debug, Clone)]
struct Identity<const N: usize>;

impl<const N: usize> Transformation for Identity<N> {
    fn in_dim(&self) -> usize {
        N
    }

    fn out_dim(&self) -> usize {
        N
    }

    fn fwd(&self, i: &[f64], o: &mut [f64]) -> Result<()> {
        o.copy_from_slice(i);
        Ok(())
    }

    fn bwd(&self, i: &[f64], o: &mut [f64]) -> Result<()> {
        o.copy_from_slice(i);
        Ok(())
    }
}

pub fn identity<const N: usize>() -> impl Transformation {
    Identity::<N>
}

pub type ToOrd = (usize, f64);

/// [CoordScaling] is used to normalize/denormalize coordinates by:
/// - reordering coordinates
/// - converting coordinates units
/// - zeroing input coordinates in excess or output coordinates in excess.
///
/// Normalized geocentric coordinates are (x, y, z) measured in metres.
/// Normalized geographic coordinates are (lon, lat(, height)) measured in radians (and metres).
#[derive(Debug, Clone)]
pub struct CoordScaling(pub Vec<ToOrd>);

impl Transformation for CoordScaling {
    fn in_dim(&self) -> usize {
        self.0.len()
    }

    fn out_dim(&self) -> usize {
        3
    }

    fn fwd(&self, input: &[f64], output: &mut [f64]) -> Result<()> {
        output.copy_from_slice(&[0.; 3]);
        for (nix, (ix, f)) in self.0.iter().enumerate() {
            output[nix] = input[*ix] * f;
        }
        Ok(())
    }

    fn bwd(&self, input: &[f64], output: &mut [f64]) -> Result<()> {
        for (nix, (ix, f)) in self.0.iter().enumerate() {
            output[*ix] = input[nix] / f;
        }
        Ok(())
    }
}

/// `GeogToGeoc` converts **normalized geographic coordinates** to **normalized geocentric coordinates**
/// between a [GeographicCrs] and a [GeocentricCrs] that share the same [GeodeticDatum].
///
/// As mentioned in the EPSG guidance note 7 part 2, paragraph 4.1.1, this transformation
/// first transform to/from non-greenwich base geographic coordinates into greenwich-base ones
/// before/after transformation to/from geocentric coordinates.
///
/// This is epsg:9602.
#[derive(Debug, Clone, Copy)]
pub struct GeodToGeoc {
    ellipsoid: Ellipsoid,
    prime_meridian: PrimeMeridian,
}

impl GeodToGeoc {
    pub fn new(datum: &GeodeticDatum) -> Self {
        Self {
            ellipsoid: datum.ellipsoid(),
            prime_meridian: datum.prime_meridian(),
        }
    }
}

impl Transformation for GeodToGeoc {
    fn in_dim(&self) -> usize {
        3
    }

    fn out_dim(&self) -> usize {
        3
    }

    /// Convert geographic coordinates to geocentric coordinates.
    fn fwd(&self, input: &[f64], output: &mut [f64]) -> Result<()> {
        let mut llh_cpy = [0.0; 3];
        llh_cpy.copy_from_slice(input);
        llh_cpy[0] = self.prime_meridian.convert_lon_to_gw(input[0]);
        self.ellipsoid.llh_to_xyz(&llh_cpy, output);
        Ok(())
    }

    /// Convert geocentric coordinates to geographic coordinates.
    fn bwd(&self, input: &[f64], output: &mut [f64]) -> Result<()> {
        self.ellipsoid.xyz_to_llh(input, output);
        output[0] = self.prime_meridian.convert_lon_from_gw(output[0]);
        Ok(())
    }
}

/// The `DatumShift` transforms **normalized goecentric coordinates** between 2 [GeocentricCrs] whose
/// [GeodeticDatum] are related by a simple translation of the origin.
///
/// This epsg:1031.
#[derive(Debug, Clone, Copy)]
pub struct DatumShift {
    t: Vec3,
}

impl DatumShift {
    /// Create a new instance with the given translation in metres.
    ///
    /// # Arguments
    ///
    /// - tx :the X coordinates of the origin of source CRS in the target CRS
    /// - ty: the Y coordinates of the origin of source CRS in the target CRS
    /// - tz: the Z coordinates of the origin of source CRS in the target CRS
    pub fn new(tx: f64, ty: f64, tz: f64) -> Self {
        DatumShift {
            t: Vec3::new([tx, ty, tz]),
        }
    }
}

impl Transformation for DatumShift {
    fn in_dim(&self) -> usize {
        3
    }

    fn out_dim(&self) -> usize {
        3
    }

    fn fwd(&self, xyz: &[f64], s_xyz: &mut [f64]) -> Result<()> {
        s_xyz[0] = xyz[0] + self.t[0];
        s_xyz[1] = xyz[1] + self.t[1];
        s_xyz[2] = xyz[2] + self.t[2];
        Ok(())
    }

    fn bwd(&self, s_xyz: &[f64], xyz: &mut [f64]) -> Result<()> {
        xyz[0] = s_xyz[0] - self.t[0];
        xyz[1] = s_xyz[1] - self.t[1];
        xyz[2] = s_xyz[2] - self.t[2];
        Ok(())
    }
}

#[derive(Debug, Clone, Copy)]
pub enum RotationConvention {
    /// Rotation convention where positive angle rotations are clockwise when seen from
    /// the origin of the CS in the positive direction of the axis of rotation,
    /// or counter-clockwise when looking at the origin along the axis of rotation.
    PositionVector,
    /// Rotation convention opposite of the [PositionVector] convention.
    CoordinateFrame,
}

impl Default for RotationConvention {
    fn default() -> Self {
        RotationConvention::PositionVector
    }
}

/// [Helmert7Params] transforms **normalized geocentric coordinates** between 2 [GeocentricCrs]
/// whose [GeodeticDatum] are related by rotation around the origin, translation and scale change,
/// using the Helmert 7-parameters coordinates transformation
///
/// This is epsg:1033 (with PositionVector convention) or
/// epsg:1032 (with the CoordinateFrame convention).
#[derive(Debug, Clone, Copy)]
pub struct Helmert7Params {
    rot: Mat3,
    inv_rot: Mat3,
    t: Vec3,
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
            RotationConvention::PositionVector => Mat3::with_rows([
                Vec3::new([0.0, -rz, ry]),
                Vec3::new([rz, 0.0, -rx]),
                Vec3::new([-ry, rx, 0.0]),
            ]),
            RotationConvention::CoordinateFrame => Mat3::with_rows([
                Vec3::new([0.0, rz, -ry]),
                Vec3::new([-rz, 0.0, rx]),
                Vec3::new([ry, -rx, 0.0]),
            ]),
        };
        let t = Vec3::new([tx, ty, tz]);
        Self {
            rot: Mat3::one() + rot,
            inv_rot: Mat3::one() - rot,
            t,
            ds_ppm,
        }
    }
}

impl Transformation for Helmert7Params {
    fn in_dim(&self) -> usize {
        3
    }

    fn out_dim(&self) -> usize {
        3
    }

    fn fwd(&self, input: &[f64], output: &mut [f64]) -> Result<()> {
        let s = Vec3::new([input[0], input[1], input[2]]);
        let x = self.rot * (1.0 + self.ds_ppm * 1e-6) * s + self.t;
        output.copy_from_slice(&x);
        Ok(())
    }

    fn bwd(&self, input: &[f64], output: &mut [f64]) -> Result<()> {
        let t: Vec3 = Vec3::new([input[0], input[1], input[2]]);
        let s = self.inv_rot * (1.0 - self.ds_ppm * 1e-6) * t - self.t;
        output.copy_from_slice(&s);
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use crate::crs::geodetic::{GeodeticAxes, GeodeticCrs};
    use crate::id::Id;
    use approx::assert_relative_eq;

    use crate::geodesy::GeodeticDatum;
    use crate::transformation::{DatumShift, Helmert7Params, RotationConvention, Transformation};

    #[test]
    fn normalization() {
        let latlondeg = GeodeticCrs::new(
            Id::name("WGS84"),
            GeodeticDatum::default(),
            GeodeticAxes::NorthWest {
                angle_unit: 1.0_f64.to_radians(),
            },
        );
        let t = latlondeg.axes().normalization();
        assert_eq!(t.in_dim(), 2);
        assert_eq!(t.out_dim(), 3);

        let i = [10.0, 110.0];
        let mut o = [0.0; 3];
        t.fwd(&i, &mut o).unwrap();
        assert_eq!(o, [-110.0f64.to_radians(), 10.0f64.to_radians(), 0.0]);
    }

    #[test]
    fn geog_to_geoc_fwd() {
        unimplemented!()
    }

    #[test]
    fn geog_to_geoc_bwd() {
        unimplemented!()
    }

    #[test]
    fn datum_shift_fwd() {
        let t = DatumShift::new(84.87, 96.49, 116.95);
        let source_xyz = [3_771_793.97, 140_253.34, 5_124_304.35];
        let mut xyz = [0.; 3];
        let expected_xyz = [3_771_878.84, 140_349.83, 5_124_421.30];
        t.fwd(&source_xyz, &mut xyz).unwrap();
        assert_relative_eq!(xyz[0], expected_xyz[0], epsilon = 1e-6);
        assert_relative_eq!(xyz[1], expected_xyz[1], epsilon = 1e-6);
        assert_relative_eq!(xyz[2], expected_xyz[2], epsilon = 1e-6);
    }

    #[test]
    fn datum_shift_bwd() {
        let t = DatumShift::new(84.87, 96.49, 116.95);
        let xs = [3_771_793.97, 140_253.34, 5_124_304.35];
        let mut xyz = [0.; 3];
        let xt = [3_771_878.84, 140_349.83, 5_124_421.30];
        t.bwd(&xt, &mut xyz).unwrap();
        assert_relative_eq!(xyz[0], xs[0], epsilon = 1e-6);
        assert_relative_eq!(xyz[1], xs[1], epsilon = 1e-6);
        assert_relative_eq!(xyz[2], xs[2], epsilon = 1e-6);
    }

    #[test]
    fn test_helmert_fwd() {
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
    fn test_helmert_bwd() {
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
