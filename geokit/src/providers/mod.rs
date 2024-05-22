use crate::{
    crs::Crs,
    geodesy::{Ellipsoid, GeodeticDatum, PrimeMeridian},
    operation::Operation,
};
use crate::crs::{GeocentricAxes, GeodeticAxes, ProjectedAxes, ProjectionSpec};
use crate::operation::{identity, Inv};
use crate::operation::conversion::{GeogToGeoc, Normalization};

pub trait GeodesyProvider {
    fn ellipsoid_ids(&self) -> impl Iterator<Item = &str>;
    fn ellipsoid(&self, id: &str) -> Option<Ellipsoid>;

    fn prime_meridian_ids(&self) -> impl Iterator<Item = &str>;
    fn prime_meridian(&self, id: &str) -> Option<PrimeMeridian>;

    fn datum_ids(&self) -> impl Iterator<Item = &str>;
    fn datum(&self, id: &str) -> Option<GeodeticDatum>;
}

pub trait CrsProvider {
    fn crs_ids(&self) -> impl Iterator<Item = &str>;
    fn crs(&self, id: &str) -> Option<Crs>;
}

pub type TransformationPath = (Box<dyn Operation>, Box<dyn Operation>);

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum TransformationProviderError {
    NoTransformationPath,
}

pub trait TransformationProvider {
    fn transformation(
        &self,
        src: &Crs,
        dst: &Crs,
    ) -> Result<TransformationPath, TransformationProviderError>;
}

pub struct DefaultTransformationProvider;

impl DefaultTransformationProvider {
    pub fn conversion(&self, src: &Crs, dst: &Crs) -> TransformationPath {
        debug_assert_eq!(src.datum(), dst.datum());
        match src {
            Crs::Geocentric {
                datum,
                axes: src_axes,
                ..
            } => self.convert_from_geocentric(datum, src_axes, dst),
            Crs::Geographic {
                datum,
                axes: src_axes,
                ..
            } => self.convert_from_geographic(datum, src_axes, dst),
            Crs::Projected {
                datum,
                axes: src_axes,
                projection: src_projection,
                ..
            } => self.convert_from_projected(datum, src_axes, src_projection, dst),
        }
    }

    fn convert_from_geocentric(
        &self,
        datum: &GeodeticDatum,
        _src_axes: &GeocentricAxes,
        dst: &Crs,
    ) -> TransformationPath {
        match dst {
            Crs::Geocentric { .. } => (identity(3, 3).boxed(), identity(3, 3).boxed()),
            Crs::Geographic { axes: dst_axes, .. } => {
                let dst_to_src = Normalization::from(*dst_axes).and_then(GeogToGeoc::new(datum));
                let src_to_dst = Inv(dst_to_src.clone());
                (src_to_dst.boxed(), dst_to_src.boxed())
            }
            Crs::Projected {
                axes: dst_axes,
                projection: dst_projection,
                ..
            } => {
                let dst_to_src = Normalization::from(*dst_axes)
                    .and_then(Inv(dst_projection.projection(datum.ellipsoid())))
                    .and_then(GeogToGeoc::new(datum));
                // skipped identity geoc normalization
                let src_to_dst = Inv(dst_to_src.clone());
                (src_to_dst.boxed(), dst_to_src.boxed())
            }
        }
    }

    fn convert_from_geographic(
        &self,
        datum: &GeodeticDatum,
        src_axes: &GeodeticAxes,
        dst: &Crs,
    ) -> TransformationPath {
        match dst {
            Crs::Geocentric { .. } => {
                let src_to_dst = Normalization::from(*src_axes).and_then(GeogToGeoc::new(datum));
                let dst_to_src = Inv(src_to_dst.clone());
                (src_to_dst.boxed(), dst_to_src.boxed())
            }
            Crs::Geographic { axes: dst_axes, .. } => {
                // Same datum
                let src_to_dst =
                    Normalization::from(*src_axes).and_then(Inv(Normalization::from(*dst_axes)));
                let dst_to_src = Inv(src_to_dst.clone());
                (src_to_dst.boxed(), dst_to_src.boxed())
            }
            Crs::Projected {
                axes: dst_axes,
                projection: dst_projection,
                ..
            } => {
                let src_to_dst = Normalization::from(*src_axes)
                    .and_then(dst_projection.projection(datum.ellipsoid()))
                    .and_then(Inv(Normalization::from(*dst_axes)));
                let dst_to_src = Inv(src_to_dst.clone());
                (src_to_dst.boxed(), dst_to_src.boxed())
            }
        }
    }

    fn convert_from_projected(
        &self,
        datum: &GeodeticDatum,
        src_axes: &ProjectedAxes,
        src_projection: &ProjectionSpec,
        dst: &Crs,
    ) -> TransformationPath {
        match dst {
            Crs::Geocentric { .. } => {
                let src_to_dst = Normalization::from(*src_axes)
                    .and_then(Inv(src_projection.projection(datum.ellipsoid())))
                    .and_then(GeogToGeoc::new(datum));
                let dst_to_src = Inv(src_to_dst.clone());
                (src_to_dst.boxed(), dst_to_src.boxed())
            }
            Crs::Geographic { axes: dst_axes, .. } => {
                let dst_to_src = Normalization::from(*dst_axes)
                    .and_then(src_projection.projection(datum.ellipsoid()))
                    .and_then(Inv(Normalization::from(*src_axes)));
                let src_to_dst = Inv(dst_to_src.clone());
                (src_to_dst.boxed(), dst_to_src.boxed())
            }
            Crs::Projected {
                axes: dst_axes,
                projection: dst_projection,
                ..
            } => {
                let src_to_dst = Normalization::from(*src_axes)
                    .and_then(Inv(src_projection.projection(datum.ellipsoid())))
                    .and_then(dst_projection.projection(datum.ellipsoid()))
                    .and_then(Inv(Normalization::from(*dst_axes)));
                let dst_to_src = Inv(src_to_dst.clone());
                (src_to_dst.boxed(), dst_to_src.boxed())
            }
        }
    }
}
impl TransformationProvider for DefaultTransformationProvider {
    fn transformation(
        &self,
        src: &Crs,
        dst: &Crs,
    ) -> Result<TransformationPath, TransformationProviderError> {
        if src.datum() == dst.datum() {
            return Ok(self.conversion(src, dst));
        }
        if src.ref_datum_id() != dst.ref_datum_id() {
            return Err(TransformationProviderError::NoTransformationPath);
        }
        let (src_to_ref, ref_to_src) = src.to_ref_geoc();
        let (dst_to_ref, ref_to_dst) = dst.to_ref_geoc();
        Ok((
            src_to_ref.and_then(ref_to_dst).boxed(),
            dst_to_ref.and_then(ref_to_src).boxed(),
        ))
    }
}
