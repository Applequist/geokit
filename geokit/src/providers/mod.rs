use crate::{
    crs::Crs,
    geodesy::{Ellipsoid, GeodeticDatum, PrimeMeridian},
    operation::{Chain, Operation},
};

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

pub type TransformationPath = (
    Chain<Box<dyn Operation>, Box<dyn Operation>>,
    Chain<Box<dyn Operation>, Box<dyn Operation>>,
);

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

impl TransformationProvider for DefaultTransformationProvider {
    fn transformation(
        &self,
        src: &Crs,
        dst: &Crs,
    ) -> Result<TransformationPath, TransformationProviderError> {
        if src.ref_datum_id() != dst.ref_datum_id() {
            return Err(TransformationProviderError::NoTransformationPath);
        }
        let (src_to_ref, ref_to_src) = src.to_wgs84_geoc();
        let (dst_to_ref, ref_to_dst) = dst.to_wgs84_geoc();
        Ok((
            src_to_ref.and_then(ref_to_dst),
            dst_to_ref.and_then(ref_to_src),
        ))
    }
}
