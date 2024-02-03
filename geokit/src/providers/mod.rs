use crate::{
    crs::Crs,
    geodesy::{Ellipsoid, GeodeticDatum, PrimeMeridian},
    operation::Operation,
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
    fn crs(&self, id: &str) -> Option<Box<dyn Crs>>;
}

pub type TransformationPath = (Box<dyn Operation>, Box<dyn Operation>);
pub enum TransformationProviderError {
    NoTransformationPath,
}

pub trait TransformationProvider<S, D> {
    fn transformation(from: &S, to: &D) -> Result<TransformationPath, TransformationProviderError>;
}

pub mod default;
