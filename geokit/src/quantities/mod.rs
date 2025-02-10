use crate::math::Float;

pub trait Convertible {
    type Unit;
    fn val(self, unit: Self::Unit) -> Float;
}

pub mod angle;
pub mod length;
pub mod scale;
