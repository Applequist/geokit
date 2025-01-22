/// TODO: Make config dependent
pub type Float = f64;
pub const ZERO: Float = 0.0;
pub const PI_2: Float = std::f64::consts::FRAC_PI_2;
pub const PI_4: Float = std::f64::consts::FRAC_PI_4;
pub const PI: Float = std::f64::consts::PI;
pub const TAU: Float = std::f64::consts::TAU;

pub mod complex;
pub mod polynomial;

pub mod utils;
