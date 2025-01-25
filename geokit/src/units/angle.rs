use crate::math::{Float, PI};

/// An angle unit as a to-radians angle converter.
/// An angle expressed in this unit is converted to radians as follow:
/// ```
/// use geokit::units::angle::AngleUnit;
/// use geokit::math::PI;
/// let deg = AngleUnit(PI, 180.0);
/// let angle_deg = 1.0;
/// let angle_rad = angle_deg * deg.rad_per_unit();
/// ```
#[derive(Debug, Copy, Clone, PartialEq)]
pub struct AngleUnit(pub Float, pub Float);

impl AngleUnit {
    #[inline]
    pub const fn rad_per_unit(&self) -> Float {
        self.0 / self.1
    }
}

pub const RAD: AngleUnit = AngleUnit(1.0, 1.0);
pub const DEG: AngleUnit = AngleUnit(PI, 180.0);
pub const GRAD: AngleUnit = AngleUnit(PI, 200.0);
pub const SEC: AngleUnit = AngleUnit(PI, 648_000.0);
