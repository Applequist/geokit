use std::f64::consts::PI;

/// An angle unit as a to-radians angle converter.
/// An angle expressed in this unit is converted to radians as follow:
/// ```
/// let deg = AngleUnit(PI, 180.0);
/// let angle_deg = 1.0;
/// leet angle_rad = angle_deg * deg.0 / deg.1;
/// ```
#[derive(Debug, Copy, Clone, PartialEq)]
pub struct AngleUnit(pub f64, pub f64);

impl AngleUnit {
    #[inline]
    pub const fn rad_per_unit(&self) -> f64 {
        self.0 / self.1
    }
}

pub const RAD: AngleUnit = AngleUnit(1.0, 1.0);
pub const DEG: AngleUnit = AngleUnit(PI, 180.0);
pub const GRAD: AngleUnit = AngleUnit(PI, 200.0);
pub const SEC: AngleUnit = AngleUnit(PI, 648_000.0);
