//! The cs modules provide types used to define and work with various **coordinates system** (or CS):
//! - system of axes, including their order, direction and units
//! - value types used as coordinates in these systems.

use crate::math::fp::Float;
use smallvec::SmallVec;

/// The type for (denormalized) tolerance.
type Tolerance = SmallVec<[Float; 3]>;

pub mod azimuth;
pub mod cartesian;
pub mod geodetic;
pub mod s1;
