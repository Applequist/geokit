use crate::cs::Coord;

/// [Pos] values represent **direct positions** via their [coordinates](crate::cs::Coord)
/// in a [Crs](crate::crs::Crs)'s coordinates system.
pub type Pos = Coord;

pub mod curve;
