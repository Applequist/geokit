pub type Coord3D = (f64, f64, f64);

/// `CoordSpace` determines a partial order on CRS.
/// 2 `CoordSpace`s are comparable if there is a **coordinates conversion path** between each other.
/// Which one is less than the other has to do with **transformation path**: coordinates are
/// usually converted to the lowest CRS kind comparable to both source and destination CRS, then a
/// datum transformation is performed if necessary.
#[derive(Debug, Copy, Clone, PartialEq)]
pub enum CoordSpace {
    Geocentric,
    Geodetic,
    Topocentric,
    Projected,
}

impl PartialOrd for CoordSpace {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        match (self, other) {
            (CoordSpace::Geocentric, CoordSpace::Geocentric) => Some(std::cmp::Ordering::Equal),
            (CoordSpace::Geocentric, _) => Some(std::cmp::Ordering::Less),

            (CoordSpace::Geodetic, CoordSpace::Geocentric) => Some(std::cmp::Ordering::Greater),
            (CoordSpace::Geodetic, CoordSpace::Geodetic) => Some(std::cmp::Ordering::Equal),
            (CoordSpace::Geodetic, CoordSpace::Projected) => Some(std::cmp::Ordering::Less),
            (CoordSpace::Geodetic, _) => None, // not comparable to Topocentric or Geodetic2D

            (CoordSpace::Topocentric, CoordSpace::Geocentric) => Some(std::cmp::Ordering::Greater),
            (CoordSpace::Topocentric, CoordSpace::Topocentric) => Some(std::cmp::Ordering::Equal),
            (CoordSpace::Topocentric, _) => None,

            (CoordSpace::Projected, CoordSpace::Geocentric) => Some(std::cmp::Ordering::Greater),
            (CoordSpace::Projected, CoordSpace::Geodetic) => Some(std::cmp::Ordering::Greater),
            (CoordSpace::Projected, CoordSpace::Projected) => Some(std::cmp::Ordering::Equal),
            (CoordSpace::Projected, _) => None,
        }
    }
}
