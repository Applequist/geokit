use crate::geodesy::GeodeticDatum;
use crate::id::Id;
use crate::transformation::Transformation;

/// A [`CoordSpace`] determines a class of coordinates, eg geocentric, geodetic...
/// Each [`CoordSpace`] has **normalized coordinates**:
/// * [`CoordSpace::Geocentric`] coordinates are normalized: x, y and z are expressed **in
/// meters**,
/// * Normalized [`CoordSpace::Geodetic`] coordinates are:
///   1. longitude in radians positive east **of the Greenwich prime meridian**,
///   2. latitude in radians positive north of the equatorial plane,
///   3. ellipsoidal height in meters positive up.
/// * Normalized [`CoordSpace::Topocentric`] coordinates are:
///   1. u in meters positive east of the origin
///   2. v in meters positive north of the origin
///   3. z in meters positive above the plane normal to the ellipsoidal height up direction at the
///      origin.
/// * Normalized [`CoordSpace::Projected`] coordinates are:
///   1. easting in meters positive east of projection origin
///   2. northing in meters positive north of projection origin
///   3. ellipsoidal height in meters positive up (directly taken from normalized geodetic
///      coordinates)
///
/// [`CoordSpace`]s are also partially sorted.
/// Two [`CoordSpace`] are comparable if there is a **coordinates conversion path** between each other.
/// Which one is less than the other has to do with **transformation path**: coordinates are
/// usually converted to the *lowest* CRS kind comparable to both source and destination CRS, then a
/// datum transformation is performed if necessary.
#[derive(Debug, Clone, Copy, PartialEq)]
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

pub type LoweringTransformation = (
    Box<dyn Crs>,
    Box<dyn Transformation>,
    Box<dyn Transformation>,
);

/// [`Crs`] defines common attributes of CRS and methods to find **coordinates transformation
/// pathes** between them.
/// For the latter, CRS are partially ordered based on their kind: 2 CRS are comparable if and only
/// if there is a **coordinates conversion path** between them (provided they use the same
/// [`GeodeticDatum`]). To find a **coordinates transformation path** between 2 CRS, we look for
/// a CRS that is less or equal than both source and target CRS and for which there exist a datum
/// transformation between the source and target datum or to and from a common reference datum.
pub trait Crs {
    /// Return the identifier of this CRS.
    fn id(&self) -> &Id;

    /// Return the coordinates dimension of this CRS.
    fn dim(&self) -> usize;

    /// Return the kind of this CRS.
    fn kind(&self) -> CoordSpace;

    /// Return the datum used by this CRS.
    fn datum(&self) -> &GeodeticDatum;

    // TODO: Add a method to tell if the CRS is normalized for its [`kind`][Self::kind()].
    // fn is_normalized(&self) -> bool

    /// Return the **normalized CRS** derived from this one with a lower [`CoordSpace`] if any
    /// and the conversions to convert coordinates from this CRS to and from the lowered CRS.
    fn lowered(&self) -> Option<LoweringTransformation>;
}

pub mod geocentric;
pub mod geodetic;
pub mod topocentric;
