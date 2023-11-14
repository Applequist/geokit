use crate::coord::CoordSpace;
use crate::geodesy::GeodeticDatum;
use crate::id::Id;

pub enum LowerTransformation {
    TO,
    FROM,
}

/// [`Crs`] defines common attributes of CRS and methods to find **coordinates transformation
/// pathes** between them.
/// For the latter, CRS are partially ordered based on their kind: 2 CRS are comparable if and only
/// if there is a **coordinates conversion path** between them (provided they use the same
/// [`GeodeticDatum`]). To find a **coordinates transformation path** between 2 CRS, we look for
/// a CRS that is less or equal than both source and target CRS and for which there exist a datum
/// transformation between the source and target datum or to and from a common reference datum.
pub trait Crs {
    /// Returns unique identifier of this CRS.
    /// The returned `id` has the following format:
    ///
    /// ```ignore
    /// <id> := [<authority> ':'] <code>
    /// <authority> := String
    /// <code> := String
    /// ```
    fn id(&self) -> &Id;

    /// Returns the coordinates dimension of this CRS.
    fn dim(&self) -> usize;

    /// Returns the kind of this CRS.
    fn kind(&self) -> CoordSpace;

    /// Returns the datum used by this CRS.
    fn datum(&self) -> &GeodeticDatum;

    /// Returns the lower CRS derived from this one if any and the tranformation to or from it.
    ///
    /// # Arguments
    ///
    /// * `transformation`: determine the *direction* of the returned transformation:
    ///   * [`LowerTransformation::TO`] returns the to-lower tranformation,
    ///   * [`LowerTransformation::FROM`] returns the from-lower transformation.
    /// FIX: Replace String in the return type by a Box<dyn Transformation> once the transformation path finding algorithm works
    fn lower(&self, transformation: LowerTransformation) -> Option<(Box<dyn Crs>, String)>;
}

pub mod geocentric;
pub mod geodetic;
pub mod topocentric;
