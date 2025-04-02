use super::Crs;
use crate::cs::Coord;
use crate::cs::cartesian::projected::ProjectedAxes;
use crate::quantities::length::Length;
use crate::{geodesy::GeodeticDatum, projections::ProjectionSpec};
use smol_str::SmolStr;

/// A [ProjectedCrs] is a **2D/3D cartesian coordinates reference system** derived from a
/// 2D or 3D [Geographic] Crs using a projection.
/// Coordinates are made up of easting, northing and optionally ellipsoidal height in various order,
/// direction and quantity.
#[derive(Debug, Clone, PartialEq)]
pub struct ProjectedCrs {
    pub id: SmolStr,
    pub datum: GeodeticDatum,
    pub axes: ProjectedAxes,
    pub projection: ProjectionSpec,
}

impl Crs for ProjectedCrs {
    fn id(&self) -> &str {
        &self.id
    }

    fn dim(&self) -> usize {
        self.axes.dim()
    }

    fn dist(&self, a: &Coord, b: &Coord) -> Result<Length, &'static str> {
        assert_eq!(a.len(), self.dim());
        assert_eq!(b.len(), self.dim());
        let na = self.axes.normalize(a);
        let nb = self.axes.normalize(b);
        Ok(na.dist_to(nb))
    }
}

#[cfg(test)]
mod tests {
    use super::ProjectedCrs;
    use crate::{
        crs::Crs,
        cs::cartesian::projected::{ProjectedAxes, ProjectedTolerance},
        geodesy::geodetic_datum::consts::WGS84,
        projections::ProjectionSpec,
        units::length::M,
    };
    use approx::assert_abs_diff_eq;
    use std::any::{Any, TypeId};

    #[test]
    fn type_id() {
        let crs: &dyn Any = &ProjectedCrs {
            id: "UTM zone1".into(),
            datum: WGS84,
            axes: ProjectedAxes::EastNorthUp {
                horiz_unit: M,
                height_unit: M,
            },
            projection: ProjectionSpec::UTMNorth { zone: 1 },
        };

        assert_eq!(crs.type_id(), TypeId::of::<ProjectedCrs>());
        assert!(crs.is::<ProjectedCrs>());
    }

    #[test]
    fn approx_eq() {
        let crs = ProjectedCrs {
            id: "UTM Zon1".into(),
            datum: WGS84,
            axes: ProjectedAxes::EastNorth { horiz_unit: M },
            projection: ProjectionSpec::UTMNorth { zone: 1 },
        };

        let tol = crs.axes.denormalize_tol(ProjectedTolerance::small());
        assert!(crs.approx_eq(&[1., 1. + 1e-3], &[1., 1.], &tol));
        assert!(!crs.approx_eq(&[1., 1. + 1e-1], &[1., 1.], &tol));
    }

    #[test]
    fn dist() {
        let crs = ProjectedCrs {
            id: "UTM 1".into(),
            datum: WGS84,
            axes: ProjectedAxes::EastNorth { horiz_unit: M },
            projection: ProjectionSpec::UTMNorth { zone: 1 },
        };

        let d = crs.dist(&[0., 0.], &[1., 1.]);
        assert!(d.is_ok());
        assert_abs_diff_eq!(d.unwrap(), 2.0f64.sqrt() * M);
    }
}
