use super::Crs;
use crate::cs::cartesian::projected::{ProjectedAxes, ProjectedTolerance};
use crate::quantities::length::Length;
use crate::units::length::M;
use crate::{geodesy::GeodeticDatum, math::fp::Float, projections::ProjectionSpec};
use approx::AbsDiffEq;
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

    fn dist(&self, a: &[Float], b: &[Float]) -> Result<Length, &'static str> {
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
        units::length::M,
    };
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
            projection: crate::projections::ProjectionSpec::UTMNorth { zone: 1 },
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
            projection: crate::projections::ProjectionSpec::UTMNorth { zone: 1 },
        };

        let tol = crs.axes.denormalize_tol(ProjectedTolerance::small());
        assert!(crs.approx_eq(&[1., 1. + 1e-3], &[1., 1.], &tol));
        assert!(!crs.approx_eq(&[1., 1. + 1e-1], &[1., 1.], &tol));
    }

    #[test]
    fn dist() {
        todo!()
    }
}
