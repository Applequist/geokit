use super::Crs;
use crate::{
    cs::{
        cartesian::geocentric::XYZ,
        geodetic::{GeodeticAxes, GeodeticTolerance},
    },
    geodesy::{
        geodesics::{vincenty::VincentyGeodesicSolver, GeodesicSolver},
        GeodeticDatum,
    },
    math::fp::Float,
    quantities::length::Length,
    transformations::{ToXYZTransformation, ToXYZTransformationProvider, TransformationError},
    units::length::M,
};
use approx::AbsDiffEq;
use smol_str::SmolStr;

/// A [GeographicCrs] is a **2D or 3D geodetic coordinates reference system** in which
/// coordinates are made up of longitude, latitude and optionally ellipsoidal height in various order, direction
/// and quantity.
#[derive(Debug, Clone, PartialEq)]
pub struct GeographicCrs {
    pub id: SmolStr,
    pub datum: GeodeticDatum,
    pub axes: GeodeticAxes,
}

impl Crs for GeographicCrs {
    type Tolerance = GeodeticTolerance;

    fn id(&self) -> &str {
        &self.id
    }

    fn dim(&self) -> usize {
        self.axes.dim()
    }

    fn approx_eq(&self, a: &[Float], b: &[Float], tol: Self::Tolerance) -> bool {
        // Convert the errors margin into this CRS units.
        let (angle_unit, height_unit) = match self.axes {
            GeodeticAxes::EastNorthUp {
                angle_unit,
                height_unit,
            }
            | GeodeticAxes::NorthEastUp {
                angle_unit,
                height_unit,
            } => (angle_unit, height_unit),
            GeodeticAxes::EastNorth { angle_unit }
            | GeodeticAxes::NorthEast { angle_unit }
            | GeodeticAxes::NorthWest { angle_unit } => (angle_unit, M),
        };
        let converted_tol = tol.convert_to(angle_unit, height_unit);

        // Take coordinates order into account
        match self.axes {
            GeodeticAxes::EastNorthUp {
                angle_unit: _,
                height_unit: _,
            } => {
                a[0].abs_diff_eq(&b[0], converted_tol.lon.0)
                    && a[1].abs_diff_eq(&b[1], converted_tol.lat.0)
                    && a[2].abs_diff_eq(&b[2], converted_tol.height.0)
            }
            GeodeticAxes::NorthEastUp {
                angle_unit: _,
                height_unit: _,
            } => {
                a[0].abs_diff_eq(&b[0], converted_tol.lat.0)
                    && a[1].abs_diff_eq(&b[1], converted_tol.lon.0)
                    && a[2].abs_diff_eq(&b[2], converted_tol.height.0)
            }
            GeodeticAxes::EastNorth { angle_unit: _ } => {
                a[0].abs_diff_eq(&b[0], converted_tol.lon.0)
                    && a[1].abs_diff_eq(&b[1], converted_tol.lat.0)
            }
            GeodeticAxes::NorthEast { angle_unit: _ } => {
                a[0].abs_diff_eq(&b[0], converted_tol.lat.0)
                    && a[1].abs_diff_eq(&b[1], converted_tol.lon.0)
            }
            GeodeticAxes::NorthWest { angle_unit: _ } => {
                a[0].abs_diff_eq(&b[0], converted_tol.lat.0)
                    && a[1].abs_diff_eq(&b[1], converted_tol.lon.0)
            }
        }
    }

    /// Computes a distance between the 2 given coordinates.
    ///
    /// When `self.dim() == 2`, this method computes the **geodesic distance** between the 2 points.
    /// when `self.dim() == 3`, this method computes the **cartesian distance** between the 2
    /// points.
    fn dist(&self, a: &[Float], b: &[Float]) -> Result<Length, &'static str> {
        let llha = self.axes.normalize(a);
        let llhb = self.axes.normalize(b);
        if self.dim() == 2 {
            let solver = VincentyGeodesicSolver::new(self.datum.ellipsoid());
            solver
                .solve_inverse((llha.lon, llha.lat), (llhb.lon, llhb.lat))
                .map(|g| g.s)
        } else {
            let xyza = self.datum.llh_to_xyz(llha);
            let xyzb = self.datum.llh_to_xyz(llhb);
            Ok(xyza.dist_to(xyzb))
        }
    }
}

impl ToXYZTransformationProvider for GeographicCrs {
    fn to_xyz_transformation<'a>(&self) -> Box<dyn ToXYZTransformation + 'a> {
        Box::new(self.clone())
    }
}

impl ToXYZTransformation for GeographicCrs {
    fn to_xyz(&self, coords: &[Float]) -> Result<XYZ, TransformationError> {
        let llh = self.axes.normalize(coords);
        Ok(self.datum.llh_to_xyz(llh))
    }

    fn from_xyz(&self, xyz: XYZ, coords: &mut [Float]) -> Result<(), TransformationError> {
        let llh = self.datum.xyz_to_llh(xyz);
        Ok(self.axes.denormalize(llh, coords))
    }
}

#[cfg(test)]
mod tests {
    use approx::assert_abs_diff_eq;

    use super::GeographicCrs;
    use crate::{
        crs::Crs,
        cs::geodetic::GeodeticTolerance,
        geodesy::{ellipsoid::consts::WGS84, prime_meridian::consts::GREENWICH, GeodeticDatum},
        units::{angle::DEG, length::M},
    };
    use std::any::{Any, TypeId};

    #[test]
    fn type_id() {
        let crs: &dyn Any = &GeographicCrs {
            id: "WGS84".into(),
            datum: GeodeticDatum::new("WGS84", WGS84, GREENWICH),
            axes: crate::cs::geodetic::GeodeticAxes::EastNorthUp {
                angle_unit: DEG,
                height_unit: M,
            },
        };

        assert_eq!(crs.type_id(), TypeId::of::<GeographicCrs>());
        assert!(crs.is::<GeographicCrs>());
    }

    #[test]
    fn approx_eq() {
        let crs = GeographicCrs {
            id: "WGS84".into(),
            datum: GeodeticDatum::new("WGS84", WGS84, GREENWICH),
            axes: crate::cs::geodetic::GeodeticAxes::EastNorthUp {
                angle_unit: DEG,
                height_unit: M,
            },
        };

        assert!(crs.approx_eq(
            &[0., 0., 0.],
            &[0., 1e-13f64.to_degrees(), 0.],
            GeodeticTolerance::tiny()
        ));

        assert!(!crs.approx_eq(
            &[0., 0., 0.],
            &[0., 1e-10f64.to_degrees(), 0.],
            GeodeticTolerance::tiny()
        ));
        assert!(!crs.approx_eq(&[0., 0., 0.], &[0., 0., 1e-3], GeodeticTolerance::tiny()));
    }

    #[test]
    fn dist() {
        // Compute linear cartesian distance
        let crs3d = GeographicCrs {
            id: "WGS84".into(),
            datum: GeodeticDatum::new("WGS84", WGS84, GREENWICH),
            axes: crate::cs::geodetic::GeodeticAxes::EastNorthUp {
                angle_unit: DEG,
                height_unit: M,
            },
        };

        let d = crs3d.dist(&[0., 0., 0.], &[0.01, 0.01, 10.]);
        assert!(d.is_ok());
        assert_abs_diff_eq!(d.unwrap(), 1569.0678 * M);

        // Compute geodesic distance on the ellipsoid
        let crs2d = GeographicCrs {
            id: "WGS84".into(),
            datum: GeodeticDatum::new("WGS84", WGS84, GREENWICH),
            axes: crate::cs::geodetic::GeodeticAxes::EastNorth { angle_unit: DEG },
        };

        let d = crs2d.dist(&[0., 0.], &[0.01, 0.01]);
        assert!(d.is_ok());
        assert_abs_diff_eq!(d.unwrap(), 1569.0347 * M);
    }
}
