use super::Crs;
use crate::{
    cs::{Coord, geodetic::GeodeticAxes},
    geodesy::{
        GeodeticDatum,
        geodesics::{GeodesicSolver, vincenty::VincentyGeodesicSolver},
    },
    quantities::length::Length,
};
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
    fn id(&self) -> &str {
        &self.id
    }

    fn dim(&self) -> usize {
        self.axes.dim()
    }

    /// Computes a distance between the 2 given coordinates.
    ///
    /// When `self.dim() == 2`, this method computes the **geodesic distance** between the 2 points.
    /// when `self.dim() == 3`, this method computes the **cartesian distance** between the 2
    /// points.
    fn dist(&self, a: &Coord, b: &Coord) -> Result<Length, &'static str> {
        assert_eq!(a.len(), self.dim());
        assert_eq!(b.len(), self.dim());
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

#[cfg(test)]
mod tests {
    use super::GeographicCrs;
    use crate::{
        crs::Crs,
        cs::geodetic::{GeodeticAxes, GeodeticTolerance},
        geodesy::geodetic_datum::consts::WGS84,
        units::{angle::DEG, length::M},
    };
    use approx::assert_abs_diff_eq;
    use std::any::{Any, TypeId};

    #[test]
    fn type_id() {
        let crs: &dyn Any = &GeographicCrs {
            id: "WGS84".into(),
            datum: WGS84,
            axes: GeodeticAxes::EastNorthUp {
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
            datum: WGS84,
            axes: GeodeticAxes::EastNorthUp {
                angle_unit: DEG,
                height_unit: M,
            },
        };

        let tol = crs.axes.denormalize_tol(GeodeticTolerance::tiny());

        assert!(crs.approx_eq(&[0., 0., 0.], &[0., 1e-13f64.to_degrees(), 0.], &tol));
        assert!(!crs.approx_eq(&[0., 0., 0.], &[0., 1e-10f64.to_degrees(), 0.], &tol));
        assert!(!crs.approx_eq(&[0., 0., 0.], &[0., 0., 1e-3], &tol));
    }

    #[test]
    fn dist() {
        // Compute linear cartesian distance
        let crs3d = GeographicCrs {
            id: "WGS84".into(),
            datum: WGS84,
            axes: GeodeticAxes::EastNorthUp {
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
            datum: WGS84,
            axes: GeodeticAxes::EastNorth { angle_unit: DEG },
        };

        let d = crs2d.dist(&[0., 0.], &[0.01, 0.01]);
        assert!(d.is_ok());
        assert_abs_diff_eq!(d.unwrap(), 1569.0347 * M);
    }
}
