/// A `PrimeMeridian` defines the origin of longitudes.
/// It is defined by its longitude with respect to the Greenwich meridian,
/// expressed **in radians** and positive eastward.
#[derive(Debug, Copy, Clone, PartialEq)]
pub struct PrimeMeridian {
    greenwich_longitude: f64,
}

impl PrimeMeridian {
    /// Creates a new [`PrimeMeridian`] with the given Greenwich longitude **in radians, positive
    /// East of Greenwich**.
    ///
    /// # Panics
    ///
    /// Panics if the `greenwich_longitude` is not in (-pi, pi].
    pub fn new(greenwich_longitude: f64) -> Self {
        assert!(
            greenwich_longitude > -std::f64::consts::PI
                && greenwich_longitude <= std::f64::consts::PI,
            "Expected greenwich_longitude in (-pi, pi]. Got {}",
            greenwich_longitude
        );
        Self {
            greenwich_longitude,
        }
    }
}

impl Default for PrimeMeridian {
    /// Returns the Greenwich prime meridian as default.
    fn default() -> Self {
        PrimeMeridian {
            greenwich_longitude: 0.0,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::PrimeMeridian;

    #[test]
    #[should_panic = "Expected greenwich_longitude in (-pi, pi]"]
    fn longitude_lt_mpi() {
        let _p = PrimeMeridian::new(-std::f64::consts::PI - f64::EPSILON);
    }

    #[test]
    #[should_panic = "Expected greenwich_longitude in (-pi, pi]"]
    fn longitude_eq_mpi() {
        let _p = PrimeMeridian::new(-std::f64::consts::PI);
    }

    #[test]
    #[should_panic = "Expected greenwich_longitude in (-pi, pi]"]
    fn longitude_gt_pi() {
        // FIX: shoud panic when only add 1.0 * f64::EPSILON!
        let _p = PrimeMeridian::new(std::f64::consts::PI + 2.0 * f64::EPSILON);
    }

    #[test]
    fn copy() {
        let pm = PrimeMeridian::new(2.23_f64.to_radians());
        let cpy = pm;
        assert_eq!(pm, cpy);
        let _s = cpy.greenwich_longitude;
    }

    #[test]
    fn parital_eq() {
        let pm = PrimeMeridian::new(2.23_f64.to_radians());
        let cpy = pm;

        assert!(pm.eq(&cpy));
        assert!(!pm.ne(&cpy));

        let pm2 = PrimeMeridian::new(0.0);
        assert_ne!(pm, pm2);
    }

    #[test]
    fn default() {
        let pm = PrimeMeridian::default();
        assert_eq!(pm.greenwich_longitude, 0.0);
    }
}
