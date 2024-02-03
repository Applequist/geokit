/// A `PrimeMeridian` defines the origin of longitudes.
/// It is defined by its longitude with respect to the Greenwich meridian,
/// expressed **in radians** and positive eastward.
#[derive(Debug, Copy, Clone, PartialEq)]
pub struct PrimeMeridian {
    gw_lon_rad: f64,
}

impl PrimeMeridian {
    /// Create a new [`PrimeMeridian`] with the given Greenwich longitude **in radians, positive
    /// East of Greenwich**.
    ///
    /// # Panics
    ///
    /// Panics if the `greenwich_longitude` is not in (-pi, pi].
    pub fn with_gw_lon(gw_lon_rad: f64) -> Self {
        assert!(
            gw_lon_rad > -std::f64::consts::PI && gw_lon_rad <= std::f64::consts::PI,
            "Expected greenwich_longitude in (-pi, pi]. Got {}",
            gw_lon_rad
        );
        Self::new(gw_lon_rad)
    }

    /// Create a new [PrimeMeridian] with the given Greenwich longitude **in (-pi, pi] radians**
    /// positive east of Greenwich.
    pub(crate) const fn new(gw_lon_rad: f64) -> Self {
        Self { gw_lon_rad }
    }

    /// Return this prime meridian's longitude in radians, positive east of the Greenwich
    /// prime meridian.
    #[inline]
    pub fn lon(&self) -> f64 {
        self.gw_lon_rad
    }

    /// Convert a **normalized longitude** with the prime meridian as origin into
    /// a **normalized longitude** with the Greenwich prime meridian as origin.
    #[inline]
    pub(crate) fn convert_lon_to_gw(&self, lon: f64) -> f64 {
        // FIX: What if we cross the antimeridian?
        lon + self.gw_lon_rad
    }

    /// Convert a **normalized longitude** with the Greenwich prime meridian as origin into
    /// a **normalized longitude** with this prime meridian as origin.
    #[inline]
    pub(crate) fn convert_lon_from_gw(&self, lon: f64) -> f64 {
        // FIX: What if we cross the antimeridian?
        lon - self.gw_lon_rad
    }
}

impl Default for PrimeMeridian {
    /// Return the Greenwich prime meridian as default.
    fn default() -> Self {
        PrimeMeridian { gw_lon_rad: 0.0 }
    }
}

#[cfg(test)]
mod tests {
    use super::PrimeMeridian;

    #[test]
    #[should_panic = "Expected greenwich_longitude in (-pi, pi]"]
    fn longitude_lt_mpi() {
        let _p = PrimeMeridian::with_gw_lon(-std::f64::consts::PI - f64::EPSILON);
    }

    #[test]
    #[should_panic = "Expected greenwich_longitude in (-pi, pi]"]
    fn longitude_eq_mpi() {
        let _p = PrimeMeridian::with_gw_lon(-std::f64::consts::PI);
    }

    #[test]
    #[should_panic = "Expected greenwich_longitude in (-pi, pi]"]
    fn longitude_gt_pi() {
        // FIX: shoud panic when only add 1.0 * f64::EPSILON!
        let _p = PrimeMeridian::with_gw_lon(std::f64::consts::PI + 2.0 * f64::EPSILON);
    }

    #[test]
    fn copy() {
        let pm = PrimeMeridian::new(2.23_f64.to_radians());
        let cpy = pm;
        assert_eq!(pm, cpy);
        let _s = cpy.gw_lon_rad;
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
        assert_eq!(pm.gw_lon_rad, 0.0);
    }
}
