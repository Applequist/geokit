use smol_str::SmolStr;

/// A `PrimeMeridian` defines the origin of longitudes.
/// It is defined by its longitude with respect to the Greenwich meridian,
/// expressed **in radians** and positive eastward.
#[derive(Debug, Clone)]
pub struct PrimeMeridian {
    name: SmolStr,
    gw_lon_rad: f64,
}

impl PrimeMeridian {
    /// Create a new [`PrimeMeridian`] with the given Greenwich longitude **in radians, positive
    /// East of Greenwich**.
    ///
    /// # Panics
    ///
    /// Panics if the `greenwich_longitude` is not in (-pi, pi].
    pub fn new(name: &str, gw_lon_rad: f64) -> Self {
        assert!(
            gw_lon_rad > -std::f64::consts::PI && gw_lon_rad <= std::f64::consts::PI,
            "Expected greenwich_longitude in (-pi, pi]. Got {}",
            gw_lon_rad
        );
        Self {
            name: SmolStr::new(name),
            gw_lon_rad,
        }
    }

    /// Create a new [PrimeMeridian] with the given Greenwich longitude **in (-pi, pi] radians**
    /// positive east of Greenwich.
    ///
    /// Panics if `name.len() > 23` or `gw_lon_rad <= -PI` or `gw_lon_rad > PI`
    #[inline(always)]
    pub(crate) const fn new_static(name: &'static str, gw_lon_rad: f64) -> Self {
        Self {
            name: SmolStr::new_static(name),
            gw_lon_rad,
        }
    }

    pub fn name(&self) -> &str {
        self.name.as_str()
    }

    /// Return this prime meridian's longitude in radians, positive east of the Greenwich
    /// prime meridian.
    #[inline]
    pub fn lon(&self) -> f64 {
        self.gw_lon_rad
    }
}

impl PartialEq for PrimeMeridian {
    fn eq(&self, other: &Self) -> bool {
        self.gw_lon_rad == other.gw_lon_rad
    }
}

pub mod consts {

    use super::PrimeMeridian;

    macro_rules! deg {
        ($d:expr) => {
            $d * std::f64::consts::PI / 180.0
        };
    }

    macro_rules! prime_meridians {
        ( $( $name:ident => ($id:literal, lon = $lon:expr) ),+ ) => {

            $(pub const $name: PrimeMeridian = PrimeMeridian::new_static($id, $lon);)+
        }
    }

    prime_meridians! {
        GREENWICH => ("Greenwich", lon = 0.),
        LISBON => ("Lisbon", lon = deg!(-9.131906111111)),
        PARIS => ("Paris", lon = deg!(2.337229166667)),
        BOGOTA => ("Bogota", lon = deg!(-74.080916666667)),
        MADRID => ("Madrid", lon = deg!(-3.687938888889)),
        ROME => ("Rome", lon = deg!(12.452333333333)),
        BERN => ("Bern", lon = deg!(7.439583333333)),
        JAKARTA => ("Jakarta", lon = deg!(106.807719444444)),
        FERRO => ("Ferro", lon = deg!(-17.666666666667)),
        BRUSSELS => ("Brussels", lon = deg!(4.367975)),
        STOCKHOLM => ("Stockholm", lon = deg!(18.058277777778)),
        ATHENS => ("Athens", lon = deg!(23.7163375)),
        OSLO => ("Oslo", lon = deg!(10.722916666667)),
        COPENHAGEN => ("Copenhagen", lon = deg!(12.57788))
    }
}

#[cfg(test)]
mod tests {
    use super::PrimeMeridian;

    #[test]
    #[should_panic = "Expected greenwich_longitude in (-pi, pi]"]
    fn longitude_lt_mpi() {
        let _p = PrimeMeridian::new("LT lower bound", -std::f64::consts::PI - f64::EPSILON);
    }

    #[test]
    #[should_panic = "Expected greenwich_longitude in (-pi, pi]"]
    fn longitude_eq_mpi() {
        let _p = PrimeMeridian::new("EQ lower bound", -std::f64::consts::PI);
    }

    #[test]
    #[should_panic = "Expected greenwich_longitude in (-pi, pi]"]
    fn longitude_gt_pi() {
        // FIX: shoud panic when only add 1.0 * f64::EPSILON!
        let _p = PrimeMeridian::new("GT upper bound", std::f64::consts::PI + 2.0 * f64::EPSILON);
    }

    #[test]
    fn clone() {
        let pm = PrimeMeridian::new("Paris", 2.23_f64.to_radians());
        let cpy = pm.clone();
        assert_eq!(pm, cpy);
        let _s = cpy.gw_lon_rad;
    }

    #[test]
    fn parital_eq() {
        let pm = PrimeMeridian::new("Paris", 2.23_f64.to_radians());
        let cpy = PrimeMeridian::new("PM", 2.23_f64.to_radians());

        assert!(pm.eq(&cpy));
        assert!(!pm.ne(&cpy));

        let pm2 = PrimeMeridian::new("Greenwich", 0.0);
        assert_ne!(pm, pm2);
    }
}
