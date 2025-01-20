use smol_str::SmolStr;

use crate::cs::geodetic::Lon;

/// A `PrimeMeridian` defines the origin of longitudes.
/// It is defined by its longitude with respect to the Greenwich meridian,
/// expressed **in radians** and positive eastward.
#[derive(Debug, Clone)]
pub struct PrimeMeridian {
    name: SmolStr,
    lon: f64,
}

impl PrimeMeridian {
    /// Create a new [`PrimeMeridian`] with the given Greenwich longitude **positive
    /// East of Greenwich**.
    pub fn new(name: &str, gw_lon: Lon) -> Self {
        Self {
            name: SmolStr::new(name),
            lon: gw_lon.normalize().rad(),
        }
    }

    /// Create a new [PrimeMeridian] with the given Greenwich longitude **in (-pi, pi] radians**
    /// positive east of Greenwich.
    #[inline(always)]
    pub(crate) const fn new_static(name: &'static str, gw_lon_rad: f64) -> Self {
        Self {
            name: SmolStr::new_static(name),
            lon: gw_lon_rad,
        }
    }

    pub fn name(&self) -> &str {
        self.name.as_str()
    }

    /// Return this prime meridian's Greenwich longitude **in (-pi..pi] radians, positive east of the Greenwich
    /// prime meridian**.
    #[inline]
    pub fn lon(&self) -> f64 {
        self.lon
    }
}

impl PartialEq for PrimeMeridian {
    fn eq(&self, other: &Self) -> bool {
        self.lon == other.lon
    }
}

/// Well known prime meridians.
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
    use std::f64::consts::PI;

    use crate::{
        cs::geodetic::Lon,
        quantity::angle::units::{DEG, RAD},
    };

    use super::PrimeMeridian;

    #[test]
    fn longitude_eq_mpi() {
        let _p = PrimeMeridian::new("EQ lower bound", Lon::new(-PI * RAD));
        assert_eq!(_p.lon(), PI);
    }

    #[test]
    fn clone() {
        let pm = PrimeMeridian::new("Paris", Lon::new(2.23 * DEG));
        let cpy = pm.clone();
        assert_eq!(pm, cpy);
        let _s = cpy.lon;
    }

    #[test]
    fn parital_eq() {
        let pm = PrimeMeridian::new("Paris", Lon::new(2.23 * DEG));
        let cpy = PrimeMeridian::new("PM", Lon::new(2.23 * DEG));

        assert!(pm.eq(&cpy));
        assert!(!pm.ne(&cpy));

        let pm2 = PrimeMeridian::new("Greenwich", Lon::new(0.0 * DEG));
        assert_ne!(pm, pm2);
    }
}
