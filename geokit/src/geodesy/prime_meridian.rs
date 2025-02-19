use crate::cs::geodetic::Lon;
use derive_more::derive::Display;
use smol_str::SmolStr;

/// A [PrimeMeridian] defines the origin of longitudes.
/// It is represented by its longitude with respect to the Greenwich meridian,
/// positive eastward.
#[derive(Debug, Clone, Display)]
#[display("(name = {}, lon = {})", name, lon.angle().dms())]
pub struct PrimeMeridian {
    name: SmolStr,
    lon: Lon,
}

impl PrimeMeridian {
    /// Create a new [`PrimeMeridian`] with the given Greenwich longitude **positive
    /// East of Greenwich**.
    pub fn new(name: &str, gw_lon: Lon) -> Self {
        Self {
            name: SmolStr::new(name),
            lon: gw_lon.normalized(),
        }
    }

    /// Create a new [PrimeMeridian] with the given Greenwich longitude **in (-pi, pi] radians**
    /// positive east of Greenwich.
    #[inline(always)]
    pub(crate) const fn new_static(name: &'static str, gw_lon: Lon) -> Self {
        Self {
            name: SmolStr::new_static(name),
            lon: gw_lon,
        }
    }

    pub fn name(&self) -> &str {
        self.name.as_str()
    }

    /// Return this prime meridian's Greenwich longitude **in (-pi..pi] radians, positive east of the Greenwich
    /// prime meridian**.
    #[inline]
    pub fn lon(&self) -> Lon {
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
    use crate::cs::geodetic::Lon;
    use crate::quantities::angle::Angle;
    use crate::units::angle::DEG;

    macro_rules! prime_meridians {
        ( $( $name:ident => ($id:literal, lon = $lon:expr) ),+ ) => {

            $(
                pub const $name: PrimeMeridian = PrimeMeridian::new_static(
                    $id,
                    Lon::const_new(Angle::new($lon, DEG))
                );
            )+
        }
    }

    prime_meridians! {
        GREENWICH => ("Greenwich", lon = 0.),
        LISBON => ("Lisbon", lon = -9.131906111111),
        PARIS => ("Paris", lon = 2.337229166667),
        BOGOTA => ("Bogota", lon = -74.080916666667),
        MADRID => ("Madrid", lon = -3.687938888889),
        ROME => ("Rome", lon = 12.452333333333),
        BERN => ("Bern", lon = 7.439583333333),
        JAKARTA => ("Jakarta", lon = 106.807719444444),
        FERRO => ("Ferro", lon = -17.666666666667),
        BRUSSELS => ("Brussels", lon = 4.367975),
        STOCKHOLM => ("Stockholm", lon = 18.058277777778),
        ATHENS => ("Athens", lon = 23.7163375),
        OSLO => ("Oslo", lon = 10.722916666667),
        COPENHAGEN => ("Copenhagen", lon = 12.57788)
    }
}

#[cfg(test)]
mod tests {
    use std::f64::consts::PI;

    use crate::{
        cs::geodetic::Lon,
        geodesy::prime_meridian::consts::PARIS,
        units::angle::{DEG, RAD},
    };

    use super::PrimeMeridian;

    #[test]
    fn longitude_eq_mpi() {
        let _p = PrimeMeridian::new("EQ lower bound", Lon::new(-PI * RAD));
        assert_eq!(_p.lon(), Lon::MAX);
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

    #[test]
    fn pm_display() {
        assert_eq!(
            format!("{}", PARIS),
            "(name = Paris, lon =    2° 20′ 14.02500000″)"
        );
    }
}
