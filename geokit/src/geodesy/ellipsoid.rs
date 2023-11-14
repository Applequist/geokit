/// The parameters defining an [Ellipsoid].
#[derive(Debug, Copy, Clone, PartialEq)]
enum EllipsoidParams {
    AB { a: f64, b: f64 },
    AInvf { a: f64, invf: f64 },
}

/// An `Ellipsoid` is a mathematical surface used as a model of the Earth surface.
#[derive(Debug, Copy, Clone, PartialEq)]
pub struct Ellipsoid {
    params: EllipsoidParams,
}

impl Ellipsoid {
    /// Create a new ellipsoid using semi-major and semi-minor axes **in meters**.
    ///
    /// # Panics
    ///
    /// If `b` is negative or zero or if `a` is less then `b`.
    pub fn from_ab(a: f64, b: f64) -> Self {
        assert!(b > 0., "Expected b > 0. Got {}", b);
        assert!(a >= b, "Expected a >= b. Got a = {}, b = {}", a, b);
        Self {
            params: EllipsoidParams::AB { a, b },
        }
    }

    /// Create a new ellipsoid using semi-major axis **in meters** and inverse flattening.
    ///
    /// # Panics
    ///
    /// if `a` is negative or zero or if `invf` is not greater than 1.
    pub fn from_ainvf(a: f64, invf: f64) -> Self {
        assert!(a > 0., "Expected a > 0. Got {}", a);
        assert!(invf > 1., "Expected invf > 1. Got {}", invf);
        Self {
            params: EllipsoidParams::AInvf { a, invf },
        }
    }
}

impl Default for Ellipsoid {
    /// Returns the WGS84 (epsg:7030) ellipsoid as default.
    fn default() -> Self {
        Ellipsoid::from_ainvf(6_378_137.0, 298.257_223_563)
    }
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    #[should_panic = "Expected b > 0"]
    fn negative_b() {
        let _ellipsoid = Ellipsoid::from_ab(1.0, -0.5);
    }

    #[test]
    #[should_panic = "Expected b > 0"]
    fn zero_b() {
        let _e = Ellipsoid::from_ab(1.0, 0.0);
    }

    #[test]
    #[should_panic = "Expected a >= b"]
    fn a_less_then_b() {
        let _ellipsoid = Ellipsoid::from_ab(1.0, 2.0);
    }

    #[test]
    #[should_panic = "Expected a > 0"]
    fn negative_a() {
        let _ellipsoid = Ellipsoid::from_ainvf(-0.1, 297.0);
    }

    #[test]
    #[should_panic = "Expected invf > 1"]
    fn small_invf() {
        let _ellipsoid = Ellipsoid::from_ainvf(1., 0.5);
    }

    #[test]
    #[should_panic = "Expected invf > 1"]
    fn unit_invf() {
        let _ellipsoid = Ellipsoid::from_ainvf(1., 1.);
    }

    #[test]
    fn copy() {
        let e = Ellipsoid::from_ab(1., 0.9);
        let cpy = e;
        assert_eq!(e, cpy);
        let _s = e.params;
    }

    #[test]
    fn partial_eq() {
        let e = Ellipsoid::from_ab(1., 0.9);
        let eq = e;
        assert!(e.eq(&eq));
        assert!(!e.ne(&eq));

        let e2 = Ellipsoid::from_ab(1.01, 0.9);
        assert!(!e.eq(&e2));
        assert!(e.ne(&e2));

        let e3 = Ellipsoid::from_ainvf(1., 10.0);
        assert!(!e.eq(&e3));
        assert!(e.ne(&e3));
    }

    #[test]
    fn default() {
        let wgs84 = Ellipsoid::default();
        assert_eq!(
            wgs84.params,
            EllipsoidParams::AInvf {
                a: 6_378_137.0,
                invf: 298.257_223_563
            }
        );
    }
}
