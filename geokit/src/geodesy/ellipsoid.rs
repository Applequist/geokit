/// An `Ellipsoid` is a mathematical surface used as a model of the Earth surface.
#[derive(Debug, Copy, Clone, PartialEq)]
pub struct Ellipsoid {
    /// The semi major axis length **in meters**.
    a: f64,
    /// The semi minor axis length **in meters**.
    b: f64,
    /// The inverse flattening: `a / (a - b)``. INFINITY if a == b.
    invf: f64,
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
            a,
            b,
            invf: if a > b { a / (a - b) } else { f64::INFINITY },
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
            a,
            b: a * (1. - 1. / invf),
            invf,
        }
    }

    /// Return whether the ellipsoid is actually a sphere or not.
    pub fn is_spherical(&self) -> bool {
        self.invf.is_infinite()
    }

    /// Return the semi-major axis length in meters.
    #[inline]
    pub fn a(&self) -> f64 {
        self.a
    }

    /// Return the semi-major axis squared.
    #[inline]
    pub fn a_sq(&self) -> f64 {
        self.a * self.a
    }

    /// Return the semi-minor axis length in meters.
    #[inline]
    pub fn b(&self) -> f64 {
        self.b
    }

    /// Return the semi minor axis squared.
    #[inline]
    pub fn b_sq(&self) -> f64 {
        self.b * self.b
    }

    /// Return the inverse flattening: `a / (a - b)`.
    #[inline]
    pub fn invf(&self) -> f64 {
        self.invf
    }

    /// Return the first eccentricity squared: `(a^2 -b^2) / a^2`
    #[inline]
    pub fn e_sq(&self) -> f64 {
        (self.a_sq() - self.b_sq()) / self.a_sq()
    }

    /// Return the radius of curvature **in meters** in the east-west direction
    /// at the given latitude **in radians**.
    pub fn prime_vertical_radius(&self, lat: f64) -> f64 {
        self.a / (1.0 - self.e_sq() * lat.sin().powi(2)).sqrt()
    }

    /// Return the radius of curvature **in meters** in the north-south direction
    /// at the given latitude **in radians**.
    pub fn prime_meridional_radius(&self, lat: f64) -> f64 {
        self.a * (1.0 - self.e_sq()) / (1.0 - self.e_sq() * lat.sin().powi(2)).powf(1.5)
    }

    /// Convert **normalized geodetic coordinates** (lon in rad, lat in rad, height in meters)
    /// into **normalized geocentric coordinates** (x, y, z) all in meters.
    pub fn llh_to_xyz(&self, llh: &[f64], xyz: &mut [f64]) {
        let lon = llh[0];
        let lat = llh[1];
        let h = llh[2];

        let v = self.prime_vertical_radius(lat);
        let (sin_lon, cos_lon) = lon.sin_cos();
        let (sin_lat, cos_lat) = lat.sin_cos();

        xyz[0] = (v + h) * cos_lat * cos_lon;
        xyz[1] = (v + h) * cos_lat * sin_lon;
        xyz[2] = (v * (1.0 - self.e_sq()) + h) * sin_lat;
    }

    /// Convert **normalized geocentric coordinates** (x, y, z) in meters
    /// into **normalized geodetic coordinates** (lon in rad, lat in rad, height in meters)
    /// using Heiskanen and Moritz iterative method.
    pub fn xyz_to_llh(&self, xyz: &[f64], llh: &mut [f64]) {
        let x = xyz[0];
        let y = xyz[1];
        let z = xyz[2];

        let a2 = self.a_sq();
        let b2 = self.b_sq();
        let e2 = self.e_sq();

        let lon = y.atan2(x);

        let p = x.hypot(y);
        let mut lat = z.atan2(p * (1.0 - e2));
        let (sin_lat, cos_lat) = lat.sin_cos();
        let n = a2 / (a2 * cos_lat * cos_lat + b2 * sin_lat * sin_lat).sqrt();
        let mut h = p / cos_lat - n;
        loop {
            let next_lat = z.atan2(p * (1.0 - e2 * n / (n + h)));
            let (sin_nlat, cos_nlat) = next_lat.sin_cos();
            let next_n = a2 / ((a2 * cos_nlat * cos_nlat) + b2 * sin_nlat * sin_nlat).sqrt();
            let next_h = p / cos_nlat - next_n;
            let delta_lat = (lat - next_lat).abs();
            let delta_h = (h - next_h).abs();
            lat = next_lat;
            h = next_h;
            if delta_lat < 0.5e-5 && delta_h < 0.5e-3 {
                break;
            }
        }

        llh[0] = lon;
        llh[1] = lat;
        llh[2] = h;
    }
}

impl Default for Ellipsoid {
    /// Return the WGS84 (EPSG:7030) ellipsoid as default.
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
        let _a = e.a;
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
        assert_eq!(wgs84.a, 6_378_137.0);
        assert_eq!(wgs84.invf, 298.257_223_563);
    }
}
