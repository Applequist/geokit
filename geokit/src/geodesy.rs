use derive_more::Display;

/// The parameters defining an [Ellipsoid].
#[derive(Debug, Display)]
enum EllipsoidParams {
    #[display("a = {} m, b = {} m", a, b)]
    AB { a: f64, b: f64 },
    #[display("a = {} m, invf = {}", a, invf)]
    AInvf { a: f64, invf: f64 },
}

/// An `Ellipsoid` is a mathematical surface used as a model of the Earth surface.
#[derive(Debug, Display)]
#[display("Ellipsoid: ({})", params)]
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
        assert!(b > 0, "Expected b > 0. Got {}", b);
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
        assert!(a > 0, "Expected a > 0. Got {}", a);
        assert!(invf > 1, "Expected invf > 1. Got {}", invf);
        Self {
            params: EllipsoidParams::AInvf { a, invf },
        }
    }
}

/// A `PrimeMeridian` defines the origin of longitudes.
/// It is defined by its longitude with respect to the Greenwich meridian,
/// expressed **in radians** and positive eastward.
#[derive(Debug, Display)]
#[display("PrimeMeridian: (greenwich_longitude = {} rad)", greenwich_longitude)]
pub struct PrimeMeridian {
    greenwich_longitude: f64,
}

impl PrimeMeridian {
    /// Create a new prime meridian at the given greenwich longitude **normalized in (-pi, pi].
    pub fn new(greenwich_longitude: f64) -> Self {
        assert!(
            greenwich_longitude > std::f64::consts::PI
                && greenwich_longitude <= std::f64::consts::PI,
            "Expected greenwich_longitude in (-pi, pi]. Got {}",
            greenwich_longitude
        );
        Self {
            greenwich_longitude,
        }
    }
}

/// A `datum` is the information required to fix a coordinate system to an object.
/// A `GeodeticDatum` is a `datum` describing the relationship of an ellipsoidal model of the Earth
/// with the real Earth.
/// It is defined by an [Ellipsoid] and a [PrimeMeridian].
#[derive(Debug, Display)]
#[display(
    "GeodeticDatum: (ellipsoid = {}, prime_meridian = {})",
    ellipsoid,
    prime_meridian
)]
pub struct GeodeticDatum {
    pub ellipsoid: Ellipsoid,
    pub prime_meridian: PrimeMeridian,
}
