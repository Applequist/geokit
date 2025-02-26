/// [Crs] is the root trait for Coordinates reference systems.
pub trait Crs {
    fn id(&self) -> &str;
}

pub mod geocentric;
pub mod geographic;
pub mod projected;
