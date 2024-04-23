/// Define a scaling factor in **part per million**
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct PPM(pub f64);

impl PPM {
    /// Return the *unitless* scaling factor
    #[inline]
    pub fn factor(&self) -> f64 {
        1.0 + self.0 * 1e-6
    }

    #[inline]
    pub fn inv_factor(&self) -> f64 {
        1.0 - self.0 * 1e-6
    }
}