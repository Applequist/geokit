use num::{Float, One, Zero};
use std::fmt::Debug;
use std::ops::{Add, Deref, DerefMut, Div, Index, IndexMut, Mul, Neg, Rem, Sub};

use crate::traits::Scalar;

#[derive(Debug, Clone, Copy, PartialEq)]
#[repr(transparent)]
pub struct Vector<T: Scalar, const N: usize> {
    coef: [T; N],
}

impl<T: Scalar, const N: usize> From<[T; N]> for Vector<T, N> {
    fn from(value: [T; N]) -> Self {
        Self::new(value)
    }
}

impl<T: Scalar, const N: usize> Vector<T, N> {
    pub fn new(coef: [T; N]) -> Self {
        Self { coef }
    }

    /// Create a new vector with all coefficients zeroed except the n-th one set to the given
    /// value.
    pub fn nth(n: usize, val: T) -> Self {
        let coef = std::array::from_fn(|i| if i == n { val } else { T::zero() });
        Self { coef }
    }

    /// Return the sum of all the vector coefficients.
    pub fn sum(self) -> T {
        self.coef.into_iter().fold(T::zero(), |acc, c| acc + c)
    }

    /// Return the inner product of 2 vectors.
    pub fn inner(self, other: Self) -> T {
        (self * other).sum()
    }

    /// Return the squared norm of this vector.
    pub fn norm_sq(self) -> T {
        (self * self).sum()
    }
}

impl<T: Scalar + Float, const N: usize> Vector<T, N> {
    pub fn norm(&self) -> T {
        self.norm_sq().sqrt()
    }
}

impl<T: Scalar, const N: usize> Default for Vector<T, N> {
    /// The default vector is the vector with all zero coefficients.
    fn default() -> Self {
        Self {
            coef: [T::zero(); N],
        }
    }
}

impl<T: Scalar, const N: usize> Deref for Vector<T, N> {
    type Target = [T];

    fn deref(&self) -> &Self::Target {
        &self.coef
    }
}

impl<T: Scalar, const N: usize> DerefMut for Vector<T, N> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.coef
    }
}

impl<T: Scalar, const N: usize> Index<usize> for Vector<T, N> {
    type Output = T;

    fn index(&self, index: usize) -> &Self::Output {
        &self.coef[index]
    }
}

impl<T: Scalar, const N: usize> IndexMut<usize> for Vector<T, N> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.coef[index]
    }
}

impl<T: Scalar, const N: usize> Zero for Vector<T, N> {
    fn zero() -> Self {
        Vector::new([T::zero(); N])
    }

    fn is_zero(&self) -> bool {
        self.iter().all(|c| c.is_zero())
    }
}

impl<T: Scalar, const N: usize> One for Vector<T, N> {
    fn one() -> Self {
        Vector::new([T::one(); N])
    }
}

impl<T: Scalar, const N: usize> Add for Vector<T, N> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        let mut coef = [T::zero(); N];
        for (ix, (a, b)) in self.iter().zip(rhs.iter()).enumerate() {
            coef[ix] = *a + *b;
        }
        Vector::new(coef)
    }
}

impl<T: Scalar + Neg, const N: usize> Neg for Vector<T, N>
where
    <T as Neg>::Output: Scalar,
{
    type Output = Vector<<T as Neg>::Output, N>;

    fn neg(self) -> Self::Output {
        let coef = std::array::from_fn(|i| -self[i]);
        Vector::<<T as Neg>::Output, N>::new(coef)
    }
}

impl<T: Scalar, const N: usize> Sub for Vector<T, N> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        let mut coef = [T::zero(); N];
        for (ix, (a, b)) in self.iter().zip(rhs.iter()).enumerate() {
            coef[ix] = *a - *b;
        }
        Vector::new(coef)
    }
}

impl<T: Scalar, const N: usize> Mul<T> for Vector<T, N> {
    type Output = Self;

    fn mul(self, rhs: T) -> Self::Output {
        let mut coef = [T::zero(); N];
        for ix in 0..N {
            coef[ix] = self[ix] * rhs;
        }
        Self { coef }
    }
}

impl<T: Scalar, const N: usize> Mul for Vector<T, N> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        let mut coef = [T::zero(); N];
        for (ix, (a, b)) in self.iter().zip(rhs.iter()).enumerate() {
            coef[ix] = *a * *b;
        }
        Vector::new(coef)
    }
}

impl<T: Scalar, const N: usize> Div<T> for Vector<T, N> {
    type Output = Self;

    fn div(self, rhs: T) -> Self::Output {
        let mut coef = [T::zero(); N];
        for ix in 0..N {
            coef[ix] = self[ix] / rhs;
        }
        Vector::new(coef)
    }
}

impl<T: Scalar, const N: usize> Div for Vector<T, N> {
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        let mut coef = [T::zero(); N];
        for (ix, (a, b)) in self.iter().zip(rhs.iter()).enumerate() {
            coef[ix] = *a / *b;
        }
        Vector::new(coef)
    }
}

impl<T: Scalar, const N: usize> Rem for Vector<T, N> {
    type Output = Vector<T, N>;

    fn rem(self, rhs: Self) -> Self::Output {
        let coef = std::array::from_fn(|i| self[i].rem(rhs[i]));
        Self { coef }
    }
}

pub type Vec2<T = f64> = Vector<T, 2>;
pub type Vec3<T = f64> = Vector<T, 3>;
pub type Vec4<T = f64> = Vector<T, 4>;

#[cfg(test)]
mod tests {
    use std::f64::EPSILON;

    use num::{One, Zero};

    use crate::linalg::vector::{Vec2, Vec3};

    #[test]
    fn debug() {
        assert_eq!(
            format!("{:?}", Vec3::<f64>::zero()),
            "Vector { coef: [0.0, 0.0, 0.0] }"
        );
    }

    #[test]
    fn eq() {
        let u = Vec3::new([0., 1., 2.]);
        assert_eq!(u, u);
        let v = Vec3::new([0. + EPSILON, 1., 2.]);
        assert!(!u.eq(&v));
        assert!(!v.eq(&u));
        assert!(u.ne(&v));
        assert!(v.ne(&u));
    }

    #[test]
    fn nth() {
        let i = Vec3::nth(0, 1.);
        assert_eq!(i, Vec3::new([1., 0., 0.]));
    }

    #[test]
    fn sum() {
        let one = Vec3::new([1., 2., 3.]);
        assert_eq!(one.sum(), 6.0);
    }

    #[test]
    fn inner() {
        let u = Vec3::nth(1, 1.);
        let v = Vec3::new([1., 2., 3.]);
        assert_eq!(u.inner(v), 2.);
    }

    #[test]
    fn norm_sq() {
        let v = Vec3::new([1., 2., 3.]);
        assert_eq!(v.norm_sq(), 14.);
    }

    #[test]
    fn norm() {
        let v = Vec3::new([3., 0., 4.]);
        assert_eq!(v.norm(), 5.)
    }

    #[test]
    fn default() {
        let v: Vec3<f64> = Vec3::default();
        assert_eq!(v, Vec3::zero());
    }

    #[test]
    fn deref() {
        let v = Vec3::new([1., 2., 3.]);
        assert_eq!(*v, [1., 2., 3.]);
    }

    #[test]
    fn deref_mut() {
        let mut v = Vec3::zero();
        v[0] = 1.0;
        assert_eq!(v, Vec3::nth(0, 1.0));
    }

    #[test]
    fn index() {
        let v = Vec3::<u8>::new([1, 2, 3]);
        assert_eq!(v[0], 1);
        assert_eq!(v[1], 2);
        assert_eq!(v[2], 3);
    }

    #[test]
    fn index_mut() {
        let mut v: Vec3<u8> = Vec3::default();
        v[1] = 1;
        assert_eq!(*v, [0, 1, 0]);
    }

    #[test]
    fn zero() {
        let z: Vec3<u8> = Vec3::zero();
        let v = Vec3::new([1, 2, 3]);
        assert_eq!(v + z, v);
        assert_eq!(z + v, v);
    }

    #[test]
    fn one() {
        let one = Vec3::one();
        let v = Vec3::new([1, 2, 3]);
        assert_eq!(v * one, v);
        assert_eq!(one * v, v);
    }

    #[test]
    fn add() {
        let u = Vec3::new([0., 1., 2.]);
        let v = Vec3::new([2., 1., 0.]);
        assert_eq!(u + v, Vec3::new([2., 2., 2.]));
    }

    #[test]
    fn neg() {
        let u = Vec3::new([1, 2, 3]);
        let v = -u;
        assert_eq!(*v, [-1, -2, -3]);
        assert_eq!(u + v, Vec3::zero());
        assert_eq!(v + u, Vec3::zero());
    }

    #[test]
    fn sub() {
        let u = Vec3::new([1, 2, 3]);
        let v = Vec3::new([3, 2, 1]);
        assert_eq!(u - v, Vec3::new([-2, 0, 2]));
    }

    #[test]
    fn mul_scalar() {
        let v = Vec3::new([1, 2, 3]);
        assert_eq!(*(v * 2), [2, 4, 6]);
    }

    #[test]
    fn mul_vec() {
        let u = Vec3::nth(1, 1);
        let v = Vec3::new([1, 2, 3]);
        assert_eq!(*(u * v), [0, 2, 0]);
    }

    #[test]
    fn div_scalar() {
        let u = Vec3::new([1, 2, 5]);
        assert_eq!(*(u / 2), [0, 1, 2]);
    }

    #[test]
    fn div_vec() {
        let u = Vec3::new([2, 3, 5]);
        let v = Vec3::new([2, 3, 5]);
        assert_eq!(u / v, Vec3::one());
    }

    #[test]
    fn rem() {
        let u = Vec2::new([16, 8]);
        let size = Vec2::new([16, 16]);
        assert_eq!(*(u % size), [0, 8]);
    }
}
