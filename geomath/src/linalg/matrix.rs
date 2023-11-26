use std::ops::{Add, Index, IndexMut, Mul, Sub};

use super::vector::Vector;
use crate::traits::Scalar;
use num::{One, Zero};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Matrix<T: Scalar, const R: usize, const C: usize> {
    rows: [Vector<T, C>; R],
}

impl<T: Scalar, const R: usize, const C: usize> Matrix<T, R, C> {
    fn from_row_array(rows: [[T; C]; R]) -> Self {
        let rows = rows.map(Vector::from);
        Self { rows }
    }

    pub fn with_rows(rows: [Vector<T, C>; R]) -> Self {
        Self { rows }
    }

    fn from_col_array(cols: [[T; R]; C]) -> Self {
        let rows = std::array::from_fn(|r| {
            let row = std::array::from_fn(|c| cols[c][r]);
            Vector::new(row)
        });
        Self { rows }
    }

    pub fn with_cols(cols: [Vector<T, R>; C]) -> Self {
        let rows = std::array::from_fn(|r| {
            let row = std::array::from_fn(|c| cols[c][r]);
            Vector::new(row)
        });
        Self { rows }
    }

    pub fn row(&self, i: usize) -> Vector<T, C> {
        self.rows[i]
    }

    pub fn col(&self, i: usize) -> Vector<T, R> {
        let coef = std::array::from_fn(|r| self.row(r)[i]);
        Vector::new(coef)
    }
}

impl<T: Scalar, const R: usize, const C: usize> Default for Matrix<T, R, C> {
    fn default() -> Self {
        Self::zero()
    }
}

impl<T: Scalar, const R: usize, const C: usize> Index<(usize, usize)> for Matrix<T, R, C> {
    type Output = T;

    fn index(&self, index: (usize, usize)) -> &Self::Output {
        &self.rows[index.0][index.1]
    }
}

impl<T: Scalar, const R: usize, const C: usize> IndexMut<(usize, usize)> for Matrix<T, R, C> {
    fn index_mut(&mut self, index: (usize, usize)) -> &mut Self::Output {
        &mut self.rows[index.0][index.1]
    }
}

impl<T: Scalar, const R: usize, const C: usize> Zero for Matrix<T, R, C> {
    fn zero() -> Self {
        Self {
            rows: [Vector::zero(); R],
        }
    }

    fn is_zero(&self) -> bool {
        self.rows.iter().all(|r| r.is_zero())
    }
}

impl<T: Scalar, const N: usize> One for Matrix<T, N, N> {
    fn one() -> Self {
        let rows = std::array::from_fn(|i| Vector::nth(i, T::one()));
        Self { rows }
    }

    fn set_one(&mut self) {
        for r in 0..N {
            for c in 0..N {
                self[(r, c)] = if r == c { T::one() } else { T::zero() }
            }
        }
    }
}

impl<T: Scalar, const R: usize, const C: usize> Add for Matrix<T, R, C> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        let rows = std::array::from_fn(|i| self.row(i) + rhs.row(i));
        Self { rows }
    }
}

impl<T: Scalar, const R: usize, const C: usize> Sub for Matrix<T, R, C> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        let rows = std::array::from_fn(|i| self.row(i) - rhs.row(i));
        Self { rows }
    }
}

impl<T: Scalar, const R: usize, const C: usize> Mul<T> for Matrix<T, R, C> {
    type Output = Matrix<T, R, C>;

    fn mul(self, rhs: T) -> Self::Output {
        let rows = self.rows.map(|r| r * rhs);
        Self { rows }
    }
}

impl<T: Scalar, const R: usize, const C: usize> Mul<Vector<T, C>> for Matrix<T, R, C> {
    type Output = Vector<T, R>;

    fn mul(self, rhs: Vector<T, C>) -> Self::Output {
        let coef = std::array::from_fn(|i| self.row(i).inner(rhs));
        Vector::new(coef)
    }
}

impl<T: Scalar, const R: usize, const C: usize, const D: usize> Mul<Matrix<T, C, D>>
    for Matrix<T, R, C>
{
    type Output = Matrix<T, R, D>;

    fn mul(self, rhs: Matrix<T, C, D>) -> Self::Output {
        let rows = std::array::from_fn(|r| std::array::from_fn(|c| self.row(r).inner(rhs.col(c))));
        Matrix::<T, R, D>::from_row_array(rows)
    }
}

pub type Mat3<T = f64> = Matrix<T, 3, 3>;
pub type Mat4<T = f64> = Matrix<T, 4, 4>;

#[cfg(test)]
mod tests {
    use crate::linalg::{matrix::Mat3, vector::Vec3};

    #[test]
    fn debug() {
        let m = Mat3::from_row_array([[1, 2, 3], [4, 5, 6], [7, 8, 9]]);
        assert_eq!(format!("{m:?}"), "Matrix { rows: [Vector { coef: [1, 2, 3] }, Vector { coef: [4, 5, 6] }, Vector { coef: [7, 8, 9] }] }");
    }

    #[test]
    fn eq() {
        let m1 = Mat3::with_rows([
            Vec3::new([1, 2, 3]),
            Vec3::new([4, 5, 6]),
            Vec3::new([7, 8, 9]),
        ]);
        let m2 = Mat3::from_row_array([[1, 2, 3], [4, 5, 6], [7, 8, 9]]);
        assert!(m1.eq(&m2));
        assert!(!m1.ne(&m2));
    }

    #[test]
    fn transpose() {
        unimplemented!();
    }

    #[test]
    fn default() {
        unimplemented!();
    }

    #[test]
    fn index() {
        unimplemented!();
    }

    #[test]
    fn index_mut() {
        unimplemented!()
    }

    #[test]
    fn zero() {
        unimplemented!();
    }

    #[test]
    fn one() {
        unimplemented!();
    }

    #[test]
    fn add() {
        unimplemented!();
    }

    #[test]
    fn neg() {
        unimplemented!();
    }

    #[test]
    fn sub() {
        unimplemented!();
    }

    #[test]
    fn mul_scalar() {
        unimplemented!();
    }

    #[test]
    fn mul_vec() {
        unimplemented!();
    }

    #[test]
    fn mul_mat() {
        unimplemented!();
    }

    #[test]
    fn div_scalar() {
        unimplemented!();
    }
}
