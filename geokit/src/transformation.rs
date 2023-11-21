use smallvec::SmallVec;
use std::fmt::Debug;
use thiserror::Error;

/// Value returned when [`Transformation::apply()`] cannot be completed.
/// `index` indicates the index of **coordinates** (as opposed to the index of individual floats)
/// in the transformation input where the error occured.
/// `reason` is a human readable reason.
/// TODO: Do we want an enum for reason instead ?
#[derive(Error, Debug)]
#[error("Transformation failed at input index {}: {}", index, reason)]
pub struct TransformationError {
    index: usize,
    reason: String,
}

/// The return type for [`Transformation::apply()`] and [`Transformation::apply_seq()`].
pub type Result<T> = std::result::Result<T, TransformationError>;

/// A trait for coordinates transformations.
pub trait Transformation: Debug {
    /// Return the dimension of the input coordinates.
    fn in_dim(&self) -> usize;
    /// Return the dimension of the output coordinates.
    fn out_dim(&self) -> usize;

    /// Apply this transformation once to the coordinates from the `input` slice and store the transformed
    /// coordinates in the `output` slice.
    ///
    /// # Panics
    ///
    /// This function will panic if the `input` length is less than [`in_dim`][Self::in_dim()] or
    /// if `output` length is less than [`out_dim`][Self::out_dim()].
    fn apply(&self, input: &[f64], output: &mut [f64]) -> Result<()>;

    /// Apply the transformation on coordinates in the `input` slice, storing the transformed
    /// coordinates into the `output` slice, and returning the number of successful applications when either:
    /// * the input slice has been fully processed
    /// * or the output slice has been filled,
    ///
    /// # Errors
    ///
    /// If an error occurs while transforming the sequence, this function returns an [`TransformationError`] containing:
    /// * the index of the transformation application, eg 0 on 1st application, 1 on 2nd...
    /// * the reason for the error
    ///
    /// If the index is greater that 0, the previous transformation results are avalaible in `output`.
    fn apply_seq(&self, input: &[f64], output: &mut [f64]) -> Result<usize> {
        let input_chunks = input.chunks_exact(self.in_dim());
        let output_chunks = output.chunks_exact_mut(self.out_dim());
        let mut index: usize = 0;
        for (i, o) in input_chunks.zip(output_chunks) {
            match self.apply(i, o) {
                Ok(_) => index += 1,
                Err(TransformationError { index: _, reason }) => {
                    return Err(TransformationError { index, reason })
                }
            };
        }
        Ok(index)
    }

    fn boxed(&self) -> Box<dyn Transformation>;

    fn and_then<T>(self, t: T) -> Chain<Self, T>
    where
        T: Transformation,
        Self: Sized,
    {
        Chain {
            first: self,
            then: t,
        }
    }
}

/// A chained [`Transformation`].
#[derive(Debug)]
pub struct Chain<A, B> {
    first: A,
    then: B,
}

impl<A, B> Transformation for Chain<A, B>
where
    A: Transformation + Clone + 'static,
    B: Transformation + Clone + 'static,
{
    fn in_dim(&self) -> usize {
        self.first.in_dim()
    }

    fn out_dim(&self) -> usize {
        self.then.out_dim()
    }

    fn apply(&self, input: &[f64], output: &mut [f64]) -> Result<()> {
        let mut middle: SmallVec<[f64; 3]> = SmallVec::from_elem(0., self.first.out_dim());
        self.first.apply(input, &mut middle)?;
        self.then.apply(&middle, output)
    }

    fn boxed(&self) -> Box<dyn Transformation> {
        Box::new(Chain {
            first: self.first.clone(),
            then: self.then.clone(),
        })
    }
}

impl<A, B> Clone for Chain<A, B>
where
    A: Clone,
    B: Clone,
{
    fn clone(&self) -> Self {
        Chain {
            first: self.first.clone(),
            then: self.then.clone(),
        }
    }
}

impl Clone for Box<dyn Transformation> {
    fn clone(&self) -> Box<dyn Transformation> {
        Transformation::boxed(&**self)
    }
}

#[derive(Debug, Clone, Copy)]
pub struct Identity<const N: usize>;

impl<const N: usize> Transformation for Identity<N> {
    fn in_dim(&self) -> usize {
        N
    }

    fn out_dim(&self) -> usize {
        N
    }

    fn apply(&self, input: &[f64], output: &mut [f64]) -> Result<()> {
        output.copy_from_slice(input);
        Ok(())
    }

    fn boxed(&self) -> Box<dyn Transformation> {
        Box::new(Identity::<N>)
    }
}

pub type ToOrd = (usize, f64);

#[derive(Debug, Clone)]
pub struct CoordScaling(pub Vec<ToOrd>);

impl Transformation for CoordScaling {
    fn in_dim(&self) -> usize {
        self.0.len()
    }

    fn out_dim(&self) -> usize {
        3
    }

    fn apply(&self, input: &[f64], output: &mut [f64]) -> Result<()> {
        output.copy_from_slice(&[0.; 3]);
        for (ox, (ix, f)) in self.0.iter().enumerate() {
            output[ox] = input[*ix] * f;
        }
        Ok(())
    }

    fn boxed(&self) -> Box<dyn Transformation> {
        Box::new(self.clone())
    }
}

impl From<(usize, f64)> for CoordScaling {
    fn from((dim, value): (usize, f64)) -> Self {
        let mut ord: Vec<ToOrd> = Vec::new();
        for ix in 0..dim {
            ord.push((ix, value))
        }
        Self(ord)
    }
}

#[cfg(test)]
mod tests {

    use super::*;
    use std::cell::Cell;

    #[derive(Debug, Clone)]
    struct Dummy<const N: usize> {
        max: usize,
        count: Cell<usize>,
    }

    impl<const N: usize> Transformation for Dummy<N> {
        fn in_dim(&self) -> usize {
            N
        }
        fn out_dim(&self) -> usize {
            N
        }
        fn apply(&self, input: &[f64], output: &mut [f64]) -> Result<()> {
            if self.count.get() >= self.max {
                Err(TransformationError {
                    index: 0,
                    reason: "Boom!".into(),
                })
            } else {
                output.copy_from_slice(input);
                let c = self.count.get();
                self.count.set(c + 1);
                Ok(())
            }
        }

        fn boxed(&self) -> Box<dyn Transformation> {
            Box::new(Dummy::<N> {
                max: self.max,
                count: Cell::new(0),
            })
        }
    }

    #[test]
    #[should_panic]
    fn apply_panic() {
        let t = Identity::<2>;
        let _ = t.apply(&[1., 2.], &mut [0.]);
    }

    #[test]
    fn apply_seq_short_output() {
        let t = Identity::<2>;
        let r = t.apply_seq(&[1., 1., 2., 2.], &mut [0., 0.]);
        assert!(matches!(r, Ok(1)), "Expected 1");
    }

    #[test]
    fn apply_seq_short_input() {
        let t = Identity::<2>;
        let r = t.apply_seq(&[1., 1., 2., 2.], &mut [0.; 10]);
        assert!(matches!(r, Ok(2)), "Expected 2");
    }

    #[test]
    fn apply_seq_error_index() {
        let t = Dummy::<2> {
            max: 2,
            count: Cell::new(0),
        };
        let i = &[1., 1., 2., 2., 3., 3.];
        let o = &mut [0.; 10];
        let r = t.apply_seq(i, o);
        assert!(matches!(
            r,
            Err(TransformationError {
                index,
                reason,
            }) if index == 2 && &reason == "Boom!"
        ));
    }

    #[test]
    fn box_clone() {
        let mut out = [0.; 3];

        let boxed_transfo: Box<dyn Transformation> = Box::new(Identity::<3>);
        let _cloned_as_transfo: Box<dyn Transformation> = boxed_transfo.clone();
        _cloned_as_transfo.apply(&[1., 2., 3.], &mut out).unwrap();
        assert_eq!(
            out,
            [1., 2., 3.],
            "Expected {:?}. Got {:?}.",
            [1., 2., 3.],
            out
        );
    }

    #[test]
    fn chain() {
        let first = Dummy::<2> {
            max: 2,
            count: Cell::new(0),
        };
        let then = Dummy::<2> {
            max: 2,
            count: Cell::new(0),
        };
        let chain = first.and_then(then);

        let input = [0.; 4];
        let mut output = [0.; 4];
        let count = chain.apply_seq(&input, &mut output).unwrap();
        assert_eq!(count, 2);
    }
}
