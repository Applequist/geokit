use core::f64;
use std::fmt::Debug;

use dyn_clone::DynClone;
use thiserror::Error;

#[derive(Error, Debug)]
pub enum OperationError {
    #[error("Domain Error")]
    DomainError,
    #[error("Not Invertible")]
    NotInvertible,
    #[error("Sequence Error")]
    SequenceTransformError(usize, Box<OperationError>),
}

/// The return type for [`Transformation::apply()`] and [`Transformation::apply_seq()`].
pub type Result<T> = std::result::Result<T, OperationError>;

/// Base trait for coordinate operations that may be invertible.
/// By default, operations are considered not invertible.
pub trait Operation: DynClone {
    /// Returns the input coordinates dimension of the forward operation.
    fn in_dim(&self) -> usize;

    /// Return the output coordinates dimension of the forward operation.
    fn out_dim(&self) -> usize;

    /// Perform the forward operation on the given input.
    fn apply_fwd(&self, input: &[f64], output: &mut [f64]) -> Result<()>;

    /// Perform the backward, eg inverse, operation on the given input.
    fn apply_bwd(&self, _input: &[f64], _output: &mut [f64]) -> Result<()> {
        Err(OperationError::NotInvertible)
    }

    fn boxed(self) -> Box<dyn Operation>
    where
        Self: Sized + 'static,
    {
        Box::new(self)
    }

    /// Chain this operation with the given one.
    /// It returns a new operation that perform this operation,
    /// followed by the given one.
    fn and_then<T>(self, t: T) -> Chain<Self, T>
    where
        T: Operation,
        Self: Sized,
    {
        Chain {
            first: self,
            then: t,
        }
    }

    /// Apply the fwd transformation into the given mutable slice.
    fn fwd(&self, input: &[f64], output: &mut [f64]) -> Result<usize> {
        assert!(
            input.len() % self.in_dim() == 0,
            "For forward transformation, the input sequence length must be a multiple of the operation's input dimension."
        );
        assert!(
            output.len() % self.out_dim() == 0,
            "For forward transformation, the output sequence length must be a multiple of the operation's output dimension."
        );
        let input_coords = input.chunks_exact(self.in_dim());
        let output_coords = output.chunks_exact_mut(self.out_dim());
        let mut ix: usize = 0;
        for (i, o) in input_coords.zip(output_coords) {
            match self.apply_fwd(i, o) {
                Ok(_) => {
                    ix += 1;
                }
                Err(e) => return Err(OperationError::SequenceTransformError(ix, Box::new(e))),
            }
        }
        Ok(ix)
    }

    /// Apply the bwd transformation into the given mutable slice.
    fn bwd(&self, input: &[f64], output: &mut [f64]) -> Result<usize> {
        assert!(
            input.len() % self.out_dim() == 0,
            "For backward transformation, the input sequence length must be a multiple of the operation's output dimension."
        );
        assert!(
            output.len() % self.in_dim() == 0,
            "For backward transformation, the output sequence length must be a multiple of the operation's input dimension."
        );
        let input_coords = input.chunks_exact(self.out_dim());
        let output_coords = output.chunks_exact_mut(self.in_dim());
        let mut ix: usize = 0;
        for (i, o) in input_coords.zip(output_coords) {
            match self.apply_bwd(i, o) {
                Ok(_) => {
                    ix += 1;
                }
                Err(e) => return Err(OperationError::SequenceTransformError(ix, Box::new(e))),
            }
        }
        Ok(ix)
    }
    /// Apply the fwd transformation to the input and stores the output into a newly allocated vector.
    fn fwd_new(&self, input: &[f64]) -> Result<Vec<f64>> {
        assert!(
            input.len() % self.in_dim() == 0,
            "For forward transformation, the input sequence length must be a multiple of the operation's input dimension."
        );
        let capacity = input.len() / self.in_dim() * self.out_dim();
        let mut output = vec![0.0; capacity];
        self.fwd(input, &mut output)?;
        Ok(output)
    }

    fn bwd_new(&self, input: &[f64]) -> Result<Vec<f64>> {
        assert!(input.len() % self.out_dim() == 0,
            "For backward transformation, the input sequence length must be a multiple of the operation's output dimension.",
        );
        let capacity = input.len() / self.out_dim() * self.in_dim();
        let mut output = vec![0.0; capacity];
        self.bwd(input, &mut output)?;
        Ok(output)
    }
}

dyn_clone::clone_trait_object!(Operation);

impl Operation for Box<dyn Operation> {
    fn in_dim(&self) -> usize {
        self.as_ref().in_dim()
    }

    fn out_dim(&self) -> usize {
        self.as_ref().out_dim()
    }

    fn apply_fwd(&self, input: &[f64], output: &mut [f64]) -> Result<()> {
        self.as_ref().apply_fwd(input, output)
    }

    fn apply_bwd(&self, input: &[f64], output: &mut [f64]) -> Result<()> {
        self.as_ref().apply_bwd(input, output)
    }
}

/// A *dummy* operation that simply copies its input into the output, eg a no-op operation.
#[derive(Debug, Clone)]
struct Identity<const N: usize>;

impl<const N: usize> Operation for Identity<N> {
    fn in_dim(&self) -> usize {
        N
    }

    fn out_dim(&self) -> usize {
        N
    }

    #[inline]
    fn apply_fwd(&self, i: &[f64], o: &mut [f64]) -> Result<()> {
        o.copy_from_slice(i);
        Ok(())
    }

    #[inline]
    fn apply_bwd(&self, input: &[f64], output: &mut [f64]) -> Result<()> {
        self.apply_fwd(input, output)?;
        Ok(())
    }
}

pub fn identity<const N: usize>() -> impl Operation {
    Identity::<N>
}

/// A [DynOperation] wrapper that selects the backward, eg inverse, operation.
#[derive(Debug, Clone)]
pub struct Inv<T>(pub T);

impl<T> Operation for Inv<T>
where
    T: Operation + Clone,
{
    fn in_dim(&self) -> usize {
        self.0.out_dim()
    }

    fn out_dim(&self) -> usize {
        self.0.in_dim()
    }

    fn apply_fwd(&self, input: &[f64], output: &mut [f64]) -> Result<()> {
        self.0.apply_bwd(input, output)
    }

    fn apply_bwd(&self, input: &[f64], output: &mut [f64]) -> Result<()> {
        self.0.apply_fwd(input, output)
    }
}

/// A chained coordinates operation.
#[derive(Debug, Clone)]
pub struct Chain<A, B> {
    first: A,
    then: B,
}

impl<A, B> Operation for Chain<A, B>
where
    A: Operation + Clone,
    B: Operation + Clone,
{
    fn in_dim(&self) -> usize {
        self.first.in_dim()
    }

    fn out_dim(&self) -> usize {
        self.then.out_dim()
    }

    fn apply_fwd(&self, i: &[f64], o: &mut [f64]) -> Result<()> {
        let mut os = [0.0; 3];
        self.first.apply_fwd(i, &mut os[0..self.first.out_dim()])?;
        self.then.apply_fwd(&os, o)
    }

    fn apply_bwd(&self, i: &[f64], o: &mut [f64]) -> Result<()> {
        let mut os = [0.0; 3];
        self.then.apply_bwd(i, &mut os[0..self.then.in_dim()])?;
        self.first.apply_bwd(&os, o)
    }
}

pub mod conversion;
pub mod transformation;

#[cfg(test)]
mod tests {
    use crate::operation::Operation;

    use super::identity;

    #[test]
    fn test_apply_seq() {
        let input = [1.0; 6];
        let mut output = [0.0; 6];
        identity::<3>().fwd(&input, &mut output).unwrap();
        assert_eq!(input, output);
    }
}
