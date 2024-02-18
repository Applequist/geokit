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
    fn fwd(&self, input: &[f64], output: &mut [f64]) -> Result<()>;

    /// Perform the backward, eg inverse, operation on the given input.
    fn bwd(&self, _input: &[f64], _output: &mut [f64]) -> Result<()> {
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
}

dyn_clone::clone_trait_object!(Operation);

impl Operation for Box<dyn Operation> {
    fn in_dim(&self) -> usize {
        self.as_ref().in_dim()
    }

    fn out_dim(&self) -> usize {
        self.as_ref().out_dim()
    }

    fn fwd(&self, input: &[f64], output: &mut [f64]) -> Result<()> {
        self.as_ref().fwd(input, output)
    }

    fn bwd(&self, input: &[f64], output: &mut [f64]) -> Result<()> {
        self.as_ref().bwd(input, output)
    }
}

pub fn apply_seq<T: Operation>(t: T, input: &[f64], output: &mut [f64]) -> Result<usize> {
    let input_coords = input.chunks_exact(t.in_dim());
    let output_coords = output.chunks_exact_mut(t.out_dim());
    let mut ix: usize = 0;
    for (i, o) in input_coords.zip(output_coords) {
        match t.fwd(i, o) {
            Ok(_) => {
                ix += 1;
            }
            Err(e) => return Err(OperationError::SequenceTransformError(ix, Box::new(e))),
        }
    }
    Ok(ix)
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
    fn fwd(&self, i: &[f64], o: &mut [f64]) -> Result<()> {
        o.copy_from_slice(i);
        Ok(())
    }

    #[inline]
    fn bwd(&self, input: &[f64], output: &mut [f64]) -> Result<()> {
        self.fwd(input, output)?;
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

    fn fwd(&self, input: &[f64], output: &mut [f64]) -> Result<()> {
        self.0.bwd(input, output)
    }

    fn bwd(&self, input: &[f64], output: &mut [f64]) -> Result<()> {
        self.0.fwd(input, output)
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

    fn fwd(&self, i: &[f64], o: &mut [f64]) -> Result<()> {
        let mut os = [0.0; 3];
        self.first.fwd(i, &mut os[0..self.first.out_dim()])?;
        self.then.fwd(&os, o)
    }

    fn bwd(&self, i: &[f64], o: &mut [f64]) -> Result<()> {
        let mut os = [0.0; 3];
        self.then.bwd(i, &mut os[0..self.then.in_dim()])?;
        self.first.bwd(&os, o)
    }
}

pub mod conversion;
pub mod transformation;
