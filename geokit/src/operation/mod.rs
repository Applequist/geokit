use core::f64;
use std::fmt::Debug;

use dyn_clone::DynClone;
use smallvec::SmallVec;
use thiserror::Error;

/// Value returned when an operation cannot be completed.
/// `index` indicates the index of **coordinates** (as opposed to the index of individual floats)
/// in the operation input where the error occured.
/// `reason` is a human readable reason.
/// TODO: Do we want an enum for reason instead ?
#[derive(Error, Debug)]
#[error("Operation failed at input index {}: {}", index, reason)]
pub struct OperationError {
    index: usize,
    reason: String,
}

/// The return type for [`Transformation::apply()`] and [`Transformation::apply_seq()`].
pub type Result<T> = std::result::Result<T, OperationError>;

/// Base trait for coordinate operations that may be invertible.
/// By default, operations are considered not invertible.
pub trait DynOperation {
    /// Return the input coordinates dimension of the forward operation.
    fn fwd_in_dim(&self) -> usize;

    /// Return the output coordinates dimension of the forward operation.
    fn fwd_out_dim(&self) -> usize;

    /// Perform the forward operation on the given input.
    fn fwd(&self, input: &[f64], output: &mut [f64]) -> Result<()>;

    /// Return whether this operation is invertible.
    fn is_invertible(&self) -> bool {
        false
    }

    /// Perform the backward, eg inverse, operation on the given input.
    fn bwd(&self, _input: &[f64], _output: &mut [f64]) -> Result<()> {
        Err(OperationError {
            index: 0,
            reason: "Operation is not invertible".into(),
        })
    }
}

/// Base trait for **unidirectional** transformation.
pub trait Operation: DynClone {
    /// Return the input coordinates dimension.
    fn in_dim(&self) -> usize;

    /// Return the output coordinates dimension.
    fn out_dim(&self) -> usize;

    /// Perform the operation on the given input.
    fn apply(&self, input: &[f64], output: &mut [f64]) -> Result<()>;

    /// Apply the operation on coordinates in the `input` slice, storing the transformed
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
                Err(OperationError { index: _, reason }) => {
                    return Err(OperationError { index, reason })
                }
            };
        }
        Ok(index)
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

/// A [DynOperation] wrapper that selects the forward operation.
pub struct Fwd<T>(pub T);

impl<T> Debug for Fwd<T>
where
    T: Debug,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_tuple("Fwd").field(&self.0).finish()
    }
}

impl<T> Clone for Fwd<T>
where
    T: Clone,
{
    fn clone(&self) -> Self {
        Self(self.0.clone())
    }
}

impl<T> Operation for Fwd<T>
where
    T: DynOperation + Clone,
{
    fn in_dim(&self) -> usize {
        self.0.fwd_in_dim()
    }

    fn out_dim(&self) -> usize {
        self.0.fwd_out_dim()
    }

    fn apply(&self, input: &[f64], output: &mut [f64]) -> Result<()> {
        self.0.fwd(input, output)
    }
}

/// A [DynOperation] wrapper that selects the backward, eg inverse, operation.
pub struct Bwd<T>(pub T);

impl<T> Debug for Bwd<T>
where
    T: Debug,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_tuple("Bwd").field(&self.0).finish()
    }
}

impl<T> Clone for Bwd<T>
where
    T: Clone,
{
    fn clone(&self) -> Self {
        Self(self.0.clone())
    }
}

impl<T> Operation for Bwd<T>
where
    T: DynOperation + Clone,
{
    fn in_dim(&self) -> usize {
        self.0.fwd_out_dim()
    }

    fn out_dim(&self) -> usize {
        self.0.fwd_in_dim()
    }

    fn apply(&self, input: &[f64], output: &mut [f64]) -> Result<()> {
        self.0.bwd(input, output)
    }
}

/// A chained coordinates operation.
pub struct Chain<A, B> {
    first: A,
    then: B,
}

impl<A, B> Debug for Chain<A, B>
where
    A: Debug,
    B: Debug,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Chain")
            .field("first", &self.first)
            .field("then", &self.then)
            .finish()
    }
}

impl<A, B> Clone for Chain<A, B>
where
    A: Clone,
    B: Clone,
{
    fn clone(&self) -> Self {
        Self {
            first: self.first.clone(),
            then: self.then.clone(),
        }
    }
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

    fn apply(&self, i: &[f64], o: &mut [f64]) -> Result<()> {
        let mut os: SmallVec<[f64; 3]> = SmallVec::from_elem(0.0, self.first.out_dim());
        self.first.apply(i, &mut os)?;
        self.then.apply(&os, o)
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

    fn apply(&self, i: &[f64], o: &mut [f64]) -> Result<()> {
        o.copy_from_slice(i);
        Ok(())
    }
}

pub fn identity<const N: usize>() -> impl Operation {
    Identity::<N>
}

pub mod conversion;
pub mod transformation;