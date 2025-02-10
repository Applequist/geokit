use crate::{
    quantities::Convertible,
    units::{
        angle::{AngleUnit, DEG},
        length::{LengthUnit, M},
    },
};
use approx::AbsDiffEq;
use std::fmt::Display;

/// Check that two quantities of kind [Angle] are equal within the given `err` bound.
pub fn approx_angle_eq<T>(res: T, exp: T, err: <T as AbsDiffEq>::Epsilon, err_msg: &str) -> bool
where
    T: AbsDiffEq + Display,
    <T as AbsDiffEq>::Epsilon: Copy + Convertible<Unit = AngleUnit>,
{
    if !res.abs_diff_eq(&exp, err) {
        println!("{}: {} != {} +/- {:e} deg", err_msg, res, exp, err.val(DEG));
        false
    } else {
        true
    }
}

/// Check that two quantities of kind [length] are equal within the given `err` bound.
pub fn approx_length_eq<T>(res: T, exp: T, err: <T as AbsDiffEq>::Epsilon, err_msg: &str) -> bool
where
    T: AbsDiffEq + Display,
    <T as AbsDiffEq>::Epsilon: Copy + Convertible<Unit = LengthUnit>,
{
    if !res.abs_diff_eq(&exp, err) {
        println!("{}: {} != {} +/- {:e} m", err_msg, res, exp, err.val(M));
        false
    } else {
        true
    }
}
