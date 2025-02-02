use super::Float;

pub fn remainder(val: Float, r: Float) -> Float {
    val - (val / r).round_ties_even() * r
}

pub fn iter_fn<const D: usize>(
    init: [Float; D],
    deltas: &dyn Fn([Float; D]) -> [Float; D],
    epsilon: [Float; D],
) -> [Float; D] {
    let mut res = init;
    loop {
        let ds = deltas(res);
        for i in 0..D {
            res[i] += ds[i];
        }
        if ds.iter().zip(epsilon.iter()).all(|(d, e)| d.abs() < *e) {
            break;
        }
    }
    res
}

/// Compute `sum(0, N, c[k] * F[k](x))` where `F[k](x) = sin(2. * k * x)`
/// using Clenshaw's recurrence formula.
/// Since `F[k+1](x) = 2*cos(2 * x) * F[k](x) - F[k-1](x)`,
/// `alpha[k](x) = 2 * cos(2 * x)` and `beta[k](x) = -1`
/// And the sum is:
/// `beta[1](x) * F[0](x) * y[2] + F[1](x) * y[1] + F[0](x) * c[0]`
/// Where
/// `y[N+2] = y[N+1] = 0` and
/// `y[k] = alpha[n](x) * y[k+1] + beta[k+1](x) * y[k+2] + c[k]`
pub fn sin_sum(c_k: &[Float], x: Float) -> Float {
    let alpha = 2. * (2. * x).cos();
    let beta = -1.;

    let mut y_n_2 = 0.;
    let mut y_n_1 = 0.;
    let mut y_n = 0.;
    for ck in c_k.iter().rev() {
        y_n_2 = y_n_1;
        y_n_1 = y_n;
        y_n = alpha * y_n_1 + beta * y_n_2 + ck;
    }
    // F[0](x) = 0.
    // F[1](x) = sin(2x)
    (2. * x).sin() * y_n_1
}

#[cfg(test)]
mod test {
    #[test]
    fn wrap() {
        use super::remainder;
        assert_eq!(remainder(-180.0, 360.0), -180.0);
        assert_eq!(remainder(0.0, 360.0), 0.0);
        assert_eq!(remainder(180.0, 360.0), 180.0);
        assert_eq!(remainder(360.0, 360.0), 0.0);
    }

    #[test]
    fn sin_sum() {
        unimplemented!()
    }
}
