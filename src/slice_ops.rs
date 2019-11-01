use std::f32;

#[inline]
/// Performs element-wise multiplication and places the result into `lhs`
///
/// Terminates at the end of the shorter sequence.
pub fn mul(lhs: &mut [f32], rhs: &[f32]) {
    for (l, r) in lhs.iter_mut().zip(rhs.iter()) {
        *l *= r;
    }
}

#[inline]
/// Performs element-wise addition and places the result into `lhs`
///
/// Terminates at the end of the shorter sequence.
pub fn add(lhs: &mut [f32], rhs: &[f32]) {
    for (l, r) in lhs.iter_mut().zip(rhs.iter()) {
        *l += r;
    }
}

#[inline]
/// Performs element-wise substraction and places the result into `lhs`
///
/// Terminates at the end of the shorter sequence.
pub fn sub(lhs: &mut [f32], rhs: &[f32]) {
    for (l, r) in lhs.iter_mut().zip(rhs.iter()) {
        *l -= r;
    }
}

#[inline]
/// Multiplies the given sequence element-wise with the given constant factor
pub fn scale(seq: &mut [f32], factor: f32) {
    for i in seq.iter_mut() {
        *i *= factor;
    }
}

#[inline]
/// Calculates the magnitude of the given sequence, same as sqrt(`mag2(seq)`)
pub fn mag(seq: &[f32]) -> f32 {
    mag2(seq).sqrt()
}

#[inline]
/// Calculates the squared magnitude of the given sequence, i.e. `seq[0] * seq[0] + seq[1] * seq[1] + ...`
pub fn mag2(seq: &[f32]) -> f32 {
    let mut sum = 0.0;
    for i in seq.iter() {
        sum += i * i;
    }
    sum
}

#[inline]
/// Calculates the standard dot product of the two sequences.
///
/// Terminates at the end of the shorter sequence.
pub fn dot(lhs: &[f32], rhs: &[f32]) -> f32 {
    let mut sum = 0.0;
    for (i1, i2) in lhs.iter().zip(rhs.iter()) {
        sum += i1 * i2;
    }
    sum
}
