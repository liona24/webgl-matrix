/// The base Vector trait
///
/// Note that vector operations are permitted on slices.
/// This is useful for WebGl vertex storages, where many vectors are stored in
/// a single array. The caveat obviously is the potential for errors. Sadly
/// there is no secure alternative for handling this use-case.
pub trait Vector {
    type VectorType;

    /// Create a vector filled with zeros
    fn zeros() -> Self::VectorType;

    /// Create a vector filled with ones
    fn ones() -> Self::VectorType;

    /// Perform element-wise multiplication with the given right-hand-side operand
    fn mul(&self, rhs: &[f32]) -> Self::VectorType;

    /// Perform element-wise addition with the given right-hand-side operand
    fn add(&self, rhs: &[f32]) -> Self::VectorType;
    /// Perform element-wise substraction with the given right-hand-side operand
    fn sub(&self, rhs: &[f32]) -> Self::VectorType;

    /// Scale the vector elment-wise by the given constant
    fn scale(&self, factor: f32) -> Self::VectorType;

    /// Calculate the magnitude of this vector
    fn mag(&self) -> f32;

    /// Calculate the squared magnitude of this vector
    fn mag2(&self) -> f32;

    /// Calculate the dot product of this vector and the given right-hand-side operand
    fn dot(&self, rhs: &[f32]) -> f32;
}

#[cfg(any(feature = "Matrix4", feature = "Matrix3"))]
/// Adds matrix operations to vector types.
pub trait MulVectorMatrix<Matrix> {
    type VectorType;

    /// Interprets `self` as a column vector and multiplies the given matrix
    /// from the left-hand-side, i.e. `lhs * self`
    fn mul_matrix_left(&self, lhs: &Matrix) -> Self::VectorType;

    /// Interprets `self` as a row vector and multiplies the given matrix
    /// from the right-hand-side, i.e. `self * rhs`
    fn mul_matrix(&self, rhs: &Matrix) -> Self::VectorType;
}

macro_rules! impl_vector {
    ($type:ty, $n:expr) => {
        impl Vector for $type {
            type VectorType = $type;
            fn zeros() -> $type {
                [0.; $n]
            }

            fn ones() -> $type {
                [1.; $n]
            }

            fn mul(&self, rhs: &[f32]) -> $type {
                let mut dst = *self;
                mul(&mut dst, rhs);
                dst
            }

            fn add(&self, rhs: &[f32]) -> $type {
                let mut dst = *self;
                add(&mut dst, rhs);
                dst
            }

            fn sub(&self, rhs: &[f32]) -> $type {
                let mut dst = *self;
                sub(&mut dst, rhs);
                dst
            }

            fn scale(&self, factor: f32) -> $type {
                let mut dst = *self;
                scale(&mut dst, factor);
                dst
            }

            fn mag(&self) -> f32 {
                mag(self)
            }

            fn mag2(&self) -> f32 {
                mag2(self)
            }

            fn dot(&self, rhs: &[f32]) -> f32 {
                dot(self, rhs)
            }
        }
    };
}
