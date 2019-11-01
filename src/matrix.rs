/// The base Matrix trait
pub trait Matrix {
    type MatrixType;
    type VectorType;

    /// Create a matrix filled with zeros
    fn zeros() -> Self::MatrixType;

    /// Create a matrix filled with ones
    fn ones() -> Self::MatrixType;

    /// Create the identity matrix
    fn identity() -> Self::MatrixType;

    /// Copy values to another matrix
    fn copy_to(&self, dst: &mut Self::MatrixType);

    /// Compute the transpose of this matrix
    fn transpose(&mut self) -> &mut Self::MatrixType;

    /// Perform matrix-multiplication with the given right-hand-side operand
    fn mul(&mut self, rhs: &Self::MatrixType) -> &mut Self::MatrixType;

    /// Multiplies this matrix with the given right-hand-side vector, i.e. `Matrix * rhs`
    ///
    /// Depending on dimensionality, the homogenous coordinate can be omitted,
    /// if so, it will be assumed to be equal to 1.
    fn mul_vector(&self, rhs: &[f32]) -> Self::VectorType;

    /// Multiplies the given row vector with this matrix, i.e. `lhs * Matrix`
    ///
    /// Depending on dimensionality, the homogenous coordinate can be omitted,
    /// if so, it will be assumed to be equal to 1.
    fn mul_vector_left(&self, lhs: &[f32]) -> Self::VectorType;

    /// Perform element-wise addition with the given right-hand-side operand
    fn add(&mut self, rhs: &Self::MatrixType) -> &mut Self::MatrixType;
    /// Perform element-wise substraction with the given right-hand-side operand
    fn sub(&mut self, rhs: &Self::MatrixType) -> &mut Self::MatrixType;

    /// Scale the matrix elment-wise by the given constant
    fn scale(&mut self, factor: f32) -> &mut Self::MatrixType;

    /// Compute the inverse of this matrix. Returns `None` if it is singular.
    fn inverse(&mut self) -> Option<&mut Self::MatrixType>;

    /// Compute the determinant of this matrix.
    fn det(&self) -> f32;

    /// Compute the adjugate of this matrix
    fn adjugate(&mut self) -> &mut Self::MatrixType;

    /// Translate this matrix into the given direction
    ///
    /// Depending on dimensionality, the homogenous coordinate of `direction` can be omitted,
    /// if so, it will be assumed to be equal to 1.
    fn translate(&mut self, direction: &[f32]) -> &mut Self::MatrixType;

    /// Rotate this matrix by the given angle (radians) around the given axis
    ///
    /// Depending on dimensionality, the homogenous coordinate of `axis` can be omitted,
    /// if so, it will be assumed to be equal to 1.
    fn rotate(&mut self, angle: f32, axis: &[f32]) -> &mut Self::MatrixType;
}
