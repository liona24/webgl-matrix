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
    fn copy_to(&self, dst : &mut Self::MatrixType);

    /// Compute the transpose of this matrix
    fn transpose(self) -> Self::MatrixType;

    /// Perform matrix-multiplication with the given right-hand-side operand
    fn mul(self, rhs : &Self::MatrixType) -> Self::MatrixType;

    /// Multiplies this matrix with the given right-hand-side vector, i.e. `Matrix * rhs`
    fn mul_vector(&self, rhs : &Self::VectorType) -> Self::VectorType;

    /// Multiplies the given row vector with this matrix, i.e. `lhs * Matrix`
    fn mul_vector_left(&self, lhs : &Self::VectorType) -> Self::VectorType;

    /// Perform element-wise addition with the given right-hand-side operand
    fn add(self, rhs : &Self::MatrixType) -> Self::MatrixType;
    /// Perform element-wise substraction with the given right-hand-side operand
    fn sub(self, rhs : &Self::MatrixType) -> Self::MatrixType;

    /// Scale the matrix elment-wise by the given constant
    fn scale(self, factor : f32) -> Self::MatrixType;

    /// Compute the inverse of this matrix. Returns `None` if it is singular.
    fn inverse(self) -> Option<Self::MatrixType>;

    /// Compute the determinant of this matrix.
    fn det(&self) -> f32;

    /// Compute the adjugate of this matrix
    fn adjugate(self) -> Self::MatrixType;

    /// Translate this matrix into the given direction
    fn translate(self, direction : &Self::VectorType) -> Self::MatrixType;
}