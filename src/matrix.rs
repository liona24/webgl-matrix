pub trait Matrix {
    type Type;

    /// Create a matrix filled with zeros
    fn zeros() -> Self::Type;

    /// Create a matrix filled with ones
    fn ones() -> Self::Type;

    /// Create the identity matrix
    fn identity() -> Self::Type;

    /// Copy values to another matrix
    fn copy_to(&self, dst : &mut Self::Type);

    /// Compute the transpose of this matrix
    fn transpose(self) -> Self::Type;

    /// Perform matrix-multiplication with the given right-hand-side operand
    fn mul(self, rhs : &Self::Type) -> Self::Type;

    /// Perform element-wise addition with the given right-hand-side operand
    fn add(self, rhs : &Self::Type) -> Self::Type;
    /// Perform element-wise substraction with the given right-hand-side operand
    fn sub(self, rhs : &Self::Type) -> Self::Type;

    /// Scale the matrix elment-wise by the given constant
    fn scale(self, factor : f32) -> Self::Type;

    /// Compute the inverse of this matrix. Returns `None` if it is singular.
    fn inverse(self) -> Option<Self::Type>;

    /// Compute the determinant of this matrix.
    fn det(&self) -> f32;

    /// Compute the adjugate of this matrix
    fn adjugate(self) -> Self::Type;
}