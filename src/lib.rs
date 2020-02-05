//! This is a lightweight matrix / vector library meant for usage with WebGL.
//!
//! At the core this library only consists of the traits `Matrix` and `Vector`. All implementations are optional features and can be added as needed.
//!
//! Available features:
//! * `Matrix4`: 4x4 matrix operations (includes *Vector4*)
//! * `Matrix3`: 3x3 matrix operations (includes *Vector3*)
//! * `Vector4`: 4-dimensional vector operations
//! * `Vector3`: 3-dimensional vector operations
//! * `SliceOps`: Low level slice operations such as addition, subtraction, scaling etc.
//!
//! ## Examples
//!
//! All the types are simple arrays. You may also just use slices as operands.
//!
//! ```rust
//! use webgl_matrix::{Matrix, Vector, ProjectionMatrix, Mat4, Vec4, Mat3, Vec3};
//!
//! fn main() {
//!     // all the default operations available
//!     let mut B = [1., 2., 3.,
//!                  4., 5., 6.,
//!                  7., 8., 9.];
//!     let b = Vec3::ones();
//!     // Matrix operations are in-place
//!     B.inverse();
//!     B.transpose();
//!     // ..
//!
//!     // Some basic vector operations
//!     let c = B.mul_vector_left(&b);
//!     let mag = c.mag(); // magnitude
//!     let d = c.scale(5.);
//!     let e = c.add(&b);
//!
//!     // Or fancier transformations
//!     B.translate(&[1., 2., 3.]);
//!
//!     let A = Mat4::identity();
//!     // operate on slices
//!     let b = [1., 2., 3., 4., 5., 6., 7.];
//!
//!     // with automatic homogenous coordinate expansion
//!     let c = A.mul_vector(&b[0..=2]);
//!     // or using all four coordinates
//!     let d = A.mul_vector(&b[3..]);
//!
//!     // create projection matrices (left, right, bot, top, near, far)
//!     let P = Mat4::create_perspective_from_viewport(0., 1., 0., 1., 0.1, 10.);
//! }
//! ```

mod matrix;
#[macro_use]
mod vector;

#[cfg(feature = "Vector3")]
mod vec3;
#[cfg(feature = "Vector3")]
pub use vec3::Vec3;

#[cfg(feature = "Matrix3")]
mod mat3;
#[cfg(feature = "Matrix3")]
pub use mat3::Mat3;

#[cfg(feature = "Vector4")]
mod vec4;
#[cfg(feature = "Vector4")]
pub use vec4::Vec4;

#[cfg(feature = "Matrix4")]
mod mat4;
#[cfg(feature = "Matrix4")]
pub use mat4::{Mat4, ProjectionMatrix};

#[cfg(any(feature = "Matrix4", feature = "Matrix3"))]
pub use vector::MulVectorMatrix;

#[cfg(feature = "SliceOps")]
pub mod slice_ops;

pub mod utils;
pub use crate::matrix::Matrix;
pub use crate::vector::Vector;

pub mod prelude;
