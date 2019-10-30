mod matrix;
mod vector;

#[cfg(feature = "Matrix3")]
mod mat3;
#[cfg(feature = "Matrix3")]
pub use mat3::Mat3;

#[cfg(feature = "Matrix4")]
mod mat4;
#[cfg(feature = "Matrix4")]
mod vec4;
#[cfg(feature = "Matrix4")]
pub use mat4::{Mat4, ProjectionMatrix};
#[cfg(feature = "Matrix4")]
pub use vec4::Vec4;

pub mod utils;
pub use crate::matrix::Matrix;
pub use crate::vector::Vector;
