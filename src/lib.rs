mod matrix;

#[cfg(feature = "Matrix3")]
mod mat3;
#[cfg(feature = "Matrix4")]
pub use mat3::Mat3;

#[cfg(feature = "Matrix4")]
mod mat4;
#[cfg(feature = "Matrix4")]
pub use mat4::Mat4;

pub mod utils;
pub use crate::matrix::Matrix;
