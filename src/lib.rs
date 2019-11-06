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

#[cfg(feature = "SliceOps")]
pub mod slice_ops;

pub mod utils;
pub use crate::matrix::Matrix;
pub use crate::vector::Vector;
