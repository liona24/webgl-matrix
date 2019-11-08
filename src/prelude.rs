#[cfg(feature = "SliceOps")]
pub use crate::slice_ops;

#[cfg(feature = "Vector3")]
pub use crate::Vec3;
#[cfg(feature = "Vector4")]
pub use crate::Vec4;
#[cfg(any(feature = "Vector3", feature = "Vector4"))]
pub use crate::Vector;

#[cfg(feature = "Matrix3")]
pub use crate::Mat3;
#[cfg(feature = "Matrix4")]
pub use crate::{Mat4, ProjectionMatrix};
#[cfg(any(feature = "Matrix3", feature = "Matrix4"))]
pub use crate::{Matrix, MulVectorMatrix};
