//! Core library for kallistors.

pub mod bias;
pub mod ec;
pub mod error;
pub mod index;
pub mod io;
pub mod pseudoalign;
pub mod quant;
pub mod util;

pub use error::Error;
pub type Result<T> = std::result::Result<T, Error>;
