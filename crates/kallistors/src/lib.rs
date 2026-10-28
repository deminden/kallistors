//! Core library for kallistors.
//!
//! `index` reads and builds the on-disk graph; `pseudoalign` turns reads into
//! counted transcript equivalence classes; `quant` estimates abundance from
//! those counts. Command-specific input and output policy lives in the binary.

pub mod bias;
mod dna;
pub mod ec;
pub mod error;
pub mod index;
pub mod io;
pub mod pseudoalign;
pub mod quant;
pub mod timing;
pub mod util;

pub use error::Error;
pub type Result<T> = std::result::Result<T, Error>;

mod minimizer_hash;
