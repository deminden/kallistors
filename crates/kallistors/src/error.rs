use std::path::PathBuf;

/// Crate-wide error type.
#[derive(thiserror::Error, Debug)]
pub enum Error {
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),
    #[error("invalid format: {0}")]
    InvalidFormat(String),
    #[error("unsupported feature: {0}")]
    UnsupportedFeature(String),
    #[error("missing file: {0}")]
    MissingFile(PathBuf),
}
