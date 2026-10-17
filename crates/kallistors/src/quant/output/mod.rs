//! Kallisto-compatible abundance and run-metadata output.

mod hdf5;
mod metadata;
mod tsv;

pub use hdf5::{H5QuantOutput, write_abundance_h5};
pub use metadata::RunInfo;
pub use tsv::{
    write_abundance_tsv, write_abundance_tsv_from_parts, write_abundance_tsv_from_parts_threaded,
};
