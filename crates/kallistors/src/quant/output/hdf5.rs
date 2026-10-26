//! HDF5 datasets, string encoding and compression layout.

use std::borrow::Cow;
use std::fs;
use std::path::Path;

use hdf5_pure::{CharacterSet, DatasetBuilder, Datatype, FileBuilder, StringPadding};

use super::super::QuantResult;
use crate::{Error, Result};

const KALLISTO_BIAS_LEN: usize = 4096;

pub struct H5QuantOutput<'a> {
    pub target_names: &'a [&'a str],
    pub target_lengths: &'a [u32],
    pub result: &'a QuantResult,
    pub bootstrap_est_counts: &'a [&'a [f64]],
    pub reads_processed: u64,
    pub fragment_length_hist: Option<&'a [u32]>,
    pub bias_observed: Option<&'a [u32]>,
    pub index_version: u64,
    pub kallisto_version: &'a str,
    pub call: &'a str,
    pub start_time: &'a str,
}

pub fn write_abundance_h5(path: &Path, output: H5QuantOutput<'_>) -> Result<()> {
    let mut file = FileBuilder::new();
    let mut aux = file.create_group("aux");

    write_h5_f64(file.create_dataset("est_counts"), &output.result.est_counts);
    write_h5_i32(
        aux.create_dataset("num_bootstrap"),
        &[
            i32::try_from(output.bootstrap_est_counts.len()).map_err(|_| {
                Error::InvalidFormat("too many bootstrap samples for HDF5 output".into())
            })?,
        ],
    );
    write_h5_i32(
        aux.create_dataset("num_processed"),
        &[i32::try_from(output.reads_processed).map_err(|_| {
            Error::InvalidFormat("n_processed exceeds kallisto-compatible HDF5 range".into())
        })?],
    );
    write_h5_i32(
        aux.create_dataset("fld"),
        &output
            .fragment_length_hist
            .map(u32_slice_to_i32)
            .transpose()?
            .unwrap_or_default(),
    );
    let bias_observed = match output.bias_observed {
        Some(values) => Cow::Owned(u32_slice_to_i32(values)?),
        None => Cow::Owned(vec![1; KALLISTO_BIAS_LEN]),
    };
    write_h5_i32(aux.create_dataset("bias_observed"), &bias_observed);

    let bias_normalized = match output.result.post_bias.as_deref() {
        Some(values) => Cow::Borrowed(values),
        None => Cow::Owned(vec![1.0; KALLISTO_BIAS_LEN]),
    };
    write_h5_f64(aux.create_dataset("bias_normalized"), &bias_normalized);
    write_h5_string_vec(
        aux.create_dataset("kallisto_version"),
        &[output.kallisto_version],
    );
    write_h5_i32(
        aux.create_dataset("index_version"),
        &[i32::try_from(output.index_version).map_err(|_| {
            Error::InvalidFormat("index version exceeds kallisto-compatible HDF5 range".into())
        })?],
    );
    write_h5_string_vec(aux.create_dataset("call"), &[output.call]);
    write_h5_string_vec(aux.create_dataset("start_time"), &[output.start_time]);
    write_h5_string_vec(aux.create_dataset("ids"), output.target_names);
    write_h5_f64(
        aux.create_dataset("eff_lengths"),
        &output.result.eff_lengths,
    );
    write_h5_i32(
        aux.create_dataset("lengths"),
        &u32_slice_to_i32(output.target_lengths)?,
    );
    file.add_group(aux.finish());

    if !output.bootstrap_est_counts.is_empty() {
        let mut bootstrap = file.create_group("bootstrap");
        for (idx, counts) in output.bootstrap_est_counts.iter().enumerate() {
            write_h5_f64(bootstrap.create_dataset(&format!("bs{idx}")), counts);
        }
        file.add_group(bootstrap.finish());
    }

    let bytes = file.finish().map_err(h5_error)?;
    fs::write(path, bytes)?;
    Ok(())
}

fn write_h5_f64(dataset: &mut DatasetBuilder, values: &[f64]) {
    dataset
        .with_f64_data(values)
        .with_shape(&[values.len() as u64]);
    set_h5_chunking(dataset, values.len());
}

fn write_h5_i32(dataset: &mut DatasetBuilder, values: &[i32]) {
    dataset
        .with_i32_data(values)
        .with_shape(&[values.len() as u64]);
    set_h5_chunking(dataset, values.len());
}

fn write_h5_string_vec(dataset: &mut DatasetBuilder, values: &[&str]) {
    let width = values
        .iter()
        .map(|value| value.len())
        .max()
        .unwrap_or(0)
        .saturating_add(1)
        .max(1);
    dataset
        .with_compound_data(
            fixed_ascii_type(width),
            fixed_string_bytes(values, width),
            values.len() as u64,
        )
        .with_shape(&[values.len() as u64]);
    set_h5_chunking(dataset, values.len());
}

fn set_h5_chunking(dataset: &mut DatasetBuilder, len: usize) {
    if len > 0 {
        // Bound the compression working set and avoid three-byte dimensions,
        // which hdf5-pure 0.47 encodes with four bytes and HDF5 2.0 rejects.
        dataset
            .with_chunks(&[len.min(16_384) as u64])
            .with_deflate(1);
    }
}

fn fixed_ascii_type(width: usize) -> Datatype {
    Datatype::String {
        size: width as u32,
        padding: StringPadding::NullTerminate,
        charset: CharacterSet::Ascii,
    }
}

fn fixed_string_bytes(values: &[&str], width: usize) -> Vec<u8> {
    let mut out = vec![0u8; values.len() * width];
    for (row, value) in values.iter().enumerate() {
        let start = row * width;
        let bytes = value.as_bytes();
        out[start..start + bytes.len()].copy_from_slice(bytes);
    }
    out
}

fn u32_slice_to_i32(values: &[u32]) -> Result<Vec<i32>> {
    values
        .iter()
        .map(|&value| {
            i32::try_from(value).map_err(|_| {
                Error::InvalidFormat("value exceeds kallisto-compatible HDF5 range".into())
            })
        })
        .collect()
}

fn h5_error(err: hdf5_pure::Error) -> Error {
    Error::InvalidFormat(format!("HDF5 output failed: {err}"))
}
