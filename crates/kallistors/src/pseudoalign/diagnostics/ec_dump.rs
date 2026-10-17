//! Diagnostic EC-block export.

use std::collections::HashSet;
use std::io::Write;
use std::path::Path;

use super::super::build_bifrost_index;
use crate::Result;

pub fn write_index_ec_dump(index_path: &Path, out: &Path, limit: usize) -> Result<()> {
    write_index_ec_dump_filtered(index_path, out, limit, None)
}

pub fn write_index_ec_dump_filtered(
    index_path: &Path,
    out: &Path,
    limit: usize,
    unitig_filter: Option<&HashSet<usize>>,
) -> Result<()> {
    let index = build_bifrost_index(index_path)?;
    let mut writer = std::io::BufWriter::new(std::fs::File::create(out)?);
    writeln!(
        writer,
        "unitig_id\tunitig_pos\tkmer\tblock_id\tblock_lb\tblock_ub\tec"
    )?;
    let mut rows = 0usize;
    let mut remaining = unitig_filter.map(|v| v.len());
    for (unitig_id, blocks) in index.ec_blocks.iter().enumerate() {
        if let Some(filter) = unitig_filter
            && !filter.contains(&unitig_id)
        {
            continue;
        }
        let seq = if unitig_id < index.unitig_count() {
            index.unitig_sequence(unitig_id).expect("valid graph node")
        } else {
            let km_idx = unitig_id - index.unitig_count();
            if km_idx < index.km_unitig_count() {
                index
                    .km_unitig_sequence(km_idx)
                    .expect("valid fixed graph node")
            } else {
                continue;
            }
        };
        for (block_id, block) in blocks.iter().enumerate() {
            let pos = block.lb as usize;
            let kmer = if pos + index.k <= seq.len() {
                String::from_utf8_lossy(&seq[pos..pos + index.k]).into_owned()
            } else if pos < seq.len() {
                String::from_utf8_lossy(&seq[pos..]).into_owned()
            } else {
                "-".to_string()
            };
            let ec_text = if block.ec.is_empty() {
                "-".to_string()
            } else {
                block
                    .ec
                    .iter()
                    .map(|v| v.to_string())
                    .collect::<Vec<_>>()
                    .join(",")
            };
            writeln!(
                writer,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                unitig_id, pos, kmer, block_id, block.lb, block.ub, ec_text
            )?;
            rows += 1;
            if limit > 0 && rows >= limit {
                return Ok(());
            }
        }
        if let Some(ref mut left) = remaining
            && *left > 0
        {
            *left -= 1;
            if *left == 0 {
                return Ok(());
            }
        }
    }
    Ok(())
}
