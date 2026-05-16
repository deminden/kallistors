use anyhow::{Result, anyhow};
use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

pub fn run(index_path: &Path, out: &Path, limit: usize, unitig_list: Option<&Path>) -> Result<()> {
    let unitig_filter = if let Some(path) = unitig_list {
        let file = File::open(path)?;
        let mut set = HashSet::new();
        for line in BufReader::new(file).lines() {
            let line = line?;
            let text = line.trim();
            if text.is_empty() {
                continue;
            }
            let id: usize = text
                .parse()
                .map_err(|_| anyhow!("invalid unitig id in list: {text}"))?;
            set.insert(id);
        }
        Some(set)
    } else {
        None
    };
    Ok(kallistors::pseudoalign::write_index_ec_dump_filtered(
        index_path,
        out,
        limit,
        unitig_filter.as_ref(),
    )?)
}
