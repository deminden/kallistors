use anyhow::{Result, anyhow};
use std::path::Path;

pub fn run(index_path: &Path) -> Result<()> {
    let index = kallistors::index::Index::load(index_path)
        .map_err(|err| anyhow!("index-info failed: {err}"))?;

    let total_len: u64 = index.transcripts.iter().map(|t| t.length as u64).sum();

    println!("index_version: {}", index.index_version);
    println!("k: {}", index.k);
    if let Some(g) = index.minimizer_len {
        println!("minimizer_len: {}", g);
    }
    if let Some(unitigs) = index.unitigs {
        println!("unitigs: {}", unitigs);
    }
    if let Some(kmers) = index.kmers {
        println!("kmers: {}", kmers);
    }
    println!("transcripts: {}", index.transcripts.len());
    println!("total_length: {}", total_len);
    Ok(())
}
