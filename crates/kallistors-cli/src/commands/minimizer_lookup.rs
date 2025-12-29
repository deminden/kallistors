use std::path::PathBuf;

use anyhow::{Result, anyhow};

pub fn run(index: PathBuf, minimizer: String) -> Result<()> {
    let index = kallistors::pseudoalign::build_bifrost_index(&index)?;
    let bytes = kallistors::index::bifrost::encode_minimizer_rep(minimizer.as_bytes())
        .ok_or_else(|| anyhow!("invalid minimizer sequence"))?;
    let (lookup, debug) = index.mphf.debug_lookup(&bytes);
    let bytes_hex = bytes
        .iter()
        .map(|b| format!("{b:02x}"))
        .collect::<Vec<_>>()
        .join("");
    let positions = lookup
        .and_then(|idx| index.minz_positions.get(idx as usize))
        .map(|v| v.len())
        .unwrap_or(0);

    println!("minimizer\t{minimizer}");
    println!("mphf_key_hex\t{bytes_hex}");
    println!(
        "mphf_lookup\t{}",
        lookup
            .map(|v| v.to_string())
            .unwrap_or_else(|| "-".to_string())
    );
    println!("mphf_level\t{}", debug.level);
    println!("mphf_hash_raw\t{}", debug.hash_raw);
    println!("mphf_hash_domain\t{}", debug.hash_domain);
    println!(
        "mphf_bucket\t{}",
        debug
            .bucket
            .map(|v| v.to_string())
            .unwrap_or_else(|| "-".to_string())
    );
    println!(
        "mphf_rank\t{}",
        debug
            .rank
            .map(|v| v.to_string())
            .unwrap_or_else(|| "-".to_string())
    );
    println!(
        "mphf_final_hash_hit\t{}",
        if debug.final_hash_hit { 1 } else { 0 }
    );
    println!("positions_count\t{positions}");
    Ok(())
}
