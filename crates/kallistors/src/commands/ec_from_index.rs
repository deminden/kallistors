use anyhow::{Result, anyhow};
use std::path::Path;

pub fn run(index_path: &Path, out_path: &Path) -> Result<()> {
    let ec = kallistors::index::extract_ec_list(index_path)
        .map_err(|err| anyhow!("ec-from-index failed: {err}"))?;
    ec.write_text(out_path)
        .map_err(|err| anyhow!("ec-from-index failed: {err}"))?;

    println!("note: debug-only command (not in kallisto)");
    println!("ecs: {}", ec.classes.len());
    println!("wrote: {}", out_path.display());
    Ok(())
}
