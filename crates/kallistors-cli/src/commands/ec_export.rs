use anyhow::{anyhow, Result};
use std::path::Path;

pub fn run(ec_path: &Path, out_path: &Path) -> Result<()> {
    let ec = kallistors::ec::EcList::load_text(ec_path)
        .map_err(|err| anyhow!("ec-export failed: {err}"))?;
    ec.write_text(out_path)
        .map_err(|err| anyhow!("ec-export failed: {err}"))?;

    println!("note: debug-only command (not in kallisto)");
    println!("wrote: {}", out_path.display());
    Ok(())
}
