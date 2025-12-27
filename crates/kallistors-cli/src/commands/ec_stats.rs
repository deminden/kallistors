use anyhow::{Result, anyhow};
use std::path::Path;

pub fn run(ec_path: &Path) -> Result<()> {
    let ec = kallistors::ec::EcList::load_text(ec_path)
        .map_err(|err| anyhow!("ec-stats failed: {err}"))?;
    let (num_ecs, max_ec_size) = ec.stats();

    println!("note: debug-only command (not in kallisto)");
    println!("ecs: {}", num_ecs);
    println!("max_ec_size: {}", max_ec_size);
    Ok(())
}
