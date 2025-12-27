use anyhow::{anyhow, Result};
use std::path::Path;

pub fn run(a_path: &Path, b_path: &Path) -> Result<()> {
    let a = kallistors::ec::EcList::load_text(a_path)
        .map_err(|err| anyhow!("ec-compare failed: {err}"))?;
    let b = kallistors::ec::EcList::load_text(b_path)
        .map_err(|err| anyhow!("ec-compare failed: {err}"))?;

    println!("note: debug-only command (not in kallisto)");
    if a.classes.len() != b.classes.len() {
        println!("ecs: {} vs {}", a.classes.len(), b.classes.len());
        return Ok(());
    }

    for (i, (ac, bc)) in a.classes.iter().zip(b.classes.iter()).enumerate() {
        if ac != bc {
            println!("mismatch at ec {}", i);
            println!("a: {:?}", ac);
            println!("b: {:?}", bc);
            return Ok(());
        }
    }

    println!("match: {} ECs", a.classes.len());
    Ok(())
}
