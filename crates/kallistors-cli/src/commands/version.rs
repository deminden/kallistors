use anyhow::Result;

pub fn run() -> Result<()> {
    println!("kallistors-cli {}", env!("CARGO_PKG_VERSION"));
    Ok(())
}
