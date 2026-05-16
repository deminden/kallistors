use anyhow::Result;

pub fn run() -> Result<()> {
    println!("kallistors {}", env!("CARGO_PKG_VERSION"));
    Ok(())
}
