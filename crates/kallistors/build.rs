fn main() {
    println!("cargo::rerun-if-changed=build.rs");
    if std::env::var("CARGO_CFG_TARGET_OS").as_deref() == Ok("windows")
        && std::env::var("CARGO_CFG_TARGET_ENV").as_deref() == Ok("msvc")
    {
        // The unoptimized clap subcommand builder exceeds MSVC's 1 MiB default.
        // Reserve room for CLI startup without changing the initial committed size,
        // worker stacks, or the library's callers.
        println!("cargo::rustc-link-arg-bin=kallistors=/STACK:8388608");
    }
}
