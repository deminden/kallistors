# kallistors

A high-performance Rust implementation of kallisto-style pseudoalignment and RNA-seq transcript quantification.

Build kallisto v13-compatible indexes, quantify paired-end or single-end reads, and generate
single-cell BUS output. Quantification writes kallisto-style TSV, JSON and HDF5 results.

```bash
cargo install kallistors --locked
kallistors --help
```

For Rust applications:

```toml
[dependencies]
kallistors = "0.5.1"
```

- [Quick start and performance comparisons](https://github.com/deminden/kallistors/blob/main/README.md)
- [Command-line guide](https://github.com/deminden/kallistors/blob/main/docs/usage.md)
- [Rust API examples](https://github.com/deminden/kallistors/blob/main/docs/development.md#using-the-rust-library)
- [Supported behavior and compatibility limits](https://github.com/deminden/kallistors/blob/main/docs/compatibility.md)
- [Match-cache memory settings](https://github.com/deminden/kallistors/blob/main/docs/match-cache.md)

Licensed under the MIT License; see the included `LICENSE` file.
