use std::{collections::BTreeMap, fs, path::PathBuf, process::Command};

#[test]
#[ignore = "requires the downloaded GENCODE v50 full index"]
fn abundant_node_keys_preserve_paired_ec_counts_in_both_loaders() {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(|p| p.parent())
        .unwrap()
        .to_path_buf();
    let index = root.join("data/srr-small/gencode.v50.kallisto.idx");
    assert!(index.is_file(), "prepare the GENCODE v50 index first");
    let dir = tempfile::tempdir().unwrap();
    let mates = [
        dir.path().join("left.fastq"),
        dir.path().join("right.fastq"),
    ];
    let mut fastqs = [String::new(), String::new()];
    let mut expected = BTreeMap::<String, u64>::new();
    for line in include_str!("../fixtures/abundant_node_gencode50.fixture").lines() {
        let fields: Vec<_> = line.split('\t').collect();
        for mate in 0..2 {
            let seq = fields[mate + 1];
            fastqs[mate].push_str(&format!(
                "@{}\n{seq}\n+\n{}\n",
                fields[0],
                "I".repeat(seq.len())
            ));
        }
        *expected.entry(fields[3].into()).or_default() += 1;
    }
    for (path, contents) in mates.iter().zip(fastqs) {
        fs::write(path, contents).unwrap();
    }
    let binary = std::env::var_os("KALLISTORS_TEST_BIN")
        .unwrap_or_else(|| env!("CARGO_BIN_EXE_kallistors").into());
    for (threads, debug) in [(1, false), (8, false), (1, true), (8, true)] {
        let output = dir.path().join(format!("ec-{threads}-{debug}.tsv"));
        let mut command = Command::new(&binary);
        command
            .args(["pseudoalign", "--threads", &threads.to_string(), "--index"])
            .arg(&index)
            .arg("--reads")
            .arg(&mates[0])
            .arg("--reads2")
            .arg(&mates[1])
            .arg("--out")
            .arg(&output);
        if debug {
            command
                .arg("--debug-out")
                .arg(dir.path().join(format!("debug-{threads}.tsv")));
        }
        let result = command.output().unwrap();
        assert!(
            result.status.success(),
            "{}",
            String::from_utf8_lossy(&result.stderr)
        );
        let actual: BTreeMap<String, u64> = fs::read_to_string(output)
            .unwrap()
            .lines()
            .filter(|line| !line.starts_with('#'))
            .filter_map(|line| {
                let fields: Vec<_> = line.split('\t').collect();
                let count: u64 = fields[1].parse().unwrap();
                (count != 0).then(|| (fields[2].into(), count))
            })
            .collect();
        assert_eq!(actual, expected);
    }
    // Ordinary pseudoalign and quant use the compact EC loader. Diagnostics
    // above force the general loader; quant below checks both read scanners.
    for general in [false, true] {
        let output = dir.path().join(format!("quant-{general}"));
        let mut command = Command::new(&binary);
        command
            .args(["quant", "-t", "8", "-i"])
            .arg(&index)
            .arg("-o")
            .arg(&output)
            .args(&mates);
        if general {
            command.env("KALLISTORS_DISABLE_FAST_PATH", "1");
        }
        let result = command.output().unwrap();
        assert!(
            result.status.success(),
            "{}",
            String::from_utf8_lossy(&result.stderr)
        );
        let info = fs::read_to_string(output.join("run_info.json")).unwrap();
        assert!(info.contains("\"n_processed\": 7,"), "{info}");
        assert!(info.contains("\"n_pseudoaligned\": 7,"), "{info}");
    }
}
