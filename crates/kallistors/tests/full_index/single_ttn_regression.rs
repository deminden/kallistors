use std::fs;
use std::path::PathBuf;
use std::process::Command;

#[test]
#[ignore = "requires the Rust-built GENCODE v50 full index"]
fn reversed_unitig_backoff_preserves_complete_svil_ecs_in_both_scanners() {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .parent()
        .unwrap()
        .to_path_buf();
    let index = std::env::var_os("KALLISTORS_TEST_RUST_INDEX")
        .map(PathBuf::from)
        .unwrap_or_else(|| root.join("data/srr-small/gencode.v50.kallistors.idx"));
    assert!(
        index.is_file(),
        "prepare the Rust-built GENCODE v50 index first"
    );
    let binary = std::env::var_os("KALLISTORS_TEST_BIN")
        .unwrap_or_else(|| env!("CARGO_BIN_EXE_kallistors").into());
    let dir = tempfile::tempdir().unwrap();
    let reads = dir.path().join("reads.fq");
    fs::write(
        &reads,
        include_bytes!("../fixtures/single_svil_gencode50.fixture"),
    )
    .unwrap();
    let expected = std::collections::BTreeMap::from([
        ((331277..331313).collect::<Vec<u32>>(), 1u64),
        ((331277..331314).collect::<Vec<u32>>(), 1u64),
    ]);
    for general in [false, true] {
        let out = dir
            .path()
            .join(if general { "general.tsv" } else { "fast.tsv" });
        let mut command = Command::new(&binary);
        command
            .arg("pseudoalign")
            .arg("--index")
            .arg(&index)
            .arg("--reads")
            .arg(&reads)
            .args(["--fragment-length", "200", "--threads", "8"])
            .arg("--out")
            .arg(&out);
        if general {
            command.env("KALLISTORS_DISABLE_FAST_PATH", "1");
        }
        let result = command.output().unwrap();
        assert!(
            result.status.success(),
            "{}",
            String::from_utf8_lossy(&result.stderr)
        );
        let actual: std::collections::BTreeMap<Vec<u32>, u64> = fs::read_to_string(out)
            .unwrap()
            .lines()
            .filter(|line| !line.starts_with('#'))
            .filter_map(|line| {
                let fields: Vec<_> = line.split('\t').collect();
                let count = fields[1].parse().unwrap();
                (count > 0).then(|| {
                    (
                        fields[2].split(',').map(|v| v.parse().unwrap()).collect(),
                        count,
                    )
                })
            })
            .collect();
        assert_eq!(actual, expected, "general scanner={general}");
    }
}

#[test]
#[ignore = "requires the downloaded GENCODE v50 full index"]
fn fallback_probe_preserves_unique_ttn_and_agap_reads_in_both_scanners() {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(|p| p.parent())
        .expect("workspace root")
        .to_path_buf();
    let index = root.join("data/srr-small/gencode.v50.kallisto.idx");
    assert!(index.is_file(), "prepare the GENCODE v50 index first");
    let dir = tempfile::tempdir().expect("tempdir");
    let fixture = dir.path().join("reads.fixture");
    let mut reads = include_bytes!("../fixtures/single_ttn_gencode50.fixture").to_vec();
    reads.extend_from_slice(include_bytes!(
        "../fixtures/single_agap12p_gencode50.fixture"
    ));
    fs::write(&fixture, reads).expect("write regression reads");
    let binary = std::env::var_os("KALLISTORS_TEST_BIN")
        .unwrap_or_else(|| env!("CARGO_BIN_EXE_kallistors").into());
    for general_scanner in [false, true] {
        let out = dir
            .path()
            .join(if general_scanner { "general" } else { "fast" });
        let mut command = Command::new(&binary);
        command
            .arg("quant")
            .args(["--single", "-l", "200", "-s", "20", "-t", "8"])
            .arg("-i")
            .arg(&index)
            .arg("-o")
            .arg(&out)
            .arg(&fixture);
        if general_scanner {
            command.env("KALLISTORS_DISABLE_FAST_PATH", "1");
        }
        let result = command.output().expect("run single-end quant");
        assert!(
            result.status.success(),
            "{}",
            String::from_utf8_lossy(&result.stderr)
        );
        let info = fs::read_to_string(out.join("run_info.json")).expect("run info");
        for counter in ["n_processed", "n_pseudoaligned", "n_unique"] {
            assert!(info.contains(&format!("\"{counter}\": 3,")), "{info}");
        }
        let abundance = fs::read_to_string(out.join("abundance.tsv")).expect("abundance");
        let row = abundance
            .lines()
            .find(|line| line.starts_with("ENST00000425332.4|"))
            .expect("TTN-206 row");
        assert_eq!(
            row.split('\t')
                .nth(3)
                .expect("count")
                .parse::<f64>()
                .expect("numeric count"),
            2.0
        );
        let row = abundance
            .lines()
            .find(|line| line.starts_with("ENST00000805061.1|"))
            .expect("AGAP12P-205 row");
        assert_eq!(
            row.split('\t')
                .nth(3)
                .expect("count")
                .parse::<f64>()
                .expect("numeric count"),
            1.0
        );
    }
}
