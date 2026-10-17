use crate::helpers::{
    json_f64, json_u64, read_abundance, synthetic_transcripts, write_fasta, write_fastq_pairs,
};
use std::fs;
use std::process::Command;

#[test]
#[ignore = "requires original kallisto on PATH"]
fn paired_quantification_agrees_across_both_index_builders_and_quantifiers() {
    let dir = tempfile::tempdir().unwrap();
    let fasta = dir.path().join("transcripts.fa");
    let upstream_index = dir.path().join("upstream.idx");
    let rust_index = dir.path().join("rust.idx");
    let reads1 = dir.path().join("reads_1.fq");
    let reads2 = dir.path().join("reads_2.fq");
    let transcripts = synthetic_transcripts();
    write_fasta(&fasta, &transcripts).unwrap();
    write_fastq_pairs(&reads1, &reads2, &transcripts, 50, 90, 64).unwrap();
    assert!(
        Command::new("kallisto")
            .arg("index")
            .arg("-i")
            .arg(&upstream_index)
            .arg(&fasta)
            .status()
            .expect("original kallisto on PATH")
            .success()
    );
    kallistors::index::build_index(&rust_index, &[fasta], Default::default()).unwrap();

    let mut baseline: Option<(String, _)> = None;
    for (builder, index) in [("original", &upstream_index), ("rust", &rust_index)] {
        for (quantifier, executable) in [
            ("original", "kallisto"),
            ("rust", env!("CARGO_BIN_EXE_kallistors")),
        ] {
            let label = format!("{builder}-index-{quantifier}-quant");
            let out = dir.path().join(&label);
            let output = Command::new(executable)
                .arg("quant")
                .arg("-i")
                .arg(index)
                .arg("-o")
                .arg(&out)
                .arg(&reads1)
                .arg(&reads2)
                .output()
                .unwrap();
            assert!(
                output.status.success(),
                "{label}: {}",
                String::from_utf8_lossy(&output.stderr)
            );
            let info = fs::read_to_string(out.join("run_info.json")).unwrap();
            let abundance = read_abundance(out.join("abundance.tsv"));
            assert_eq!(json_u64(&info, "n_processed"), Some(64));
            assert_eq!(json_u64(&info, "n_pseudoaligned"), Some(64));
            assert!(
                json_u64(&info, "n_unique").unwrap() < 64,
                "fixture must include ambiguous reads"
            );
            if let Some((expected_info, expected_abundance)) = &baseline {
                for field in [
                    "n_targets",
                    "n_processed",
                    "n_pseudoaligned",
                    "n_unique",
                    "index_version",
                    "k-mer length",
                ] {
                    assert_eq!(
                        json_u64(&info, field),
                        json_u64(expected_info, field),
                        "{label}: {field}"
                    );
                }
                for field in ["p_pseudoaligned", "p_unique"] {
                    assert_eq!(
                        json_f64(&info, field),
                        json_f64(expected_info, field),
                        "{label}: {field}"
                    );
                }
                let expected_abundance: &std::collections::HashMap<String, (u32, f64, f64, f64)> =
                    expected_abundance;
                assert_eq!(abundance.len(), expected_abundance.len());
                for (target, &(length, effective, count, tpm)) in expected_abundance {
                    let row = abundance[target];
                    assert_eq!(row.0, length, "{label}: {target} length");
                    for (field, actual, expected) in [
                        ("effective length", row.1, effective),
                        ("count", row.2, count),
                        ("TPM", row.3, tpm),
                    ] {
                        // Original TSV prints six significant digits. Allow its rounding,
                        // rather than an absolute tolerance that loosens small counts.
                        let tolerance = 1e-10 + 5e-6 * expected.abs();
                        assert!(
                            (actual - expected).abs() <= tolerance,
                            "{label}: {target} {field}: {actual} != {expected}"
                        );
                    }
                }
            } else {
                baseline = Some((info, abundance));
            }
        }
    }
}
