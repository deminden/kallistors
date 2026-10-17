use std::fs::{self, File};
use std::io::BufReader;
use std::path::Path;
use std::process::Command;

use kallistors::io::FastqReader;
use kallistors::pseudoalign::{
    FragmentFilter, Strand, build_bifrost_index_with_positions,
    pseudoalign_single_end_bifrost_with_strand_and_filter,
};

use crate::helpers;

#[test]
#[ignore = "requires original kallisto on PATH"]
fn single_fragment_filter_and_overhang_match_kallisto() {
    assert!(
        helpers::kallisto_available(),
        "original kallisto must be on PATH"
    );
    let directory = tempfile::tempdir().expect("parity directory");
    let dataset = helpers::build_upstream_synthetic_index();
    let reads = directory.path().join("reads.fq");
    let mut fastq = Vec::new();
    for (tid, transcript) in dataset.transcripts.iter().enumerate() {
        // End-proximal reads distinguish the fragment filter from --single-overhang.
        for start in [0, 60, transcript.len() - 50] {
            let sequence = transcript[start..start + 50].to_vec();
            for (reverse, read) in [sequence.clone(), helpers::revcomp(&sequence)]
                .iter()
                .enumerate()
            {
                fastq.extend_from_slice(format!("@tx{tid}_{start}_{reverse}\n").as_bytes());
                fastq.extend_from_slice(read);
                fastq.extend_from_slice(b"\n+\n");
                fastq.extend(std::iter::repeat_n(b'I', read.len()));
                fastq.push(b'\n');
            }
        }
    }
    fs::write(&reads, fastq).expect("write reads");
    let index = build_bifrost_index_with_positions(&dataset.index_path, true).expect("load index");
    let mut aligned = Vec::new();
    for overhang in [false, true] {
        let mut reader = FastqReader::new(BufReader::new(File::open(&reads).unwrap()));
        let ours = pseudoalign_single_end_bifrost_with_strand_and_filter(
            &index,
            &mut reader,
            Strand::Unstranded,
            Some(FragmentFilter {
                fragment_length: 100,
                single_overhang: overhang,
            }),
        )
        .expect("single counts");
        let expected = run_kallisto(&dataset.index_path, &reads, overhang);
        assert_eq!(ours.reads_processed, 24);
        assert_eq!(
            ours.reads_aligned, expected,
            "single parity, overhang={overhang}"
        );
        assert!(expected > 0);
        aligned.push(expected);
    }
    assert!(
        aligned[1] > aligned[0],
        "fixture must exercise the overhang distinction"
    );
}

fn run_kallisto(index: &Path, reads: &Path, overhang: bool) -> u64 {
    let directory = tempfile::tempdir().expect("native quant directory");
    let out = directory.path().join("quant");
    let mut command = Command::new("kallisto");
    command
        .args(["quant", "--single", "-l", "100", "-s", "20", "-i"])
        .arg(index)
        .arg("-o")
        .arg(&out)
        .arg(reads);
    if overhang {
        command.arg("--single-overhang");
    }
    let output = command.output().expect("run original kallisto");
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let info = fs::read_to_string(out.join("run_info.json")).expect("read native metadata");
    assert_eq!(helpers::json_u64(&info, "n_processed"), Some(24));
    helpers::json_u64(&info, "n_pseudoaligned").expect("native aligned count")
}
