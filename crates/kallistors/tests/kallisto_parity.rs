use std::fs;
use std::path::Path;
use std::process::Command;

use kallistors::pseudoalign::{
    build_bifrost_index_with_positions, pseudoalign_single_end_bifrost_with_strand_and_filter,
    FragmentFilter, Strand,
};

fn find_n_pseudoaligned(run_info: &str) -> Option<u64> {
    let key = "\"n_pseudoaligned\"";
    let idx = run_info.find(key)?;
    let after = &run_info[idx + key.len()..];
    let colon = after.find(':')?;
    let rest = after[colon + 1..].trim_start();
    let mut digits = String::new();
    for c in rest.chars() {
        if c.is_ascii_digit() {
            digits.push(c);
        } else if !digits.is_empty() {
            break;
        }
    }
    digits.parse().ok()
}

fn kallisto_available() -> bool {
    Command::new("kallisto")
        .arg("version")
        .status()
        .map(|s| s.success())
        .unwrap_or(false)
}

fn run_kallisto_quant(
    index: &Path,
    reads: &Path,
    out_dir: &Path,
    single_overhang: bool,
) -> Option<u64> {
    let mut cmd = Command::new("kallisto");
    cmd.arg("quant")
        .arg("--single")
        .arg("-l")
        .arg("100")
        .arg("-s")
        .arg("20")
        .arg("-i")
        .arg(index)
        .arg("-o")
        .arg(out_dir)
        .arg(reads);
    if single_overhang {
        cmd.arg("--single-overhang");
    }
    let status = cmd.status().ok()?;
    if !status.success() {
        return None;
    }
    let run_info = fs::read_to_string(out_dir.join("run_info.json")).ok()?;
    find_n_pseudoaligned(&run_info)
}

#[test]
fn parity_with_kallisto_single_modes() {
    let index = Path::new("data/large.idx");
    let reads = Path::new("data/diverse_reads.fq");
    if !index.exists() || !reads.exists() {
        return;
    }
    if !kallisto_available() {
        return;
    }

    let tmp_base =
        std::env::temp_dir().join(format!("kallistors_kallisto_parity_{}", std::process::id()));
    let _ = fs::remove_dir_all(&tmp_base);
    fs::create_dir_all(&tmp_base).expect("failed to create temp dir");

    let index_pos = build_bifrost_index_with_positions(index, true)
        .expect("failed to load index with positions");
    let reads_file = fs::File::open(reads).expect("failed to open reads");
    let mut reader = kallistors::io::FastqReader::new(std::io::BufReader::new(reads_file));

    let res = pseudoalign_single_end_bifrost_with_strand_and_filter(
        &index_pos,
        &mut reader,
        Strand::Unstranded,
        Some(FragmentFilter {
            fragment_length: 100,
            single_overhang: false,
        }),
    )
    .expect("pseudoalign failed");
    let expected = run_kallisto_quant(index, reads, &tmp_base.join("no_overhang"), false)
        .expect("kallisto run failed");
    assert_eq!(
        res.reads_aligned, expected,
        "pseudoaligned reads mismatch for --single"
    );

    let reads_file = fs::File::open(reads).expect("failed to open reads");
    let mut reader = kallistors::io::FastqReader::new(std::io::BufReader::new(reads_file));
    let res = pseudoalign_single_end_bifrost_with_strand_and_filter(
        &index_pos,
        &mut reader,
        Strand::Unstranded,
        Some(FragmentFilter {
            fragment_length: 100,
            single_overhang: true,
        }),
    )
    .expect("pseudoalign failed");
    let expected = run_kallisto_quant(index, reads, &tmp_base.join("single_overhang"), true)
        .expect("kallisto run failed");
    assert_eq!(
        res.reads_aligned, expected,
        "pseudoaligned reads mismatch for --single-overhang"
    );
}
