use crate::helpers::revcomp;
use std::io::Cursor;

use kallistors::io::FastqReader;
use kallistors::pseudoalign::{
    Strand, build_bifrost_index, pseudoalign_single_end_bifrost_with_strand,
};

use crate::helpers;

#[test]
fn reverse_strand_aligns_only_with_revcomp() {
    let dataset = helpers::build_synthetic_index();
    let index = build_bifrost_index(&dataset.index_path).expect("index");
    assert_eq!(index.transcript_names.len(), dataset.transcripts.len());
    // Strand restricts lookup relative to stored unitigs. Their orientation is
    // chosen by graph construction and need not match the input transcript.
    let unitigs = index
        .unitigs
        .iter()
        .chain(&index.km_unitigs)
        .map(Vec::as_slice)
        .collect::<Vec<_>>();
    let read_len = index.k + 10;
    let mut read = None;
    'outer: for seq in &unitigs {
        if seq.len() < read_len {
            continue;
        }
        for start in 0..=seq.len() - read_len {
            let slice = &seq[start..start + read_len];
            let rc = revcomp(slice);
            if rc
                .windows(index.k)
                .all(|kmer| !unitigs.iter().any(|unitig| contains_subseq(unitig, kmer)))
            {
                read = Some(slice.to_vec());
                break 'outer;
            }
        }
    }
    let read = read.expect("find non-palindromic read");
    let rc = revcomp(&read);
    let mut fastq = Vec::new();
    fastq.extend_from_slice(b"@rc\n");
    fastq.extend_from_slice(&rc);
    fastq.extend_from_slice(b"\n+\n");
    fastq.extend_from_slice(&vec![b'I'; rc.len()]);
    fastq.extend_from_slice(b"\n");

    let mut rc_reader = FastqReader::new(Cursor::new(fastq.clone()));
    let res_forward =
        pseudoalign_single_end_bifrost_with_strand(&index, &mut rc_reader, Strand::Forward)
            .expect("forward");
    assert_eq!(res_forward.reads_aligned, 0);

    let mut rc_reader = FastqReader::new(Cursor::new(fastq.clone()));
    let res_reverse =
        pseudoalign_single_end_bifrost_with_strand(&index, &mut rc_reader, Strand::Reverse)
            .expect("reverse");
    assert_eq!(res_reverse.reads_aligned, 1);

    let mut rc_reader = FastqReader::new(Cursor::new(fastq));
    let res_unstranded =
        pseudoalign_single_end_bifrost_with_strand(&index, &mut rc_reader, Strand::Unstranded)
            .expect("unstranded");
    assert_eq!(res_unstranded.reads_aligned, 1);
}

fn contains_subseq(haystack: &[u8], needle: &[u8]) -> bool {
    if needle.is_empty() {
        return true;
    }
    haystack
        .windows(needle.len())
        .any(|window| window == needle)
}

#[test]
fn paired_strand_flags_select_expected_transcripts_in_both_cli_commands() {
    use std::collections::BTreeMap;
    use std::fs;
    use std::process::Command;

    let dataset = helpers::build_synthetic_index();
    let dir = tempfile::tempdir().unwrap();
    let reads1 = dir.path().join("left.fq");
    let reads2 = dir.path().join("right.fq");
    let forward = &dataset.transcripts[2];
    let reverse = &dataset.transcripts[3];
    let (left, right) = helpers::to_fastq_pairs(&[
        (forward[..50].to_vec(), revcomp(&forward[40..90])),
        (revcomp(&reverse[40..90]), reverse[..50].to_vec()),
        // A mapped mate rejected by strand must reject the pair, while an
        // actually unmapped mate still permits the correctly oriented mate.
        (forward[..50].to_vec(), forward[40..90].to_vec()),
        (forward[..50].to_vec(), vec![b'N'; 50]),
        (vec![b'N'; 50], reverse[..50].to_vec()),
    ]);
    fs::write(&reads1, left).unwrap();
    fs::write(&reads2, right).unwrap();
    for threads in [1, 8] {
        for (label, flag, targets) in [
            ("unstranded", None, vec![2, 2, 2, 3, 3]),
            ("fr", Some("--fr-stranded"), vec![2, 2]),
            ("rf", Some("--rf-stranded"), vec![3, 3]),
        ] {
            let mut expected = BTreeMap::new();
            for &id in &targets {
                *expected.entry(vec![id]).or_insert(0u64) += 1;
            }
            let ec_path = dir.path().join(format!("{label}-{threads}.tsv"));
            let mut pseudo = Command::new(env!("CARGO_BIN_EXE_kallistors"));
            pseudo
                .arg("pseudoalign")
                .arg("--index")
                .arg(&dataset.index_path)
                .arg("--reads")
                .arg(&reads1)
                .arg("--reads2")
                .arg(&reads2)
                .arg("--out")
                .arg(&ec_path)
                .arg("--threads")
                .arg(threads.to_string());
            if let Some(flag) = flag {
                pseudo.arg(flag);
            }
            let result = pseudo.output().unwrap();
            assert!(
                result.status.success(),
                "{label}/{threads}: {}",
                String::from_utf8_lossy(&result.stderr)
            );
            let actual: BTreeMap<Vec<u32>, u64> = fs::read_to_string(ec_path)
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
            assert_eq!(actual, expected, "pseudoalign {label}/{threads}");

            let out = dir.path().join(format!("quant-{label}-{threads}"));
            let mut quant = Command::new(env!("CARGO_BIN_EXE_kallistors"));
            quant
                .arg("quant")
                .arg("-i")
                .arg(&dataset.index_path)
                .arg("-o")
                .arg(&out)
                .arg("-t")
                .arg(threads.to_string())
                .arg(&reads1)
                .arg(&reads2);
            if let Some(flag) = flag {
                quant.arg(flag);
            }
            let result = quant.output().unwrap();
            assert!(
                result.status.success(),
                "{label}/{threads}: {}",
                String::from_utf8_lossy(&result.stderr)
            );
            let abundance = helpers::read_abundance(out.join("abundance.tsv"));
            for id in 0..4 {
                assert_eq!(
                    abundance[&format!("tx{id}")].2,
                    targets.iter().filter(|&&target| target == id).count() as f64,
                    "quant {label}/{threads}/tx{id}"
                );
            }
        }
    }
}

#[test]
fn strand_diagnostics_load_annotations_without_a_fragment_filter() {
    use std::fs;
    use std::process::Command;

    let dataset = helpers::build_synthetic_index();
    let dir = tempfile::tempdir().unwrap();
    let reads = dir.path().join("reads.fq");
    let names = dir.path().join("names.txt");
    fs::write(
        &reads,
        helpers::to_fastq(&[
            dataset.transcripts[2][..50].to_vec(),
            revcomp(&dataset.transcripts[3][..50]),
        ]),
    )
    .unwrap();
    fs::write(&names, "r1\nr2\n").unwrap();
    for command in ["trace-reads", "trace-compare"] {
        for flag in ["--fr-stranded", "--rf-stranded"] {
            for loader in [
                None,
                Some("--kallisto-direct-kmer"),
                Some("--kallisto-fallback"),
            ] {
                let loader_name = loader.unwrap_or("default");
                let case = format!("{command}-{flag}-{loader_name}");
                let out = dir.path().join(&case);
                let mut cmd = Command::new(env!("CARGO_BIN_EXE_kallistors"));
                cmd.arg(command)
                    .arg("--index")
                    .arg(&dataset.index_path)
                    .arg("--reads")
                    .arg(&reads)
                    .arg("--read-list")
                    .arg(&names)
                    .arg("--out")
                    .arg(&out)
                    .arg(flag);
                if let Some(loader) = loader {
                    cmd.arg(loader);
                }
                let result = cmd.output().unwrap();
                assert!(
                    result.status.success(),
                    "{case}: {}",
                    String::from_utf8_lossy(&result.stderr)
                );
                let table = fs::read_to_string(&out).unwrap();
                let mut lines = table.lines();
                let header: Vec<_> = lines.next().unwrap().split('\t').collect();
                let columns: &[&str] = if command == "trace-reads" {
                    &["ec_after_strand"]
                } else {
                    &["ec_baseline", "ec_fast"]
                };
                let rows: Vec<Vec<_>> = lines.map(|row| row.split('\t').collect()).collect();
                assert_eq!(rows.len(), 2, "{case}");
                let expected = if flag == "--fr-stranded" {
                    ["2", "-"]
                } else {
                    ["-", "3"]
                };
                for column in columns {
                    let col = header.iter().position(|v| v == column).unwrap();
                    for (row, expected) in rows.iter().zip(expected) {
                        assert_eq!(row[col], expected, "{case}/{column}");
                    }
                }
            }
        }
    }
}
