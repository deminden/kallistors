use std::assert_matches;
use std::collections::HashSet;
use std::fs;

use kallistors::Error;
use kallistors::index::{
    Index, IndexBuildOptions, bifrost::encode_minimizer_rep, build_index, build_index_with_report,
    extract_ec_list,
};
use kallistors::pseudoalign::build_bifrost_index_with_kmer;

#[test]
fn builds_loadable_v13_index_from_fasta() {
    let dir = tempfile::tempdir().expect("tempdir");
    let fasta = dir.path().join("transcripts.fa");
    let index = dir.path().join("transcripts.idx");
    fs::write(
        &fasta,
        b">tx1
ACGTGCACTGATCGTACGATCGTACGTTAGCTAGCTAGGCTAGCATCGATCGATGCTAGCTAGCTGACT
>tx2
TTGACCGTAGCTAGGATCCGATCGTACGATCGTAGCTAGCTAACGTTAGCTAGGCTACGATCGATCGT
",
    )
    .expect("write fasta");

    let report = build_index_with_report(
        &index,
        &[fasta],
        IndexBuildOptions {
            k: 31,
            g: None,
            threads: 1,
            make_unique: false,
            ec_max_size: -1,
            aa: false,
        },
    )
    .expect("build index");
    assert_eq!(report.transcripts, 2);
    assert!(report.unitigs < report.kmers);
    assert!(report.minimizers > 0);

    let meta = Index::load(&index).expect("load metadata");
    assert_eq!(meta.index_version, 13);
    assert_eq!(meta.k, 31);
    assert_eq!(meta.minimizer_len, Some(23));
    assert!(meta.unitigs.expect("unitigs") < meta.kmers.expect("kmers"));
    assert_eq!(meta.transcripts.len(), 2);
    assert_eq!(meta.transcripts[0].name, "tx1");
    assert_eq!(meta.transcripts[1].name, "tx2");

    let bifrost = build_bifrost_index_with_kmer(&index).expect("load bifrost index");
    assert_eq!(bifrost.k, 31);
    assert_eq!(bifrost.g, 23);
    assert_eq!(bifrost.transcript_names, ["tx1", "tx2"]);
    assert!(
        bifrost
            .kmer_index
            .as_ref()
            .is_some_and(|idx| !idx.ecs.is_empty())
    );

    let ecs = extract_ec_list(&index).expect("extract ecs");
    assert!(ecs.classes.iter().any(|ec| ec == &[0]));
    assert!(ecs.classes.iter().any(|ec| ec == &[1]));
}

#[test]
fn generated_index_uses_compact_mphf_cascade() {
    let dir = tempfile::tempdir().expect("tempdir");
    let fasta = dir.path().join("transcripts.fa");
    let index = dir.path().join("transcripts.idx");
    fs::write(
        &fasta,
        b">tx1
ACGTGCACTGATCGTACGATCGTACGTTAGCTAGCTAGGCTAGCATCGATCGATGCTAGCTAGCTGACT
>tx2
TTGACCGTAGCTAGGATCCGATCGTACGATCGTAGCTAGCTAACGTTAGCTAGGCTACGATCGATCGT
",
    )
    .expect("write fasta");

    build_index(
        &index,
        &[fasta],
        IndexBuildOptions {
            k: 31,
            ..IndexBuildOptions::default()
        },
    )
    .expect("build index");

    let bifrost = build_bifrost_index_with_kmer(&index).expect("load bifrost index");
    let mut keys = HashSet::new();
    for unitig in &bifrost.unitigs {
        for window in unitig.windows(bifrost.g) {
            keys.insert(encode_minimizer_rep(window).expect("valid minimizer"));
        }
    }

    assert!(keys.len() > 10);
    assert!(bifrost.mphf.size() > 0);
    assert!(
        bifrost.mphf.size() < keys.len() as u64,
        "builder should index selected minimizers, not every g-mer window"
    );
    assert!(bifrost.mphf.final_hash_entries().len() < bifrost.mphf.size() as usize);

    for rank in 0..bifrost.mphf.size() {
        assert!(!bifrost.minz_positions.get(rank as usize).is_empty());
    }
}

#[test]
fn duplicate_transcript_names_require_make_unique() {
    let dir = tempfile::tempdir().expect("tempdir");
    let fasta = dir.path().join("dupes.fa");
    fs::write(
        &fasta,
        b">tx
ACGTGCACTGATCGTACGATCGTACGTTAGCTAGCTAGGCTAGCATCGATCGATGCTAGCTAGCTGACT
>tx
TTGACCGTAGCTAGGATCCGATCGTACGATCGTAGCTAGCTAACGTTAGCTAGGCTACGATCGATCGT
",
    )
    .expect("write fasta");

    let err = build_index(
        &dir.path().join("dupes_error.idx"),
        std::slice::from_ref(&fasta),
        IndexBuildOptions {
            k: 31,
            ..IndexBuildOptions::default()
        },
    )
    .expect_err("duplicate names should fail");
    assert_matches!(err, Error::InvalidFormat(message) if message.contains("repeated name"));

    let index = dir.path().join("dupes_unique.idx");
    build_index(
        &index,
        &[fasta],
        IndexBuildOptions {
            k: 31,
            make_unique: true,
            ..IndexBuildOptions::default()
        },
    )
    .expect("build with unique names");
    let meta = Index::load(&index).expect("load metadata");
    assert_eq!(meta.transcripts[0].name, "tx");
    assert_eq!(meta.transcripts[1].name, "tx_1");
}

#[test]
fn builds_index_for_empty_and_short_transcripts_without_minimizers() {
    let dir = tempfile::tempdir().expect("tempdir");
    let fasta = dir.path().join("short.fa");
    let index = dir.path().join("short.idx");
    fs::write(
        &fasta,
        b">empty
>short
ACGTACGTACGT
>poly_a_only
AAAAAAAAAAAA
",
    )
    .expect("write fasta");

    let report = build_index_with_report(
        &index,
        &[fasta],
        IndexBuildOptions {
            k: 31,
            ..IndexBuildOptions::default()
        },
    )
    .expect("build short index");

    assert_eq!(report.transcripts, 3);
    assert_eq!(report.unitigs, 0);
    assert_eq!(report.kmers, 0);
    assert_eq!(report.minimizers, 0);

    let meta = Index::load(&index).expect("load metadata");
    assert_eq!(meta.index_version, 13);
    assert_eq!(meta.transcripts.len(), 3);
    assert_eq!(meta.transcripts[0].length, 0);
    assert_eq!(meta.transcripts[1].length, 12);
    assert_eq!(meta.transcripts[2].length, 12);

    let bifrost = build_bifrost_index_with_kmer(&index).expect("load bifrost index");
    assert!(bifrost.unitigs.is_empty());
    assert_eq!(bifrost.mphf.size(), 0);
}

#[test]
fn normalizes_ambiguous_bases_and_preserves_original_lengths() {
    let dir = tempfile::tempdir().expect("tempdir");
    let fasta = dir.path().join("ambiguous.fa");
    let index = dir.path().join("ambiguous.idx");
    fs::write(
        &fasta,
        b">rna_and_masked
acgtuNNNNrrrryyyyACGTACGTACGTACGTACGTACGTACGTACGTAAAAAAAAAA
>mixed_case
ttttccccaaaaggggNNNNACGTACGTACGTACGTACGTACGTACGTACGT
",
    )
    .expect("write fasta");

    let report = build_index_with_report(
        &index,
        &[fasta],
        IndexBuildOptions {
            k: 31,
            ..IndexBuildOptions::default()
        },
    )
    .expect("build ambiguous index");

    assert_eq!(report.transcripts, 2);
    assert!(report.kmers > 0);
    assert!(report.minimizers > 0);

    let meta = Index::load(&index).expect("load metadata");
    assert_eq!(meta.transcripts[0].length, 59);
    assert_eq!(meta.transcripts[1].length, 52);
    assert_eq!(meta.transcripts[0].name, "rna_and_masked");
    assert_eq!(meta.transcripts[1].name, "mixed_case");

    let bifrost = build_bifrost_index_with_kmer(&index).expect("load bifrost index");
    assert!(!bifrost.unitigs.is_empty());
    assert_eq!(bifrost.transcript_names, ["rna_and_masked", "mixed_case"]);
}

#[test]
fn repetitive_duplicate_transcripts_are_deduplicated_into_loadable_graph() {
    let dir = tempfile::tempdir().expect("tempdir");
    let fasta = dir.path().join("repetitive.fa");
    let index = dir.path().join("repetitive.idx");
    let mut fasta_text = String::new();
    for i in 0..24 {
        fasta_text.push_str(&format!(">dup\n{}\n", "ACGT".repeat(32 + (i % 3))));
    }
    fasta_text.push_str(">homopolymer\n");
    fasta_text.push_str(&"A".repeat(96));
    fasta_text.push('\n');
    fs::write(&fasta, fasta_text).expect("write fasta");

    let report = build_index_with_report(
        &index,
        &[fasta],
        IndexBuildOptions {
            k: 31,
            make_unique: true,
            threads: 4,
            ..IndexBuildOptions::default()
        },
    )
    .expect("build repetitive index");

    assert_eq!(report.transcripts, 25);
    assert!(report.kmers > 0);
    assert!(report.unitigs > 0);

    let meta = Index::load(&index).expect("load metadata");
    assert_eq!(meta.transcripts[0].name, "dup");
    assert_eq!(meta.transcripts[1].name, "dup_1");
    assert_eq!(meta.transcripts[23].name, "dup_23");
    assert_eq!(meta.transcripts[24].name, "homopolymer");

    let bifrost = build_bifrost_index_with_kmer(&index).expect("load bifrost index");
    assert!(!bifrost.unitigs.is_empty());
    assert!(bifrost.mphf.size() > 0);
}
