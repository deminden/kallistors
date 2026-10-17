//! BUS command configuration and sidecar-output regressions.

use super::{
    BusArgs, bus_pseudoalign_options, bus_query_quality, concat_sequence_slices, concat_slices,
    mark_technology_paired, onlist_target_count, read_groups, string_to_binary, technology_spec,
    total_len, write_batch_long_flens, write_matrix_ec, write_transcripts,
};
use kallistors::io::FastqRecord;
use kallistors::pseudoalign::StrandSpecific;

#[test]
fn bus_binary_encoding_matches_kallisto_order() {
    let mut flag = 99;
    assert_eq!(
        string_to_binary(b"ACGT", &mut flag, "test").unwrap(),
        0b00_01_10_11
    );
    assert_eq!(flag, 0);
}

#[test]
fn preset_metadata_matches_kallisto_layouts() {
    for (name, files, barcode, umi, paired) in [
        ("10xv3", 2, Some(16), Some(12), false),
        ("PARSEV3", 2, Some(24), Some(10), false),
        ("SMARTSEQ3", 4, None, Some(19), true),
    ] {
        let spec = technology_spec(name).unwrap();
        assert_eq!(spec.nfiles, files, "{name}");
        assert_eq!(total_len(&spec.bc), barcode, "{name}");
        assert_eq!(total_len(&spec.umi), umi, "{name}");
        assert_eq!(spec.paired, paired, "{name}");
        assert_eq!(spec.default_strand, Some(StrandSpecific::FR), "{name}");
    }
}

#[test]
fn bulk_is_single_end_by_default_and_pairable() {
    let mut spec = technology_spec("Bulk").unwrap();
    assert_eq!(spec.nfiles, 1);
    assert_eq!(spec.seq.len(), 1);
    assert_eq!(total_len(&spec.bc), None);
    assert_eq!(total_len(&spec.umi), None);
    assert!(!spec.paired);

    mark_technology_paired(&mut spec, "--paired").unwrap();

    assert_eq!(spec.nfiles, 2);
    assert_eq!(spec.seq.len(), 2);
    assert!(spec.paired);
}

#[test]
fn read_groups_rejects_zero_input_file_technologies() {
    let args = BusArgs {
        reads: vec!["reads.fastq".into()],
        ..BusArgs::default()
    };
    let err = match read_groups(&args, 0) {
        Ok(_) => panic!("zero-input-file technology should fail"),
        Err(err) => err.to_string(),
    };
    assert!(err.contains("at least one input file"), "{err}");
}

#[test]
fn bus_dfk_onlist_respects_explicit_and_amino_acid_modes() {
    for (aa, dfk_onlist, expected) in [
        (false, false, false),
        (false, true, true),
        (true, false, true),
    ] {
        let args = BusArgs {
            aa,
            dfk_onlist,
            ..BusArgs::default()
        };
        assert_eq!(
            bus_pseudoalign_options(&args, None).dfk_onlist,
            expected,
            "aa={aa}, explicit={dfk_onlist}"
        );
    }
}

#[test]
fn compound_preset_slices_match_kallisto_layouts() {
    let records = vec![
        FastqRecord {
            header: b"@r1".to_vec(),
            seq: b"0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz".to_vec(),
            plus: b"+".to_vec(),
            qual: vec![b'I'; 62],
        },
        FastqRecord {
            header: b"@r2".to_vec(),
            seq: b"abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789".to_vec(),
            plus: b"+".to_vec(),
            qual: vec![b'I'; 62],
        },
        FastqRecord {
            header: b"@r3".to_vec(),
            seq: b"QRST".to_vec(),
            plus: b"+".to_vec(),
            qual: vec![b'I'; 4],
        },
    ];

    let parse = technology_spec("PARSEV3").unwrap();
    assert_eq!(
        concat_slices(&records, &parse.bc, false).unwrap(),
        b"klmnopqrEFGHIJKLYZ012345"
    );
    assert_eq!(
        concat_slices(&records, &parse.umi, false).unwrap(),
        b"abcdefghij"
    );
    assert_eq!(
        concat_sequence_slices(&records, &parse, false).unwrap(),
        records[0].seq
    );

    let petri = technology_spec("PETRISEQ").unwrap();
    assert_eq!(
        concat_slices(&records, &petri.bc, false).unwrap(),
        b"789ABCDTUVWXYZopqrstuv"
    );
    assert_eq!(
        concat_slices(&records, &petri.umi, false).unwrap(),
        b"0123456"
    );
    assert_eq!(
        concat_sequence_slices(&records, &petri, false).unwrap(),
        b"abcdefghijklmnopq"
    );

    let surecell = technology_spec("SURECELL").unwrap();
    assert_eq!(
        concat_slices(&records, &surecell.bc, false).unwrap(),
        b"012345LMNOPQghijkl"
    );
    assert_eq!(
        concat_slices(&records, &surecell.umi, false).unwrap(),
        b"pqrstuvw"
    );

    let bdwta = technology_spec("BDWTA").unwrap();
    assert_eq!(
        concat_slices(&records, &bdwta.bc, false).unwrap(),
        b"012345678LMNOPQRSThijklmnop"
    );
    assert_eq!(
        concat_slices(&records, &bdwta.umi, false).unwrap(),
        b"qrstuvwx"
    );

    let vasa = technology_spec("VASA-SEQ").unwrap();
    assert_eq!(
        concat_slices(&records, &vasa.bc, false).unwrap(),
        b"6789ABCD"
    );
    assert_eq!(
        concat_slices(&records, &vasa.umi, false).unwrap(),
        b"012345"
    );
    assert_eq!(
        concat_sequence_slices(&records, &vasa, false).unwrap(),
        b"EFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"
    );

    let dropseq = technology_spec("DROPSEQ").unwrap();
    assert_eq!(
        concat_slices(&records, &dropseq.bc, false).unwrap(),
        b"0123456789AB"
    );
    assert_eq!(
        concat_slices(&records, &dropseq.umi, false).unwrap(),
        b"CDEFGHIJ"
    );

    let indrops_v1 = technology_spec("INDROPSV1").unwrap();
    assert_eq!(
        concat_slices(&records, &indrops_v1.bc, false).unwrap(),
        b"0123456789AUVWXYZab"
    );
    assert_eq!(
        concat_slices(&records, &indrops_v1.umi, false).unwrap(),
        b"ghijkl"
    );

    let indrops_v2 = technology_spec("INDROPSV2").unwrap();
    assert_eq!(
        concat_slices(&records, &indrops_v2.bc, false).unwrap(),
        b"abcdefghijkEFGHIJKL"
    );
    assert_eq!(
        concat_slices(&records, &indrops_v2.umi, false).unwrap(),
        b"QRSTUV"
    );

    let indrops_v3 = technology_spec("INDROPSV3").unwrap();
    assert_eq!(
        concat_slices(&records, &indrops_v3.bc, false).unwrap(),
        b"01234567abcdefgh"
    );
    assert_eq!(
        concat_slices(&records, &indrops_v3.umi, false).unwrap(),
        b"ijklmn"
    );
    assert_eq!(
        concat_sequence_slices(&records, &indrops_v3, false).unwrap(),
        b"QRST"
    );

    let celseq = technology_spec("CELSEQ").unwrap();
    assert_eq!(
        concat_slices(&records, &celseq.bc, false).unwrap(),
        b"01234567"
    );
    assert_eq!(
        concat_slices(&records, &celseq.umi, false).unwrap(),
        b"89AB"
    );

    let celseq2 = technology_spec("CELSEQ2").unwrap();
    assert_eq!(
        concat_slices(&records, &celseq2.bc, false).unwrap(),
        b"6789AB"
    );
    assert_eq!(
        concat_slices(&records, &celseq2.umi, false).unwrap(),
        b"012345"
    );

    let scrbseq = technology_spec("SCRBSEQ").unwrap();
    assert_eq!(
        concat_slices(&records, &scrbseq.bc, false).unwrap(),
        b"012345"
    );
    assert_eq!(
        concat_slices(&records, &scrbseq.umi, false).unwrap(),
        b"6789ABCDEF"
    );
}

#[test]
fn custom_technology_parser_accepts_kallisto_triples() {
    let spec = technology_spec("0,0,16:0,16,28:1,0,0").unwrap();
    assert_eq!(spec.nfiles, 2);
    assert_eq!(total_len(&spec.bc), Some(16));
    assert_eq!(total_len(&spec.umi), Some(12));
    assert!(!spec.paired);
}

#[test]
fn custom_technology_parser_keeps_at_least_one_input_file() {
    let spec = technology_spec("-1,-1,-1:-1,-1,-1:-1,-1,-1").unwrap();
    assert_eq!(spec.nfiles, 1);
    assert_eq!(total_len(&spec.bc), None);
    assert_eq!(total_len(&spec.umi), None);
}

#[test]
fn custom_technology_parser_rejects_malformed_sentinel_slices() {
    assert!(technology_spec("-1,0,0:-1,0,0:0,0,0,1,0,0").is_ok());

    let err = technology_spec("-1,1,4:0,0,4:1,0,0")
        .unwrap_err()
        .to_string();
    assert!(
        err.contains("sentinel technology slices must be -1,-1,-1 or -1,0,0"),
        "{err}"
    );
}

#[test]
fn custom_rx_technology_reads_umi_from_fastq_comments() {
    let spec = technology_spec("0,0,16:RX:1,0,0").unwrap();
    assert!(spec.keep_fastq_comments);
    assert!(!spec.bulk_like());
    assert_eq!(total_len(&spec.umi), None);
}

#[test]
fn custom_technology_parser_reports_colon_count_and_empty_lists() {
    let no_colons = technology_spec("0,0,16").unwrap_err().to_string();
    assert!(no_colons.contains("none found"), "{no_colons}");

    let one_colon = technology_spec("0,0,16:0,16,28").unwrap_err().to_string();
    assert!(one_colon.contains("only one found"), "{one_colon}");

    let three_colons = technology_spec("0,0,16:0,16,28:1,0,0:extra")
        .unwrap_err()
        .to_string();
    assert!(three_colons.contains("three found"), "{three_colons}");

    let empty_barcode = technology_spec(":0,16,28:1,0,0").unwrap_err().to_string();
    assert!(
        empty_barcode.contains("custom technology barcode list is empty"),
        "{empty_barcode}"
    );
}

#[test]
fn matrix_ec_uses_comma_separated_transcript_lists() {
    let dir = tempfile::tempdir().expect("tempdir");
    let path = dir.path().join("matrix.ec");

    write_matrix_ec(&path, &[vec![0, 2], vec![1]]).expect("write matrix.ec");

    assert_eq!(std::fs::read_to_string(path).unwrap(), "0\t0,2\n1\t1\n");
}

#[test]
fn transcripts_sidecar_uses_only_onlist_targets() {
    let dir = tempfile::tempdir().expect("tempdir");
    let path = dir.path().join("transcripts.txt");
    let names = vec![
        "tx0".to_string(),
        "tx1".to_string(),
        "d_list.0".to_string(),
        "tx2".to_string(),
    ];
    let onlist = vec![true, true, false, true];

    write_transcripts(&path, &names, Some(&onlist)).expect("write transcripts");

    assert_eq!(std::fs::read_to_string(path).unwrap(), "tx0\ntx1\ntx2\n");
    assert_eq!(onlist_target_count(&names, Some(&onlist)), 3);
    assert_eq!(onlist_target_count(&names, None), 4);
}

#[test]
fn all_listed_technology_names_are_accepted() {
    for name in [
        "10xv1",
        "10xv2",
        "10xv3",
        "10xv4",
        "Bulk",
        "ParseV3",
        "SmartSeq2",
        "BDWTA",
        "CELSeq",
        "CELSeq2",
        "DropSeq",
        "inDropsv1",
        "inDropsv2",
        "inDropsv3",
        "MATQSEQ",
        "PETRISEQ",
        "SCRBSeq",
        "SmartSeq3",
        "SPLiT-seq",
        "STORM-seq",
        "SureCell",
        "VASA-seq",
        "Visium",
    ] {
        let spec = technology_spec(name).unwrap_or_else(|err| {
            panic!("listed technology {name} should parse: {err}");
        });
        assert!(spec.nfiles > 0, "listed technology {name} should use reads");
        assert!(
            !spec.seq.is_empty(),
            "listed technology {name} should define cDNA slices"
        );
        assert!(
            !spec.bc.is_empty(),
            "listed technology {name} should define barcode slices"
        );
        assert!(
            !spec.umi.is_empty(),
            "listed technology {name} should define UMI slices"
        );
    }
}

#[test]
fn custom_two_sequence_slices_are_paired_only_when_requested() {
    let unpaired = technology_spec("-1,-1,-1:0,0,4:0,4,8,1,0,4").unwrap();
    assert!(!unpaired.paired);
    assert_eq!(unpaired.nfiles, 2);
    assert_eq!(unpaired.seq.len(), 2);

    let paired = technology_spec("-1,-1,-1:0,0,4:0,4,8,1,0,4%NONE%PAIRED").unwrap();
    assert!(paired.paired);
    assert_eq!(paired.nfiles, 2);
    assert_eq!(paired.seq.len(), 2);
}

#[test]
fn technology_suffix_applies_strand_and_paired_mode() {
    let forward = technology_spec("0,0,16:0,16,28:1,0,0%FORWARD").unwrap();
    assert_eq!(forward.default_strand, Some(StrandSpecific::FR));
    assert!(!forward.paired);

    let none_paired = technology_spec("0,0,16:0,16,28:1,0,0%NONE%PAIRED").unwrap();
    assert_eq!(none_paired.default_strand, None);
    assert!(none_paired.paired);

    let reverse_paired = technology_spec("0,0,16:0,16,28:1,0,0%REVERSE%PAIRED").unwrap();
    assert_eq!(reverse_paired.default_strand, Some(StrandSpecific::RF));
    assert!(reverse_paired.paired);
    assert_eq!(reverse_paired.nfiles, 3);
    assert_eq!(reverse_paired.seq.len(), 2);
    assert_eq!(reverse_paired.seq[1].file, Some(2));

    let invalid_strand = technology_spec("0,0,16:0,16,28:1,0,0%SIDEWAYS")
        .unwrap_err()
        .to_string();
    assert!(
        invalid_strand.contains("invalid technology strand suffix: SIDEWAYS"),
        "{invalid_strand}"
    );
    let invalid_pairing = technology_spec("0,0,16:0,16,28:1,0,0%FORWARD%MATED")
        .unwrap_err()
        .to_string();
    assert!(
        invalid_pairing.contains("invalid technology pairing suffix: MATED"),
        "{invalid_pairing}"
    );
    let extra = technology_spec("0,0,16:0,16,28:1,0,0%FORWARD%PAIRED%EXTRA")
        .unwrap_err()
        .to_string();
    assert!(
        extra.contains("unexpected technology suffix field: EXTRA"),
        "{extra}"
    );
}

#[test]
fn multiple_sequence_slices_are_separated_by_n() {
    let records = vec![
        FastqRecord {
            header: b"@r1".to_vec(),
            seq: b"AAAACCCC".to_vec(),
            plus: b"+".to_vec(),
            qual: b"IIIIIIII".to_vec(),
        },
        FastqRecord {
            header: b"@r2".to_vec(),
            seq: b"GGGGTTTT".to_vec(),
            plus: b"+".to_vec(),
            qual: b"IIIIIIII".to_vec(),
        },
    ];
    let spec = technology_spec("-1,-1,-1:0,0,4:0,4,8,1,0,4").unwrap();
    assert_eq!(
        concat_sequence_slices(&records, &spec, false).unwrap(),
        b"CCCCNGGGG"
    );
    assert_eq!(
        bus_query_quality(&records, &spec, false).unwrap(),
        b"IIII!IIII"
    );
}

#[test]
fn batch_long_flens_writer_preserves_each_batch_estimate() {
    let dir = tempfile::tempdir().expect("tempdir");
    let path = dir.path().join("flens.txt");
    let transcript_lengths = vec![100, 120];
    let sums = vec![vec![80, 0], vec![0, 90]];
    let counts = vec![vec![2, 0], vec![0, 3]];

    write_batch_long_flens(&path, &transcript_lengths, &sums, &counts, 31)
        .expect("write long flens");

    assert_eq!(std::fs::read_to_string(path).unwrap(), "9 89\n69 1\n");
}
