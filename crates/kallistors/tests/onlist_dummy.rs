use kallistors::index::{IndexBuildOptions, build_index};
use kallistors::io::FastqReader;
use kallistors::pseudoalign::{
    PseudoalignOptions, Strand, build_bifrost_index, pseudoalign_single_end_bifrost_with_options,
};
use std::collections::HashSet;
use std::fs;
use std::io::Cursor;

#[test]
fn deferred_dlist_cardinality_preserves_both_scanners() {
    let dir = tempfile::tempdir().expect("tempdir");
    let fasta = dir.path().join("transcripts.fa");
    let mut seed = 0x1234_5678u64;
    let sequences: Vec<Vec<u8>> = (0..2)
        .map(|_| {
            (0..150)
                .map(|_| {
                    seed ^= seed << 13;
                    seed ^= seed >> 7;
                    seed ^= seed << 17;
                    b"ACGT"[(seed & 3) as usize]
                })
                .collect()
        })
        .collect();
    fs::write(
        &fasta,
        format!(
            ">tx0\n{}\n>tx1\n{}\n",
            String::from_utf8_lossy(&sequences[0]),
            String::from_utf8_lossy(&sequences[1]),
        ),
    )
    .expect("write transcripts");
    let index_path = dir.path().join("transcripts.idx");
    build_index(&index_path, &[fasta], IndexBuildOptions::default()).expect("build index");
    let mut index = build_bifrost_index(&index_path).expect("load index");
    let kmer = &sequences[0][..index.k];
    let mut fwd = 0u64;
    let mut rev = 0u64;
    for (i, &base) in kmer.iter().enumerate() {
        let code = b"ACGT".iter().position(|&b| b == base).expect("base") as u64;
        fwd = (fwd << 2) | code;
        rev |= (3 - code) << (2 * i);
    }
    let mut reads = Vec::new();
    for (i, sequence) in sequences.iter().enumerate() {
        reads.extend_from_slice(format!("@r{i}\n").as_bytes());
        reads.extend_from_slice(&sequence[..index.k]);
        reads.extend_from_slice(b"\n+\n");
        reads.extend_from_slice(&vec![b'I'; index.k]);
        reads.push(b'\n');
    }

    for general_scanner in [false, true] {
        let options = PseudoalignOptions {
            kallisto_bifrost_find: general_scanner,
            ..PseudoalignOptions::default()
        };
        for dlist_present in [false, true] {
            index.dlist = dlist_present.then(|| HashSet::from([fwd.min(rev)]));
            let mut reader = FastqReader::new(Cursor::new(reads.clone()));
            let counts = pseudoalign_single_end_bifrost_with_options(
                &index,
                &mut reader,
                Strand::Unstranded,
                None,
                options,
            )
            .expect("pseudoalign");
            assert_eq!(counts.reads_processed, 2);
            assert_eq!(counts.reads_aligned, if dlist_present { 1 } else { 2 });
            if dlist_present {
                assert_eq!(counts.ec_list, [vec![1]]);
                assert_eq!(counts.counts, [1]);
            }
        }

        // With no on-list, no dummy hit is synthesized even when the D-list matches.
        let onlist = index.onlist.take();
        let mut reader = FastqReader::new(Cursor::new(reads.clone()));
        let counts = pseudoalign_single_end_bifrost_with_options(
            &index,
            &mut reader,
            Strand::Unstranded,
            None,
            options,
        )
        .expect("pseudoalign without on-list");
        assert_eq!(counts.reads_aligned, 2);
        index.onlist = onlist;
    }
}
