//! Payloads regressions.

use super::*;
#[test]
fn block_sweep_matches_scan_builder() {
    let row = |trid| TRInfo {
        trid,
        start: 0,
        stop: 3,
        ..TRInfo::default()
    };
    let counts = [3u32, 1, 5, 1];
    for compact in [false, true] {
        let source = vec![vec![(2usize, row(1)), (0, row(2))], vec![(2, row(3))]];
        let shards = source
            .into_iter()
            .map(|records| {
                let mut shard = TrInfoShard::new(compact);
                for (unitig, info) in records {
                    shard.push(unitig, info, counts[unitig]);
                }
                shard
            })
            .collect();
        let flat = TranscriptLocations::from_shards(shards, &counts);
        assert_eq!(flat.offsets, [0, 1, 1, 3, 3]);
        assert_eq!(flat.records, [row(2), row(1), row(3)]);
    }
    let mut wide = TrInfoShard::new(false);
    wide.push(u32::MAX as usize, row(7), 3);
    assert_eq!(
        wide.iter(&[]).collect::<Vec<_>>(),
        [(u32::MAX as usize, row(7))]
    );

    let mut grouped = BlockPositions::default();
    grouped.push([]);
    grouped.push([1, 1, 2, 2]);
    grouped.push([2, 3]);
    assert_eq!(
        grouped.iter().collect::<Vec<_>>(),
        [vec![], vec![1, 2], vec![2, 3]]
    );
    grouped.clear();
    assert_eq!(grouped.len(), 0);
    grouped.push([u32::MAX]);
    assert_eq!(grouped.iter().collect::<Vec<_>>(), [vec![u32::MAX]]);
    let kmer_count = 64;
    let mut infos = (0..48)
        .map(|i| TRInfo {
            trid: (i % 7) as u32,
            pos: (i * 11) as u32,
            start: i as u32,
            stop: (i + 5) as u32,
        })
        .collect::<Vec<_>>();
    let mut breakpoints = vec![0, kmer_count];
    for info in &infos {
        breakpoints.push(info.start);
        breakpoints.push(info.stop);
    }
    breakpoints.sort_unstable();
    breakpoints.dedup();

    let mut scan_infos = infos.clone();
    let expected =
        build_blocks_for_unitig_by_scan(&mut scan_infos, &breakpoints, kmer_count, usize::MAX);
    let actual = build_blocks_for_unitig(&mut infos, kmer_count, usize::MAX);
    assert_eq!(actual, expected);
}
