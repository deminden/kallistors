//! Transcript runs, annotation blocks and EC payload construction.

use super::*;

#[derive(Clone, Copy)]
pub(super) struct TrInfoEvent {
    pub(super) coord: u32,
    pub(super) is_start: bool,
    pub(super) trid: u32,
    pub(super) pos: u32,
}

pub(super) fn build_trinfos<L: ReferenceKmerLocations + ?Sized>(
    transcripts: &mut [PackedSequence],
    k: usize,
    locations: &L,
    paths: &[UnitigPath],
    threads: usize,
) -> Result<TrInfoShards> {
    let total_positions = transcripts
        .iter()
        .map(|tr| kmer_position_count(tr.len(), k))
        .sum::<usize>();
    if threads > 1
        && transcripts.len() > 1
        && total_positions >= PARALLEL_TRINFO_COLLECTION_MIN
        && !paths.is_empty()
    {
        let pool = IndexWorkers::new(threads)?;
        return pool.install(|| {
            transcripts
                .par_iter_mut()
                .enumerate()
                .map(|(tr_id, transcript)| {
                    let records =
                        collect_transcript_trinfos(tr_id as u32, transcript, k, locations, paths)?;
                    transcript.clear();
                    Ok(records)
                })
                .collect()
        });
    }
    let mut records = TrInfoShard::new(paths.len() <= (1 << 31));
    for (tr_id, transcript) in transcripts.iter_mut().enumerate() {
        scan_transcript_trinfos(
            tr_id as u32,
            transcript,
            k,
            locations,
            paths,
            |unitig_id, info| {
                records.push(unitig_id, info, (paths[unitig_id].seq.len() - k + 1) as u32)
            },
        )?;
        transcript.clear();
    }
    Ok(vec![records])
}

pub(super) fn collect_transcript_trinfos<L: ReferenceKmerLocations + ?Sized>(
    trid: u32,
    transcript: &PackedSequence,
    k: usize,
    locations: &L,
    paths: &[UnitigPath],
) -> Result<TrInfoShard> {
    let mut out = TrInfoShard::new(paths.len() <= (1 << 31));
    scan_transcript_trinfos(trid, transcript, k, locations, paths, |unitig_id, info| {
        out.push(unitig_id, info, (paths[unitig_id].seq.len() - k + 1) as u32);
    })?;
    Ok(out)
}

// The mapping context is shared by serial and parallel transcript collection.
pub(super) fn scan_transcript_trinfos<L: ReferenceKmerLocations + ?Sized, F>(
    trid: u32,
    transcript: &PackedSequence,
    k: usize,
    locations: &L,
    paths: &[UnitigPath],
    mut emit: F,
) -> Result<()>
where
    F: FnMut(usize, TRInfo),
{
    let seq = transcript;
    if seq.len() < k {
        return Ok(());
    }
    let mask = kmer_mask(k);
    let mut fwd = seq.word(0, k);
    let mut rev = reverse_complement_word(fwd, k, mask);
    let mut pos = 0;
    loop {
        let mut advance = 1;
        if let Some(loc) = locations.find(fwd.min(rev)) {
            let same_strand = (fwd <= rev) == loc.unitig_is_canonical();
            let unitig = &paths[loc.unitig_id()].seq;
            let offset = loc.offset();
            let remaining = seq.len() - pos - k;
            // Every interior edge of these compacted and joined unitigs has
            // exactly one successor. The matched kmer and strand therefore
            // prove the entire continuation up to the unitig/transcript end.
            // Branches end unitigs; no sampling or unchecked ambiguous jump.
            let matched = if same_strand {
                remaining.min(unitig.len() - offset - k)
            } else {
                remaining.min(offset)
            };
            #[cfg(test)]
            {
                let checked = (0..matched)
                    .take_while(|&i| {
                        let expected = if same_strand {
                            unitig[offset + k + i]
                        } else {
                            code_to_base(3 - u64::from(BASE_CODES[unitig[offset - 1 - i] as usize]))
                        };
                        seq.base(pos + k + i) == expected
                    })
                    .count();
                assert_eq!(
                    matched, checked,
                    "unitig interior must have an exact continuation"
                );
            }

            emit(
                loc.unitig_id(),
                TRInfo {
                    trid,
                    pos: pos as u32 | if same_strand { 0 } else { 0x8000_0000 },
                    start: (if same_strand {
                        offset
                    } else {
                        offset - matched
                    }) as u32,
                    stop: (if same_strand {
                        offset + matched + 1
                    } else {
                        offset + 1
                    }) as u32,
                },
            );
            advance += matched;
        }
        pos += advance;
        if pos > seq.len() - k {
            break;
        }
        if advance == 1 {
            let code = seq.code(pos + k - 1);
            fwd = shift_append_word(fwd, code, k, mask);
            rev = shift_prepend_word(rev, 3 - code, mask);
        } else {
            fwd = seq.word(pos, k);
            rev = reverse_complement_word(fwd, k, mask);
        }
    }
    Ok(())
}

pub(super) fn build_node_payloads(
    paths: &[UnitigPath],
    mut trinfos: TranscriptLocations,
    ec_threshold: usize,
    k: usize,
    threads: usize,
) -> Result<Vec<NodePayload>> {
    if trinfos.offsets.len() != paths.len() + 1 {
        return Err(Error::InvalidFormat(
            "invalid transcript-location table".into(),
        ));
    }
    let mut tail = trinfos.records.as_mut_slice();
    let mut slices = Vec::with_capacity(paths.len());
    for pair in trinfos.offsets.windows(2) {
        let (current, remaining) = tail.split_at_mut(pair[1] - pair[0]);
        slices.push(current);
        tail = remaining;
    }
    if threads > 1 && paths.len() >= 1024 {
        let pool = IndexWorkers::new(threads)?;
        return Ok(pool.install(|| {
            paths
                .par_iter()
                .zip(slices.into_par_iter())
                .enumerate()
                .map(|(id, (path, infos))| NodePayload {
                    id: id as u32,
                    blocks: build_blocks_for_unitig(
                        infos,
                        (path.seq.len().saturating_sub(k) + 1) as u32,
                        ec_threshold,
                    ),
                })
                .collect()
        }));
    }
    Ok(paths
        .iter()
        .zip(slices)
        .enumerate()
        .map(|(id, (path, infos))| NodePayload {
            id: id as u32,
            blocks: build_blocks_for_unitig(
                infos,
                (path.seq.len().saturating_sub(k) + 1) as u32,
                ec_threshold,
            ),
        })
        .collect())
}

pub(super) fn build_blocks_for_unitig(
    trinfos: &mut [TRInfo],
    kmer_count: u32,
    ec_threshold: usize,
) -> Vec<NodeBlock> {
    if trinfos.is_empty() {
        return vec![NodeBlock {
            lb: 0,
            ub: kmer_count,
            ec: Vec::new(),
            positions: BlockPositions::default(),
        }];
    }

    if kmer_count > 0
        && trinfos
            .iter()
            .all(|info| info.start == 0 && info.stop == kmer_count)
    {
        // A complete-unitig EC needs one block. Group exact positions directly,
        // reserving each transcript's actual group size rather than a tree.
        trinfos.sort_unstable_by_key(|info| (info.trid, info.pos));
        let count = 1 + trinfos
            .windows(2)
            .filter(|pair| pair[0].trid != pair[1].trid)
            .count();
        let capacity = if count <= ec_threshold { count } else { 0 };
        let (mut ec, mut positions) = (
            Vec::with_capacity(capacity),
            BlockPositions::with_capacity(capacity, if capacity > 0 { trinfos.len() } else { 0 }),
        );
        if count <= ec_threshold {
            let mut first = 0;
            while first < trinfos.len() {
                let end = first
                    + trinfos[first..].partition_point(|info| info.trid == trinfos[first].trid);
                ec.push(trinfos[first].trid);
                positions.push(trinfos[first..end].iter().map(|info| info.pos));
                first = end;
            }
        }
        return vec![NodeBlock {
            lb: 0,
            ub: kmer_count,
            ec,
            positions,
        }];
    }

    let mut breakpoints = Vec::with_capacity(trinfos.len() * 2 + 2);
    breakpoints.push(0);
    breakpoints.push(kmer_count);
    for info in trinfos.iter() {
        breakpoints.push(info.start);
        breakpoints.push(info.stop);
    }
    breakpoints.sort_unstable();
    breakpoints.dedup();

    // Small EC breakpoint sets are cheaper to scan directly; large ones use a sweep-line.
    if breakpoints.len() <= SWEEP_BLOCK_MIN_BREAKPOINTS {
        return build_blocks_for_unitig_by_scan(trinfos, &breakpoints, kmer_count, ec_threshold);
    }

    let mut events = Vec::with_capacity(trinfos.len() * 2);
    for info in trinfos.iter() {
        events.push(TrInfoEvent {
            coord: info.start,
            is_start: true,
            trid: info.trid,
            pos: info.pos,
        });
        events.push(TrInfoEvent {
            coord: info.stop,
            is_start: false,
            trid: info.trid,
            pos: info.pos,
        });
    }
    events.sort_unstable_by_key(|event| event.coord);

    let mut blocks = Vec::new();
    let mut event_idx = 0usize;
    let mut active: BTreeMap<u32, BTreeMap<u32, u32>> = BTreeMap::new();
    for win in breakpoints.windows(2) {
        let lb = win[0];
        let ub = win[1];
        let event_start = event_idx;
        while event_idx < events.len() && events[event_idx].coord == lb {
            event_idx += 1;
        }
        for event in events[event_start..event_idx]
            .iter()
            .filter(|event| !event.is_start)
        {
            remove_active_trinfo(&mut active, event.trid, event.pos);
        }
        for event in events[event_start..event_idx]
            .iter()
            .filter(|event| event.is_start)
        {
            *active
                .entry(event.trid)
                .or_default()
                .entry(event.pos)
                .or_default() += 1;
        }
        if lb == ub {
            continue;
        }
        let mut ec = Vec::new();
        let mut positions = BlockPositions::default();
        if active.len() <= ec_threshold {
            ec.reserve(active.len());
            positions = BlockPositions::with_capacity(
                active.len(),
                active.values().map(BTreeMap::len).sum(),
            );
            for (&trid, values) in &active {
                ec.push(trid);
                positions.push(values.keys().copied());
            }
        }
        push_or_merge_block(
            &mut blocks,
            NodeBlock {
                lb,
                ub,
                ec,
                positions,
            },
        );
    }

    if blocks.is_empty() {
        blocks.push(NodeBlock {
            lb: 0,
            ub: kmer_count,
            ec: Vec::new(),
            positions: BlockPositions::default(),
        });
    }
    blocks
}

pub(super) fn build_blocks_for_unitig_by_scan(
    trinfos: &mut [TRInfo],
    breakpoints: &[u32],
    kmer_count: u32,
    ec_threshold: usize,
) -> Vec<NodeBlock> {
    trinfos.sort_unstable_by_key(|info| (info.trid, info.pos, info.start, info.stop));
    let mut blocks = Vec::new();
    for win in breakpoints.windows(2) {
        let lb = win[0];
        let ub = win[1];
        if lb == ub {
            continue;
        }
        // Sorted TRInfo already groups transcript IDs and positions. A second
        // tree and per-transcript sort would reconstruct the same order.
        let mut ec = Vec::new();
        let mut positions = BlockPositions::default();
        for info in trinfos.iter() {
            if info.start <= lb && info.stop >= ub {
                if ec.last() != Some(&info.trid) {
                    ec.push(info.trid);
                    positions.start_group();
                }
                positions.push_last_unique(info.pos);
            }
        }
        if ec.len() > ec_threshold {
            ec.clear();
            positions.clear();
        }
        push_or_merge_block(
            &mut blocks,
            NodeBlock {
                lb,
                ub,
                ec,
                positions,
            },
        );
    }

    if blocks.is_empty() {
        blocks.push(NodeBlock {
            lb: 0,
            ub: kmer_count,
            ec: Vec::new(),
            positions: BlockPositions::default(),
        });
    }
    blocks
}

pub(super) fn remove_active_trinfo(
    active: &mut BTreeMap<u32, BTreeMap<u32, u32>>,
    trid: u32,
    pos: u32,
) {
    let Some(pos_counts) = active.get_mut(&trid) else {
        return;
    };
    let Some(count) = pos_counts.get_mut(&pos) else {
        return;
    };
    *count -= 1;
    if *count == 0 {
        pos_counts.remove(&pos);
    }
    if pos_counts.is_empty() {
        active.remove(&trid);
    }
}

pub(super) fn push_or_merge_block(blocks: &mut Vec<NodeBlock>, block: NodeBlock) {
    if let Some(last) = blocks.last_mut()
        && last.ub == block.lb
        && last.ec == block.ec
        && last.positions == block.positions
    {
        last.ub = block.ub;
        return;
    }
    blocks.push(block);
}
