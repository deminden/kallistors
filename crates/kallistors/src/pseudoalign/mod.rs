//! Pseudoalignment algorithms.

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, Read, Seek, SeekFrom};
use std::path::Path;

use boomphf::Mphf;

use crate::index::bifrost::{
    decode_pos_id, deserialize_roaring_to_vec, encode_minimizer_rep, read_bitcontainer_values,
    wyhash, BooPhf,
};

use crate::bias::{hexamer_to_int, BiasCounts};
use crate::index::{parse_graph_section, read_block_array_blocks};
use crate::io::ReadSource;
use crate::{Error, Result};

const KMER_BYTES_CANDIDATES: [usize; 4] = [8, 16, 24, 32];
const REP_HASH_VALS: [u64; 4] = [
    2053695854357871005,
    5073395517033431291,
    10060236952204337488,
    7783083932390163561,
];

/// A naive k-mer -> EC mapping built from a kallisto index.
pub struct KmerEcIndex {
    pub k: usize,
    pub mphf: Mphf<u64>,
    pub keys_by_index: Vec<u64>,
    pub ecs: Vec<Vec<u32>>,
    pub onlist: Option<Vec<bool>>,
}

pub struct BifrostIndex {
    pub k: usize,
    pub g: usize,
    pub unitigs: Vec<Vec<u8>>,
    pub km_unitigs: Vec<Vec<u8>>,
    pub h_kmers: Vec<Vec<u8>>,
    pub(crate) ec_blocks: Vec<Vec<crate::index::EcBlock>>,
    pub minz_positions: Vec<Vec<u64>>,
    pub mphf: BooPhf,
    pub onlist: Option<Vec<bool>>,
    pub transcript_lengths: Vec<u32>,
    pub transcript_names: Vec<String>,
    pub shade_to_color: Vec<Option<u32>>,
    pub shade_sequences: Vec<bool>,
    pub use_shade: bool,
}

/// Pseudoalignment EC counts.
pub struct EcCounts {
    pub ec_list: Vec<Vec<u32>>,
    pub counts: Vec<u32>,
    pub reads_processed: u64,
    pub reads_aligned: u64,
    pub bias: Option<BiasCounts>,
}

#[derive(Debug, Clone, Copy)]
pub struct FragmentFilter {
    pub fragment_length: u32,
    pub single_overhang: bool,
}

#[derive(Debug, Clone, Copy)]
pub struct PseudoalignOptions {
    pub min_range: usize,
    pub do_union: bool,
    pub dfk_onlist: bool,
    pub strand_specific: Option<StrandSpecific>,
    pub no_jump: bool,
    pub bias: bool,
    pub max_bias: usize,
}

impl Default for PseudoalignOptions {
    fn default() -> Self {
        Self {
            min_range: 1,
            do_union: false,
            dfk_onlist: false,
            strand_specific: None,
            no_jump: false,
            bias: false,
            max_bias: 1_000_000,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StrandSpecific {
    FR,
    RF,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum DebugFailReason {
    NoValidKmer,
    NoMinimizerHit,
    NoPositions,
    NoKmerMatch,
    EmptyEc,
    IntersectionEmpty,
    Unknown,
}

impl std::fmt::Display for DebugFailReason {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let s = match self {
            DebugFailReason::NoValidKmer => "no_valid_kmer",
            DebugFailReason::NoMinimizerHit => "no_minimizer_hit",
            DebugFailReason::NoPositions => "no_positions",
            DebugFailReason::NoKmerMatch => "no_kmer_match",
            DebugFailReason::EmptyEc => "empty_ec",
            DebugFailReason::IntersectionEmpty => "intersection_empty",
            DebugFailReason::Unknown => "unknown",
        };
        f.write_str(s)
    }
}

#[derive(Debug, Clone)]
pub struct ReadTrace {
    pub header: String,
    pub reason: DebugFailReason,
    pub kmer_pos: Option<usize>,
    pub min_pos: Option<usize>,
    pub kmer_seq: Option<String>,
    pub min_seq: Option<String>,
    pub positions: Option<usize>,
    pub sample_positions: Option<String>,
    pub used_revcomp: bool,
}

#[derive(Debug)]
pub struct DebugReport {
    pub counts: HashMap<DebugFailReason, u64>,
    pub traces: Vec<ReadTrace>,
    pub max_traces: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Strand {
    Unstranded,
    Forward,
    Reverse,
}

impl DebugReport {
    fn new(max_traces: usize) -> Self {
        Self {
            counts: HashMap::new(),
            traces: Vec::new(),
            max_traces,
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn record(
        &mut self,
        header: &[u8],
        reason: DebugFailReason,
        kmer_pos: Option<usize>,
        min_pos: Option<usize>,
        kmer_seq: Option<String>,
        min_seq: Option<String>,
        positions: Option<usize>,
        sample_positions: Option<String>,
        used_revcomp: bool,
    ) {
        *self.counts.entry(reason).or_insert(0) += 1;
        if self.traces.len() < self.max_traces {
            let header = String::from_utf8_lossy(header).into_owned();
            self.traces.push(ReadTrace {
                header,
                reason,
                kmer_pos,
                min_pos,
                kmer_seq,
                min_seq,
                positions,
                sample_positions,
                used_revcomp,
            });
        }
    }
}

#[derive(Default)]
struct ReadDebugState {
    saw_valid_kmer: bool,
    saw_mphf_hit: bool,
    saw_positions: bool,
    saw_match: bool,
    saw_ec: bool,
    intersection_empty: bool,
    first_mphf_miss: Option<(usize, usize)>,
    first_no_positions: Option<(usize, usize)>,
    first_no_match: Option<(usize, usize)>,
    first_empty_ec: Option<(usize, usize)>,
    first_no_match_positions: Option<(usize, usize, usize, String)>,
    used_revcomp: bool,
}

#[derive(Debug, Clone, Copy)]
struct MatchInfo {
    unitig_id: usize,
    unitig_pos: usize,
    read_pos: usize,
    used_revcomp: bool,
}

#[derive(Debug, Clone)]
struct Hit {
    unitig_id: usize,
    read_pos: usize,
    block_idx: usize,
    used_revcomp: bool,
}

#[derive(Debug)]
struct ReadEc {
    ec: Vec<u32>,
    best_match: Option<MatchInfo>,
    first_hit: Option<Hit>,
    hits: Vec<Hit>,
    had_offlist: bool,
    shade_union: Vec<u32>,
}

/// Build a naive k-mer -> EC map by scanning unitigs and EC blocks in the index.
pub fn build_kmer_ec_index(path: &Path) -> Result<KmerEcIndex> {
    let file = File::open(path).map_err(|_| Error::MissingFile(path.to_path_buf()))?;
    let mut reader = BufReader::new(file);

    let _index_version = read_u64_le(&mut reader)?;
    let dbg_size = read_u64_le(&mut reader)?;

    let mut k: Option<u32> = None;
    let mut kmer_bytes: Option<usize> = None;

    let sequences = if dbg_size > 0 {
        let dbg_start = reader.stream_position()?;
        let meta = parse_graph_section(&mut reader, dbg_start, dbg_size, &KMER_BYTES_CANDIDATES)?;
        k = Some(meta.k);
        kmer_bytes = Some(meta.kmer_bytes);

        reader.seek(SeekFrom::Start(dbg_start))?;
        let seqs =
            read_unitig_sequences(&mut reader, dbg_start, dbg_size, meta.k, meta.kmer_bytes)?;
        reader.seek(SeekFrom::Start(dbg_start + dbg_size))?;

        let mphf_size = read_u64_le(&mut reader)?;
        reader.seek(SeekFrom::Current(mphf_size as i64))?;
        seqs
    } else {
        Vec::new()
    };

    let dlist_size = read_u64_le(&mut reader)?;
    let _dlist_overhang = read_u64_le(&mut reader)?;
    let kmer_bytes =
        kmer_bytes.ok_or_else(|| Error::InvalidFormat("missing k-mer width".into()))?;
    let dlist_skip = dlist_size
        .checked_mul(kmer_bytes as u64)
        .ok_or_else(|| Error::InvalidFormat("d-list size overflow".into()))?;
    reader.seek(SeekFrom::Current(dlist_skip as i64))?;

    let node_count = read_u64_le(&mut reader)?;
    let k = k.ok_or_else(|| Error::InvalidFormat("missing k-mer length".into()))? as usize;
    if k > 32 {
        return Err(Error::UnsupportedFeature(
            "k-mer length > 32 not supported in MPHF path".into(),
        ));
    }

    if sequences.len() != node_count as usize {
        return Err(Error::InvalidFormat(
            "unitig and node counts mismatch".into(),
        ));
    }

    let mut map: HashMap<u64, Vec<u32>> = HashMap::new();

    for (idx, seq) in sequences.iter().enumerate() {
        let _ = idx;
        reader.seek(SeekFrom::Current(k as i64))?;
        let node_size = read_u32_le(&mut reader)? as usize;
        let mut buf = vec![0u8; node_size];
        reader.read_exact(&mut buf)?;

        let mut cur = std::io::Cursor::new(buf.as_slice());
        let _node_id = read_u32_le(&mut cur)?;
        let blocks = read_block_array_blocks(&mut cur)?;
        for block in blocks {
            let start = block.lb as usize;
            let end = block.ub as usize;
            for pos in start..end {
                if pos + k > seq.len() {
                    break;
                }
                if let Some((fwd, rev)) = encode_kmer_pair(&seq[pos..pos + k]) {
                    map.entry(fwd)
                        .and_modify(|existing| merge_sorted_unique(existing, &block.ec))
                        .or_insert_with(|| block.ec.clone());
                    if rev != fwd {
                        map.entry(rev)
                            .and_modify(|existing| merge_sorted_unique(existing, &block.ec))
                            .or_insert_with(|| block.ec.clone());
                    }
                }
            }
        }
    }

    let mut keys: Vec<u64> = map.keys().cloned().collect();
    keys.sort_unstable();
    let mphf = Mphf::new(1.7, &keys);
    let mut keys_by_index = vec![0u64; keys.len()];
    let mut ecs = vec![Vec::new(); keys.len()];
    for (key, ec) in map {
        let idx = mphf.hash(&key) as usize;
        keys_by_index[idx] = key;
        ecs[idx] = ec;
    }

    let onlist = read_onlist(&mut reader)?;

    Ok(KmerEcIndex {
        k,
        mphf,
        keys_by_index,
        ecs,
        onlist,
    })
}

pub fn build_bifrost_index(path: &Path) -> Result<BifrostIndex> {
    build_bifrost_index_with_positions(path, false)
}

pub fn build_bifrost_index_with_positions(
    path: &Path,
    load_positional_info: bool,
) -> Result<BifrostIndex> {
    let file = File::open(path).map_err(|_| Error::MissingFile(path.to_path_buf()))?;
    let mut reader = BufReader::new(file);

    let _index_version = read_u64_le(&mut reader)?;
    let dbg_size = read_u64_le(&mut reader)?;
    let dbg_start = reader.stream_position()?;

    let file_format_version = read_u64_le(&mut reader)?;
    let k = read_i32_le(&mut reader)?;
    let g = read_i32_le(&mut reader)?;
    let _ = file_format_version;

    if k <= 0 || g <= 0 {
        return Err(Error::InvalidFormat("invalid k/g".into()));
    }
    let k = k as usize;
    let g = g as usize;
    if k > 32 || g > 32 {
        return Err(Error::UnsupportedFeature(
            "k/g > 32 not supported in Bifrost path".into(),
        ));
    }

    let v_unitigs_sz = read_u64_le(&mut reader)? as usize;
    let mut unitigs = Vec::with_capacity(v_unitigs_sz);
    for _ in 0..v_unitigs_sz {
        let len = read_u64_le(&mut reader)? as usize;
        let data_sz = len.div_ceil(4);
        let mut data = vec![0u8; data_sz];
        reader.read_exact(&mut data)?;
        unitigs.push(decode_compressed_sequence(&data, len));
    }

    let km_unitigs_sz = read_u64_le(&mut reader)? as usize;
    let mut km_unitigs = Vec::with_capacity(km_unitigs_sz);
    for _ in 0..km_unitigs_sz {
        let mut buf = vec![0u8; 8];
        reader.read_exact(&mut buf)?;
        km_unitigs.push(decode_kmer(&buf, k));
    }

    let h_kmers_ccov_sz = read_u64_le(&mut reader)? as usize;
    let mut h_kmers = Vec::with_capacity(h_kmers_ccov_sz);
    for _ in 0..h_kmers_ccov_sz {
        let mut buf = vec![0u8; 8];
        reader.read_exact(&mut buf)?;
        h_kmers.push(decode_kmer(&buf, k));
    }

    let _bfg_header = read_u64_le(&mut reader)?;
    let _checksum = read_u64_le(&mut reader)?;
    let _v_unitigs_sz2 = read_u64_le(&mut reader)?;
    let _km_unitigs_sz2 = read_u64_le(&mut reader)?;
    let _h_kmers_sz2 = read_u64_le(&mut reader)?;
    let hmap_min_unitigs_sz = read_u64_le(&mut reader)? as usize;

    let nb_bmp_unitigs = read_u64_le(&mut reader)? as usize;
    let mut bmp_unitigs = Vec::with_capacity(nb_bmp_unitigs);
    for _ in 0..nb_bmp_unitigs {
        bmp_unitigs.push(read_bitcontainer_values(&mut reader)?);
    }

    let nb_bmp_km = read_u64_le(&mut reader)? as usize;
    let mut bmp_km = Vec::with_capacity(nb_bmp_km);
    for _ in 0..nb_bmp_km {
        bmp_km.push(read_bitcontainer_values(&mut reader)?);
    }

    let nb_special_minz = read_u64_le(&mut reader)? as usize;
    let nb_bmp_special = (nb_special_minz >> 32) + 1;
    let mut abundant = std::collections::HashSet::new();
    let mut overcrowded = std::collections::HashSet::new();
    for _ in 0..nb_bmp_special {
        for id in read_bitcontainer_values(&mut reader)? {
            abundant.insert(id as usize);
        }
    }
    for _ in 0..nb_bmp_special {
        for id in read_bitcontainer_values(&mut reader)? {
            overcrowded.insert(id as usize);
        }
    }

    reader.seek(SeekFrom::Start(dbg_start + dbg_size))?;

    let mphf_size = read_u64_le(&mut reader)? as usize;
    let mut mphf_buf = vec![0u8; mphf_size];
    reader.read_exact(&mut mphf_buf)?;
    let mut mphf_cur = std::io::Cursor::new(mphf_buf.as_slice());
    let mphf = BooPhf::load(&mut mphf_cur)?;

    let dlist_size = read_u64_le(&mut reader)?;
    let _dlist_overhang = read_u64_le(&mut reader)?;
    let dlist_skip = dlist_size
        .checked_mul(8)
        .ok_or_else(|| Error::InvalidFormat("d-list size overflow".into()))?;
    reader.seek(SeekFrom::Current(dlist_skip as i64))?;

    let node_count = read_u64_le(&mut reader)? as usize;
    let mut ec_blocks = Vec::with_capacity(node_count);
    for _ in 0..node_count {
        reader.seek(SeekFrom::Current(k as i64))?;
        let node_size = read_u32_le(&mut reader)? as usize;
        let mut buf = vec![0u8; node_size];
        reader.read_exact(&mut buf)?;
        let mut cur = std::io::Cursor::new(buf.as_slice());
        let _node_id = read_u32_le(&mut cur)?;
        let blocks = if load_positional_info {
            crate::index::read_block_array_blocks_with_positions(&mut cur)?
        } else {
            read_block_array_blocks(&mut cur)?
        };
        ec_blocks.push(blocks);
    }

    let mut minz_positions = vec![Vec::<u64>::new(); hmap_min_unitigs_sz];

    let mut prefix = Vec::with_capacity(unitigs.len());
    let mut acc = 0u64;
    for seq in &unitigs {
        let count = seq.len().saturating_sub(g) + 1;
        acc += count as u64;
        prefix.push(acc);
    }

    for (i, bitmap) in bmp_unitigs.iter().enumerate() {
        let base = (i as u64) << 32;
        for &pos_bmp in bitmap {
            let pos = base + pos_bmp as u64;
            let unitig_id = match prefix.binary_search(&(pos + 1)) {
                Ok(idx) => idx,
                Err(idx) => idx,
            };
            let prev = if unitig_id == 0 {
                0
            } else {
                prefix[unitig_id - 1]
            };
            let rel = (pos - prev) as usize;
            if let Some(bytes) = encode_minimizer_rep(&unitigs[unitig_id][rel..rel + g]) {
                if let Some(idx) = mphf.lookup(&bytes) {
                    minz_positions[idx as usize].push(((unitig_id as u64) << 32) | rel as u64);
                }
            }
        }
    }

    let km_glen = k.saturating_sub(g) + 1;
    for (i, bitmap) in bmp_km.iter().enumerate() {
        let base = (i as u64) << 32;
        for &pos_bmp in bitmap {
            let pos = base + pos_bmp as u64;
            let km_id = (pos / km_glen as u64) as usize;
            let km_pos = (pos % km_glen as u64) as usize;
            if km_id >= km_unitigs.len() || km_pos + g > k {
                continue;
            }
            if let Some(bytes) = encode_minimizer_rep(&km_unitigs[km_id][km_pos..km_pos + g]) {
                if let Some(idx) = mphf.lookup(&bytes) {
                    let pos_id = ((km_id as u64) << 32) | 0x8000_0000 | (km_pos as u64);
                    minz_positions[idx as usize].push(pos_id);
                }
            }
        }
    }

    for i in 0..nb_special_minz {
        let minz = read_u64_le(&mut reader)?;
        let bytes = minz.to_le_bytes();
        let is_abundant = abundant.contains(&i);
        let is_overcrowded = overcrowded.contains(&i);
        let pos_id = if is_abundant || is_overcrowded {
            let mut id = 0xffff_ffff_0000_0000 | if is_overcrowded { 0x8000_0000 } else { 0 };
            if is_abundant {
                let count = read_u32_le(&mut reader)? as u64;
                id |= count;
            }
            id
        } else {
            read_u64_le(&mut reader)?
        };
        if let Some(idx) = mphf.lookup(&bytes) {
            minz_positions[idx as usize].push(pos_id);
        }
    }

    let onlist = read_onlist_with_lengths(&mut reader)?;

    Ok(BifrostIndex {
        k,
        g,
        unitigs,
        km_unitigs,
        h_kmers,
        ec_blocks,
        minz_positions,
        mphf,
        onlist: onlist.onlist,
        transcript_lengths: onlist.lengths,
        transcript_names: onlist.names,
        shade_to_color: onlist.shade_to_color,
        shade_sequences: onlist.shade_sequences,
        use_shade: onlist.use_shade,
    })
}

/// Pseudoalign single-end reads using a naive k-mer map.
pub fn pseudoalign_single_end<R: ReadSource>(
    index: &KmerEcIndex,
    reader: &mut R,
) -> Result<EcCounts> {
    let mut ec_list: Vec<Vec<u32>> = Vec::new();
    let mut ec_map: HashMap<Vec<u32>, usize> = HashMap::new();
    let mut counts: Vec<u32> = Vec::new();
    let mut reads_processed = 0u64;
    let mut reads_aligned = 0u64;

    let mut current: Vec<u32> = Vec::new();
    let mut next: Vec<u32> = Vec::new();

    while let Some(record) = reader.next_record() {
        let record = record?;
        reads_processed += 1;
        if record.seq.len() < index.k {
            continue;
        }

        let mut has_hit = false;
        current.clear();

        for pos in 0..=record.seq.len() - index.k {
            if let Some((fwd, rev)) = encode_kmer_pair(&record.seq[pos..pos + index.k]) {
                if let Some(ec) = lookup_ec(index, fwd).or_else(|| lookup_ec(index, rev)) {
                    if !has_hit {
                        current.extend_from_slice(ec);
                        has_hit = true;
                    } else {
                        next.clear();
                        intersect_sorted(&current, ec, &mut next);
                        std::mem::swap(&mut current, &mut next);
                        if current.is_empty() {
                            break;
                        }
                    }
                }
            }
        }

        if !has_hit || current.is_empty() {
            continue;
        }

        reads_aligned += 1;
        let ec_id = match ec_map.get(&current) {
            Some(id) => *id,
            None => {
                let id = ec_list.len();
                ec_map.insert(current.clone(), id);
                ec_list.push(current.clone());
                counts.push(0);
                id
            }
        };
        counts[ec_id] += 1;
    }

    Ok(EcCounts {
        ec_list,
        counts,
        reads_processed,
        reads_aligned,
        bias: None,
    })
}

pub fn pseudoalign_paired_naive<R1: ReadSource, R2: ReadSource>(
    index: &KmerEcIndex,
    reader1: &mut R1,
    reader2: &mut R2,
) -> Result<EcCounts> {
    let mut ec_list: Vec<Vec<u32>> = Vec::new();
    let mut ec_map: HashMap<Vec<u32>, usize> = HashMap::new();
    let mut counts: Vec<u32> = Vec::new();
    let mut reads_processed = 0u64;
    let mut reads_aligned = 0u64;

    loop {
        let r1 = reader1.next_record();
        let r2 = reader2.next_record();
        match (r1, r2) {
            (None, None) => break,
            (Some(_), None) | (None, Some(_)) => {
                return Err(Error::InvalidFormat("paired FASTQ length mismatch".into()))
            }
            (Some(a), Some(b)) => {
                let a = a?;
                let b = b?;
                reads_processed += 1;
                let ec1 = ec_for_read_naive(index, &a.seq);
                let ec2 = ec_for_read_naive(index, &b.seq);
                let (Some(ec1), Some(ec2)) = (ec1, ec2) else {
                    continue;
                };
                let mut merged = Vec::new();
                intersect_sorted(&ec1, &ec2, &mut merged);
                if merged.is_empty() {
                    continue;
                }
                reads_aligned += 1;
                let ec_id = match ec_map.get(&merged) {
                    Some(id) => *id,
                    None => {
                        let id = ec_list.len();
                        ec_map.insert(merged.clone(), id);
                        ec_list.push(merged.clone());
                        counts.push(0);
                        id
                    }
                };
                counts[ec_id] += 1;
            }
        }
    }

    Ok(EcCounts {
        ec_list,
        counts,
        reads_processed,
        reads_aligned,
        bias: None,
    })
}

pub fn pseudoalign_single_end_bifrost<R: ReadSource>(
    index: &BifrostIndex,
    reader: &mut R,
) -> Result<EcCounts> {
    pseudoalign_single_end_bifrost_inner(
        index,
        reader,
        Strand::Unstranded,
        None,
        None,
        PseudoalignOptions::default(),
    )
}

pub fn pseudoalign_single_end_bifrost_debug<R: ReadSource>(
    index: &BifrostIndex,
    reader: &mut R,
    max_traces: usize,
) -> Result<(EcCounts, DebugReport)> {
    let mut report = DebugReport::new(max_traces);
    let counts = pseudoalign_single_end_bifrost_inner(
        index,
        reader,
        Strand::Unstranded,
        Some(&mut report),
        None,
        PseudoalignOptions::default(),
    )?;
    Ok((counts, report))
}

pub fn pseudoalign_single_end_bifrost_with_strand<R: ReadSource>(
    index: &BifrostIndex,
    reader: &mut R,
    strand: Strand,
) -> Result<EcCounts> {
    pseudoalign_single_end_bifrost_inner(
        index,
        reader,
        strand,
        None,
        None,
        PseudoalignOptions::default(),
    )
}

pub fn pseudoalign_single_end_bifrost_debug_with_strand<R: ReadSource>(
    index: &BifrostIndex,
    reader: &mut R,
    strand: Strand,
    max_traces: usize,
) -> Result<(EcCounts, DebugReport)> {
    let mut report = DebugReport::new(max_traces);
    let counts = pseudoalign_single_end_bifrost_inner(
        index,
        reader,
        strand,
        Some(&mut report),
        None,
        PseudoalignOptions::default(),
    )?;
    Ok((counts, report))
}

pub fn pseudoalign_single_end_bifrost_with_strand_and_filter<R: ReadSource>(
    index: &BifrostIndex,
    reader: &mut R,
    strand: Strand,
    filter: Option<FragmentFilter>,
) -> Result<EcCounts> {
    pseudoalign_single_end_bifrost_inner(
        index,
        reader,
        strand,
        None,
        filter,
        PseudoalignOptions::default(),
    )
}

pub fn pseudoalign_single_end_bifrost_debug_with_strand_and_filter<R: ReadSource>(
    index: &BifrostIndex,
    reader: &mut R,
    strand: Strand,
    max_traces: usize,
    filter: Option<FragmentFilter>,
) -> Result<(EcCounts, DebugReport)> {
    let mut report = DebugReport::new(max_traces);
    let counts = pseudoalign_single_end_bifrost_inner(
        index,
        reader,
        strand,
        Some(&mut report),
        filter,
        PseudoalignOptions::default(),
    )?;
    Ok((counts, report))
}

pub fn pseudoalign_single_end_bifrost_with_options<R: ReadSource>(
    index: &BifrostIndex,
    reader: &mut R,
    strand: Strand,
    filter: Option<FragmentFilter>,
    options: PseudoalignOptions,
) -> Result<EcCounts> {
    pseudoalign_single_end_bifrost_inner(index, reader, strand, None, filter, options)
}

pub fn pseudoalign_single_end_bifrost_debug_with_options<R: ReadSource>(
    index: &BifrostIndex,
    reader: &mut R,
    strand: Strand,
    max_traces: usize,
    filter: Option<FragmentFilter>,
    options: PseudoalignOptions,
) -> Result<(EcCounts, DebugReport)> {
    let mut report = DebugReport::new(max_traces);
    let counts = pseudoalign_single_end_bifrost_inner(
        index,
        reader,
        strand,
        Some(&mut report),
        filter,
        options,
    )?;
    Ok((counts, report))
}

pub fn pseudoalign_paired_bifrost_with_strand<R1: ReadSource, R2: ReadSource>(
    index: &BifrostIndex,
    reader1: &mut R1,
    reader2: &mut R2,
    strand: Strand,
) -> Result<EcCounts> {
    pseudoalign_paired_bifrost_inner(
        index,
        reader1,
        reader2,
        strand,
        None,
        PseudoalignOptions::default(),
    )
}

pub fn pseudoalign_paired_bifrost_debug_with_strand<R1: ReadSource, R2: ReadSource>(
    index: &BifrostIndex,
    reader1: &mut R1,
    reader2: &mut R2,
    strand: Strand,
    max_traces: usize,
) -> Result<(EcCounts, DebugReport)> {
    let mut report = DebugReport::new(max_traces);
    let counts = pseudoalign_paired_bifrost_inner(
        index,
        reader1,
        reader2,
        strand,
        Some(&mut report),
        PseudoalignOptions::default(),
    )?;
    Ok((counts, report))
}

pub fn pseudoalign_paired_bifrost_with_options<R1: ReadSource, R2: ReadSource>(
    index: &BifrostIndex,
    reader1: &mut R1,
    reader2: &mut R2,
    strand: Strand,
    options: PseudoalignOptions,
) -> Result<EcCounts> {
    pseudoalign_paired_bifrost_inner(index, reader1, reader2, strand, None, options)
}

pub fn pseudoalign_paired_bifrost_debug_with_options<R1: ReadSource, R2: ReadSource>(
    index: &BifrostIndex,
    reader1: &mut R1,
    reader2: &mut R2,
    strand: Strand,
    max_traces: usize,
    options: PseudoalignOptions,
) -> Result<(EcCounts, DebugReport)> {
    let mut report = DebugReport::new(max_traces);
    let counts = pseudoalign_paired_bifrost_inner(
        index,
        reader1,
        reader2,
        strand,
        Some(&mut report),
        options,
    )?;
    Ok((counts, report))
}

fn pseudoalign_paired_bifrost_inner<R1: ReadSource, R2: ReadSource>(
    index: &BifrostIndex,
    reader1: &mut R1,
    reader2: &mut R2,
    strand: Strand,
    mut report: Option<&mut DebugReport>,
    options: PseudoalignOptions,
) -> Result<EcCounts> {
    let mut ec_list: Vec<Vec<u32>> = Vec::new();
    let mut ec_map: HashMap<Vec<u32>, usize> = HashMap::new();
    let mut counts: Vec<u32> = Vec::new();
    let mut reads_processed = 0u64;
    let mut reads_aligned = 0u64;
    let mut bias = if options.bias {
        Some(BiasCounts::new())
    } else {
        None
    };

    loop {
        let r1 = reader1.next_record();
        let r2 = reader2.next_record();
        match (r1, r2) {
            (None, None) => break,
            (Some(_), None) | (None, Some(_)) => {
                return Err(Error::InvalidFormat("paired FASTQ length mismatch".into()))
            }
            (Some(a), Some(b)) => {
                let a = a?;
                let b = b?;
                reads_processed += 1;

                let mut dbg1 = report.as_ref().map(|_| ReadDebugState::default());
                let mut dbg2 = report.as_ref().map(|_| ReadDebugState::default());
                let ec1 = ec_for_read_bifrost(index, &a.seq, strand, dbg1.as_mut(), options);
                let ec2 = ec_for_read_bifrost(index, &b.seq, strand, dbg2.as_mut(), options);

                if ec1.is_none() || ec2.is_none() {
                    if let Some(r) = report.as_deref_mut() {
                        let (state, header) = if ec1.is_none() {
                            (dbg1, &a.header)
                        } else {
                            (dbg2, &b.header)
                        };
                        if let Some(state) = state {
                            let (reason, kmer_pos, min_pos) = debug_reason(&state);
                            let (positions, sample_positions) = state
                                .first_no_match_positions
                                .as_ref()
                                .map(|v| (Some(v.2), Some(v.3.clone())))
                                .unwrap_or((None, None));
                            r.record(
                                header,
                                reason,
                                kmer_pos,
                                min_pos,
                                None,
                                None,
                                positions,
                                sample_positions,
                                state.used_revcomp,
                            );
                        }
                    }
                    continue;
                }
                let ec1 = ec1.unwrap();
                let ec2 = ec2.unwrap();
                let mut ec1 = ec1;
                let mut ec2 = ec2;
                if let Some(mode) = options.strand_specific {
                    let comprehensive = options.do_union || options.no_jump || index.use_shade;
                    if let Some(filtered) = apply_strand_filter(
                        index,
                        &ec1,
                        mode,
                        true,
                        comprehensive,
                        &mut report,
                        &a.header,
                    ) {
                        ec1.ec = filtered;
                    }
                    if let Some(filtered) = apply_strand_filter(
                        index,
                        &ec2,
                        mode,
                        false,
                        comprehensive,
                        &mut report,
                        &b.header,
                    ) {
                        ec2.ec = filtered;
                    }
                }
                if ec1.ec.is_empty() || ec2.ec.is_empty() {
                    continue;
                }

                let use_shade = index.use_shade;
                if use_shade && options.do_union {
                    filter_shades_in_place(&mut ec1.ec, &index.shade_sequences);
                    filter_shades_in_place(&mut ec2.ec, &index.shade_sequences);
                }

                let mut merged = Vec::new();
                intersect_sorted(&ec1.ec, &ec2.ec, &mut merged);
                if options.dfk_onlist && (ec1.had_offlist || ec2.had_offlist) && !merged.is_empty()
                {
                    if let Some(onlist) = index.onlist.as_deref() {
                        let dummy = onlist.len() as u32;
                        if merged.last().copied() != Some(dummy) {
                            merged.push(dummy);
                        }
                    }
                }
                if merged.is_empty() {
                    if let Some(r) = report.as_deref_mut() {
                        let mut header = Vec::new();
                        header.extend_from_slice(&a.header);
                        header.extend_from_slice(b"|");
                        header.extend_from_slice(&b.header);
                        r.record(
                            &header,
                            DebugFailReason::IntersectionEmpty,
                            None,
                            None,
                            None,
                            None,
                            None,
                            None,
                            dbg1.as_ref().map(|d| d.used_revcomp).unwrap_or(false)
                                || dbg2.as_ref().map(|d| d.used_revcomp).unwrap_or(false),
                        );
                    }
                    continue;
                }

                if use_shade && !merged.is_empty() {
                    let mut shade_candidates = Vec::new();
                    merge_sorted_unique_vec(&mut shade_candidates, &ec1.shade_union);
                    merge_sorted_unique_vec(&mut shade_candidates, &ec2.shade_union);
                    if !shade_candidates.is_empty() {
                        for shade in shade_candidates {
                            let Some(color) =
                                index.shade_to_color.get(shade as usize).copied().flatten()
                            else {
                                continue;
                            };
                            if merged.binary_search(&color).is_ok() {
                                merged.push(shade);
                            }
                        }
                        merged.sort_unstable();
                        merged.dedup();
                    }
                }

                if let Some(bias_counts) = bias.as_mut() {
                    if bias_counts.total < options.max_bias as u64 {
                        if let Some(best_match) = ec1.best_match {
                            if let Some(hex) = bias_hexamer_for_match(index, best_match) {
                                bias_counts.record(hex);
                            }
                        }
                    }
                }

                reads_aligned += 1;
                let ec_id = match ec_map.get(&merged) {
                    Some(id) => *id,
                    None => {
                        let id = ec_list.len();
                        ec_map.insert(merged.clone(), id);
                        ec_list.push(merged.clone());
                        counts.push(0);
                        id
                    }
                };
                counts[ec_id] += 1;
            }
        }
    }

    Ok(EcCounts {
        ec_list,
        counts,
        reads_processed,
        reads_aligned,
        bias,
    })
}

fn pseudoalign_single_end_bifrost_inner<R: ReadSource>(
    index: &BifrostIndex,
    reader: &mut R,
    strand: Strand,
    mut report: Option<&mut DebugReport>,
    filter: Option<FragmentFilter>,
    options: PseudoalignOptions,
) -> Result<EcCounts> {
    let mut ec_list: Vec<Vec<u32>> = Vec::new();
    let mut ec_map: HashMap<Vec<u32>, usize> = HashMap::new();
    let mut counts: Vec<u32> = Vec::new();
    let mut reads_processed = 0u64;
    let mut reads_aligned = 0u64;
    let mut bias = if options.bias {
        Some(BiasCounts::new())
    } else {
        None
    };

    while let Some(record) = reader.next_record() {
        let record = record?;
        reads_processed += 1;

        let mut dbg = report.as_ref().map(|_| ReadDebugState::default());
        let ec = ec_for_read_bifrost(index, &record.seq, strand, dbg.as_mut(), options);
        let Some(mut read_ec) = ec else {
            if let Some(r) = report.as_deref_mut() {
                if let Some(state) = dbg {
                    let (reason, kmer_pos, min_pos) = debug_reason(&state);
                    let kmer_seq = kmer_pos.and_then(|pos| {
                        if pos + index.k <= record.seq.len() {
                            Some(
                                String::from_utf8_lossy(&record.seq[pos..pos + index.k])
                                    .into_owned(),
                            )
                        } else {
                            None
                        }
                    });
                    let min_seq = match (kmer_pos, min_pos) {
                        (Some(kp), Some(mp)) => {
                            let start = kp + mp;
                            let end = start + index.g;
                            if end <= record.seq.len() {
                                Some(String::from_utf8_lossy(&record.seq[start..end]).into_owned())
                            } else {
                                None
                            }
                        }
                        _ => None,
                    };
                    let (positions, sample_positions) = state
                        .first_no_match_positions
                        .as_ref()
                        .map(|v| (Some(v.2), Some(v.3.clone())))
                        .unwrap_or((None, None));
                    r.record(
                        &record.header,
                        reason,
                        kmer_pos,
                        min_pos,
                        kmer_seq,
                        min_seq,
                        positions,
                        sample_positions,
                        state.used_revcomp,
                    );
                } else {
                    r.record(
                        &record.header,
                        DebugFailReason::Unknown,
                        None,
                        None,
                        None,
                        None,
                        None,
                        None,
                        false,
                    );
                }
            }
            continue;
        };
        if let Some(filter) = filter {
            if !filter.single_overhang && filter.fragment_length > 0 {
                if let Some(best_match) = read_ec.best_match {
                    read_ec.ec = filter_ec_by_fragment(
                        index,
                        &read_ec.ec,
                        best_match,
                        filter.fragment_length as i64,
                    );
                }
            }
        }
        if let Some(mode) = options.strand_specific {
            let comprehensive = options.do_union || options.no_jump || index.use_shade;
            if let Some(filtered) = apply_strand_filter(
                index,
                &read_ec,
                mode,
                true,
                comprehensive,
                &mut report,
                &record.header,
            ) {
                read_ec.ec = filtered;
            }
        }
        if read_ec.ec.is_empty() {
            if let Some(r) = report.as_deref_mut() {
                if let Some(state) = dbg {
                    let (reason, kmer_pos, min_pos) = debug_reason(&state);
                    r.record(
                        &record.header,
                        reason,
                        kmer_pos,
                        min_pos,
                        None,
                        None,
                        None,
                        None,
                        state.used_revcomp,
                    );
                } else {
                    r.record(
                        &record.header,
                        DebugFailReason::Unknown,
                        None,
                        None,
                        None,
                        None,
                        None,
                        None,
                        false,
                    );
                }
            }
            continue;
        }

        if let Some(bias_counts) = bias.as_mut() {
            if bias_counts.total < options.max_bias as u64 {
                if let Some(best_match) = read_ec.best_match {
                    if let Some(hex) = bias_hexamer_for_match(index, best_match) {
                        bias_counts.record(hex);
                    }
                }
            }
        }

        reads_aligned += 1;
        let ec_id = match ec_map.get(&read_ec.ec) {
            Some(id) => *id,
            None => {
                let id = ec_list.len();
                ec_map.insert(read_ec.ec.clone(), id);
                ec_list.push(read_ec.ec.clone());
                counts.push(0);
                id
            }
        };
        counts[ec_id] += 1;
    }

    Ok(EcCounts {
        ec_list,
        counts,
        reads_processed,
        reads_aligned,
        bias,
    })
}

fn ec_for_read_naive(index: &KmerEcIndex, seq: &[u8]) -> Option<Vec<u32>> {
    if seq.len() < index.k {
        return None;
    }
    let mut current: Vec<u32> = Vec::new();
    let mut next: Vec<u32> = Vec::new();
    let mut has_hit = false;
    for pos in 0..=seq.len() - index.k {
        if let Some((fwd, rev)) = encode_kmer_pair(&seq[pos..pos + index.k]) {
            if let Some(ec) = lookup_ec(index, fwd).or_else(|| lookup_ec(index, rev)) {
                if !has_hit {
                    current.extend_from_slice(ec);
                    has_hit = true;
                } else {
                    next.clear();
                    intersect_sorted(&current, ec, &mut next);
                    std::mem::swap(&mut current, &mut next);
                    if current.is_empty() {
                        break;
                    }
                }
            }
        }
    }
    if !has_hit || current.is_empty() {
        return None;
    }
    if let Some(onlist) = index.onlist.as_deref() {
        current.retain(|&t| (t as usize) < onlist.len() && onlist[t as usize]);
    }
    if current.is_empty() {
        None
    } else {
        Some(current)
    }
}

fn ec_for_read_bifrost(
    index: &BifrostIndex,
    seq: &[u8],
    strand: Strand,
    mut dbg: Option<&mut ReadDebugState>,
    options: PseudoalignOptions,
) -> Option<ReadEc> {
    if seq.len() < index.k {
        if let Some(state) = dbg {
            state.saw_valid_kmer = false;
        }
        return None;
    }
    let mut current: Vec<u32> = Vec::new();
    let mut next: Vec<u32> = Vec::new();
    let mut rev_buf: Vec<u8> = Vec::new();
    let allow_forward = strand != Strand::Reverse;
    let allow_rev = strand != Strand::Forward;
    let mut has_hit = false;
    let diff = index.k.saturating_sub(index.g);
    let mut best_match: Option<MatchInfo> = None;
    let mut first_hit: Option<Hit> = None;
    let mut hits: Vec<Hit> = Vec::new();

    let mut pos = 0usize;
    let last_pos = seq.len() - index.k;
    while pos <= last_pos {
        let kmer = &seq[pos..pos + index.k];
        let mut matched = None;
        let min_candidates = match minimizers_for_kmer(kmer, index.g) {
            Some(v) => {
                if let Some(state) = dbg.as_deref_mut() {
                    state.saw_valid_kmer = true;
                }
                v
            }
            None => {
                pos += 1;
                continue;
            }
        };
        let mut any_mphf = false;
        let mut any_positions = false;
        for (min_bytes, min_pos) in min_candidates {
            let Some(min_idx) = index.mphf.lookup(&min_bytes) else {
                if let Some(state) = dbg.as_deref_mut() {
                    if state.first_mphf_miss.is_none() {
                        state.first_mphf_miss = Some((pos, min_pos));
                    }
                }
                continue;
            };
            any_mphf = true;
            let positions = &index.minz_positions[min_idx as usize];
            if positions.is_empty() {
                if let Some(state) = dbg.as_deref_mut() {
                    if state.first_no_positions.is_none() {
                        state.first_no_positions = Some((pos, min_pos));
                    }
                }
                continue;
            }
            any_positions = true;
            if let Some(state) = dbg.as_deref_mut() {
                if state.first_no_match_positions.is_none() {
                    state.first_no_match_positions = Some((
                        pos,
                        min_pos,
                        positions.len(),
                        format_positions_sample(positions),
                    ));
                }
            }
            let mut local_match = None;
            for &pos_id in positions {
                let (unitig_id, rel_pos, is_km) = decode_pos_id(pos_id);
                if is_km {
                    let uid = unitig_id as usize;
                    if uid >= index.km_unitigs.len() {
                        continue;
                    }
                    let km_pos = rel_pos as usize;
                    let rev_match = km_pos <= diff && min_pos == diff - km_pos;
                    let fwd_match = min_pos == km_pos;
                    let has_rev = rev_match && allow_rev;
                    let has_fwd = fwd_match && allow_forward;
                    if !has_rev && !has_fwd {
                        continue;
                    }
                    if !fill_revcomp(kmer, &mut rev_buf) {
                        continue;
                    }
                    let canonical = if rev_buf.as_slice() < kmer {
                        rev_buf.as_slice()
                    } else {
                        kmer
                    };
                    if index.km_unitigs[uid] == canonical {
                        local_match =
                            Some((index.unitigs.len() + uid, 0usize, pos, min_pos, rev_match));
                        if rev_match {
                            if let Some(state) = dbg.as_deref_mut() {
                                state.used_revcomp = true;
                            }
                        }
                        break;
                    }
                } else {
                    let uid = unitig_id as usize;
                    if uid >= index.unitigs.len() {
                        continue;
                    }
                    let rel_pos = rel_pos as isize;
                    if allow_forward {
                        let start = rel_pos - min_pos as isize;
                        if start >= 0 {
                            let start = start as usize;
                            if start + index.k <= index.unitigs[uid].len()
                                && &index.unitigs[uid][start..start + index.k] == kmer
                            {
                                local_match = Some((uid, start, pos, min_pos, false));
                                break;
                            }
                        }
                    }
                    if allow_rev && fill_revcomp(kmer, &mut rev_buf) {
                        let start = rel_pos - diff as isize + min_pos as isize;
                        if start >= 0 {
                            let start = start as usize;
                            if start + index.k <= index.unitigs[uid].len()
                                && index.unitigs[uid][start..start + index.k] == rev_buf
                            {
                                local_match = Some((uid, start, pos, min_pos, true));
                                if let Some(state) = dbg.as_deref_mut() {
                                    state.used_revcomp = true;
                                }
                                break;
                            }
                        }
                    }
                }
            }
            if let Some(hit) = local_match {
                matched = Some(hit);
                break;
            }
            if let Some(state) = dbg.as_deref_mut() {
                if state.first_no_match.is_none() {
                    state.first_no_match = Some((pos, min_pos));
                }
            }
        }
        if let Some(state) = dbg.as_deref_mut() {
            if any_mphf {
                state.saw_mphf_hit = true;
            }
            if any_positions {
                state.saw_positions = true;
            }
        }

        let Some((uid, start, kmer_pos, min_pos, used_revcomp)) = matched else {
            pos += 1;
            continue;
        };
        if let Some(state) = dbg.as_deref_mut() {
            state.saw_match = true;
        }
        if best_match.map(|m| kmer_pos < m.read_pos).unwrap_or(true) {
            best_match = Some(MatchInfo {
                unitig_id: uid,
                unitig_pos: start,
                read_pos: kmer_pos,
                used_revcomp,
            });
        }
        let Some(block_idx) = block_index_for_position(&index.ec_blocks[uid], start) else {
            pos += 1;
            continue;
        };
        let ec = &index.ec_blocks[uid][block_idx].ec;
        if ec.is_empty() {
            if let Some(state) = dbg.as_deref_mut() {
                if state.first_empty_ec.is_none() {
                    state.first_empty_ec = Some((kmer_pos, min_pos));
                }
            }
        } else if let Some(state) = dbg.as_deref_mut() {
            state.saw_ec = true;
        }
        hits.push(Hit {
            unitig_id: uid,
            read_pos: kmer_pos,
            block_idx,
            used_revcomp,
        });
        if first_hit
            .as_ref()
            .map(|h| kmer_pos < h.read_pos)
            .unwrap_or(true)
        {
            first_hit = Some(Hit {
                unitig_id: uid,
                read_pos: kmer_pos,
                block_idx,
                used_revcomp,
            });
        }
        has_hit = true;
        if !options.no_jump {
            if let Some(dist) = jump_distance_for_match(index, uid, block_idx, start, used_revcomp)
            {
                let mut next_pos = pos + dist;
                if next_pos > last_pos {
                    next_pos = last_pos;
                }
                if next_pos > pos {
                    let kmer_next = &seq[next_pos..next_pos + index.k];
                    if let Some((uid2, _start2, _rev2, block_idx2)) = match_kmer_at_pos(
                        index,
                        kmer_next,
                        allow_forward,
                        allow_rev,
                        diff,
                        &mut rev_buf,
                    ) {
                        if uid2 == uid
                            && index.ec_blocks[uid][block_idx].ec
                                == index.ec_blocks[uid2][block_idx2].ec
                        {
                            pos = next_pos;
                            continue;
                        }
                    }
                }
            }
        }
        pos += 1;
    }

    if !has_hit || hits.is_empty() {
        return None;
    }

    hits.sort_by(|a, b| {
        a.unitig_id
            .cmp(&b.unitig_id)
            .then_with(|| a.read_pos.cmp(&b.read_pos))
    });

    let mut minpos = usize::MAX;
    let mut maxpos = 0usize;
    let mut last_ec: Vec<u32> = Vec::new();
    let mut last_unitig = None;
    let mut found_nonempty = false;
    let mut current_offlist = false;
    let use_shade = index.use_shade && !options.do_union;
    let mut shade_scratch: Vec<u32> = Vec::new();
    let mut ec_scratch: Vec<u32> = Vec::new();
    for hit in &hits {
        minpos = minpos.min(hit.read_pos);
        maxpos = maxpos.max(hit.read_pos);
        let ec = &index.ec_blocks[hit.unitig_id][hit.block_idx].ec;
        let ec = if use_shade {
            filter_shades(ec, &index.shade_sequences, &mut shade_scratch);
            shade_scratch.as_slice()
        } else {
            ec
        };
        let ec_offlist = if let Some(onlist) = index.onlist.as_deref() {
            ec_has_offlist(ec, onlist)
        } else {
            false
        };
        if !found_nonempty {
            if !ec.is_empty() {
                if options.do_union {
                    current.extend_from_slice(ec);
                    current.sort_unstable();
                    current.dedup();
                } else {
                    current.extend_from_slice(ec);
                    last_ec = ec.to_vec();
                    last_unitig = Some(hit.unitig_id);
                }
                current_offlist = ec_offlist;
                found_nonempty = true;
            }
            continue;
        }
        if last_unitig == Some(hit.unitig_id) && ec == last_ec.as_slice() {
            continue;
        }
        if ec.is_empty() {
            continue;
        }
        if options.do_union {
            if options.dfk_onlist && (current_offlist || ec_offlist) {
                let dummy = index.onlist.as_ref().map(|v| v.len() as u32);
                if let Some(dummy) = dummy {
                    if current.last().copied() != Some(dummy) {
                        current.push(dummy);
                    }
                }
            }
            merge_sorted_unique_vec(&mut current, ec);
            current_offlist = current_offlist || ec_offlist;
        } else {
            if current.is_empty() {
                current.extend_from_slice(ec);
                last_ec = ec.to_vec();
                last_unitig = Some(hit.unitig_id);
                current_offlist = ec_offlist;
                continue;
            }
            let ec_slice = if options.dfk_onlist && (current_offlist || ec_offlist) {
                if let Some(onlist) = index.onlist.as_deref() {
                    let dummy = onlist.len() as u32;
                    ec_scratch.clear();
                    ec_scratch.extend_from_slice(ec);
                    if ec_scratch.last().copied() != Some(dummy) {
                        ec_scratch.push(dummy);
                    }
                    if current.last().copied() != Some(dummy) {
                        current.push(dummy);
                    }
                    ec_scratch.sort_unstable();
                    ec_scratch.dedup();
                    ec_scratch.as_slice()
                } else {
                    ec
                }
            } else {
                ec
            };
            next.clear();
            intersect_sorted(&current, ec_slice, &mut next);
            std::mem::swap(&mut current, &mut next);
            if current.is_empty() {
                if let Some(state) = dbg.as_deref_mut() {
                    state.intersection_empty = true;
                }
                return None;
            }
            last_ec = ec.to_vec();
            last_unitig = Some(hit.unitig_id);
            current_offlist = current_offlist || ec_offlist;
        }
    }

    if current.is_empty() {
        return None;
    }
    let min_range = options.min_range.max(1);
    if maxpos >= minpos && (maxpos - minpos + index.k) < min_range {
        return None;
    }
    if let Some(onlist) = index.onlist.as_deref() {
        current.retain(|&t| (t as usize) < onlist.len() && onlist[t as usize]);
        if options.dfk_onlist && current_offlist && !current.is_empty() {
            let dummy = onlist.len() as u32;
            if current.last().copied() != Some(dummy) {
                current.push(dummy);
            }
        }
    }
    if current.is_empty() {
        None
    } else {
        let mut shade_union = Vec::new();
        if index.use_shade {
            for hit in &hits {
                let ec = &index.ec_blocks[hit.unitig_id][hit.block_idx].ec;
                for &tr in ec {
                    let idx = tr as usize;
                    if idx < index.shade_sequences.len() && index.shade_sequences[idx] {
                        shade_union.push(tr);
                    }
                }
            }
            shade_union.sort_unstable();
            shade_union.dedup();
        }
        Some(ReadEc {
            ec: current,
            best_match,
            first_hit,
            hits,
            had_offlist: current_offlist,
            shade_union,
        })
    }
}

fn match_kmer_at_pos(
    index: &BifrostIndex,
    kmer: &[u8],
    allow_forward: bool,
    allow_rev: bool,
    diff: usize,
    rev_buf: &mut Vec<u8>,
) -> Option<(usize, usize, bool, usize)> {
    let min_candidates = minimizers_for_kmer(kmer, index.g)?;
    for (min_bytes, min_pos) in min_candidates {
        let Some(min_idx) = index.mphf.lookup(&min_bytes) else {
            continue;
        };
        let positions = &index.minz_positions[min_idx as usize];
        if positions.is_empty() {
            continue;
        }
        for &pos_id in positions {
            let (unitig_id, rel_pos, is_km) = decode_pos_id(pos_id);
            if is_km {
                let uid = unitig_id as usize;
                if uid >= index.km_unitigs.len() {
                    continue;
                }
                let km_pos = rel_pos as usize;
                let rev_match = km_pos <= diff && min_pos == diff - km_pos;
                let fwd_match = min_pos == km_pos;
                let has_rev = rev_match && allow_rev;
                let has_fwd = fwd_match && allow_forward;
                if !has_rev && !has_fwd {
                    continue;
                }
                if !fill_revcomp(kmer, rev_buf) {
                    continue;
                }
                let canonical = if rev_buf.as_slice() < kmer {
                    rev_buf.as_slice()
                } else {
                    kmer
                };
                if index.km_unitigs[uid] == canonical {
                    let unitig_id = index.unitigs.len() + uid;
                    let block_idx = block_index_for_position(&index.ec_blocks[unitig_id], 0)?;
                    return Some((unitig_id, 0usize, rev_match, block_idx));
                }
            } else {
                let uid = unitig_id as usize;
                if uid >= index.unitigs.len() {
                    continue;
                }
                let rel_pos = rel_pos as isize;
                if allow_forward {
                    let start = rel_pos - min_pos as isize;
                    if start >= 0 {
                        let start = start as usize;
                        if start + index.k <= index.unitigs[uid].len()
                            && &index.unitigs[uid][start..start + index.k] == kmer
                        {
                            let block_idx = block_index_for_position(&index.ec_blocks[uid], start)?;
                            return Some((uid, start, false, block_idx));
                        }
                    }
                }
                if allow_rev && fill_revcomp(kmer, rev_buf) {
                    let start = rel_pos - diff as isize + min_pos as isize;
                    if start >= 0 {
                        let start = start as usize;
                        if start + index.k <= index.unitigs[uid].len()
                            && &index.unitigs[uid][start..start + index.k] == rev_buf.as_slice()
                        {
                            let block_idx = block_index_for_position(&index.ec_blocks[uid], start)?;
                            return Some((uid, start, true, block_idx));
                        }
                    }
                }
            }
        }
    }
    None
}

fn jump_distance_for_match(
    index: &BifrostIndex,
    unitig_id: usize,
    block_idx: usize,
    start: usize,
    used_revcomp: bool,
) -> Option<usize> {
    let blocks = index.ec_blocks.get(unitig_id)?;
    let block = blocks.get(block_idx)?;
    if block.ub <= block.lb {
        return None;
    }
    let contig_start = block.lb as isize;
    let contig_len = (block.ub - block.lb) as isize;
    let start = start as isize;
    let forward = !used_revcomp;
    let dist = if forward {
        contig_len - 1 - (start - contig_start)
    } else {
        start - contig_start
    };
    if dist >= 2 {
        Some(dist as usize)
    } else {
        None
    }
}

fn ec_has_offlist(ec: &[u32], onlist: &[bool]) -> bool {
    for &t in ec {
        let idx = t as usize;
        if idx >= onlist.len() || !onlist[idx] {
            return true;
        }
    }
    false
}

fn filter_shades(ec: &[u32], shade_sequences: &[bool], out: &mut Vec<u32>) {
    out.clear();
    for &tr in ec {
        let idx = tr as usize;
        if idx < shade_sequences.len() && shade_sequences[idx] {
            continue;
        }
        out.push(tr);
    }
}

fn filter_shades_in_place(ec: &mut Vec<u32>, shade_sequences: &[bool]) {
    ec.retain(|&tr| {
        let idx = tr as usize;
        idx >= shade_sequences.len() || !shade_sequences[idx]
    });
}

fn merge_sorted_unique_vec(base: &mut Vec<u32>, extra: &[u32]) {
    if extra.is_empty() {
        return;
    }
    if base.is_empty() {
        base.extend_from_slice(extra);
        return;
    }
    let mut merged = Vec::with_capacity(base.len() + extra.len());
    let mut i = 0;
    let mut j = 0;
    while i < base.len() && j < extra.len() {
        if base[i] == extra[j] {
            merged.push(base[i]);
            i += 1;
            j += 1;
        } else if base[i] < extra[j] {
            merged.push(base[i]);
            i += 1;
        } else {
            merged.push(extra[j]);
            j += 1;
        }
    }
    if i < base.len() {
        merged.extend_from_slice(&base[i..]);
    }
    if j < extra.len() {
        merged.extend_from_slice(&extra[j..]);
    }
    merged.dedup();
    *base = merged;
}

fn debug_reason(state: &ReadDebugState) -> (DebugFailReason, Option<usize>, Option<usize>) {
    if !state.saw_valid_kmer {
        return (DebugFailReason::NoValidKmer, None, None);
    }
    if !state.saw_mphf_hit {
        let loc = state.first_mphf_miss;
        return (
            DebugFailReason::NoMinimizerHit,
            loc.map(|v| v.0),
            loc.map(|v| v.1),
        );
    }
    if !state.saw_positions {
        let loc = state.first_no_positions;
        return (
            DebugFailReason::NoPositions,
            loc.map(|v| v.0),
            loc.map(|v| v.1),
        );
    }
    if !state.saw_match {
        let loc = state.first_no_match;
        return (
            DebugFailReason::NoKmerMatch,
            loc.map(|v| v.0),
            loc.map(|v| v.1),
        );
    }
    if !state.saw_ec {
        let loc = state.first_empty_ec;
        return (DebugFailReason::EmptyEc, loc.map(|v| v.0), loc.map(|v| v.1));
    }
    if state.intersection_empty {
        return (DebugFailReason::IntersectionEmpty, None, None);
    }
    (DebugFailReason::Unknown, None, None)
}

fn intersect_sorted(a: &[u32], b: &[u32], out: &mut Vec<u32>) {
    out.clear();
    let mut i = 0;
    let mut j = 0;
    while i < a.len() && j < b.len() {
        if a[i] == b[j] {
            out.push(a[i]);
            i += 1;
            j += 1;
        } else if a[i] < b[j] {
            i += 1;
        } else {
            j += 1;
        }
    }
}

fn merge_sorted_unique(a: &mut Vec<u32>, b: &[u32]) {
    if a.is_empty() {
        a.extend_from_slice(b);
        return;
    }
    let mut merged = Vec::with_capacity(a.len() + b.len());
    let mut i = 0;
    let mut j = 0;
    while i < a.len() && j < b.len() {
        if a[i] == b[j] {
            merged.push(a[i]);
            i += 1;
            j += 1;
        } else if a[i] < b[j] {
            merged.push(a[i]);
            i += 1;
        } else {
            merged.push(b[j]);
            j += 1;
        }
    }
    if i < a.len() {
        merged.extend_from_slice(&a[i..]);
    }
    if j < b.len() {
        merged.extend_from_slice(&b[j..]);
    }
    merged.dedup();
    *a = merged;
}

fn minimizers_for_kmer(seq: &[u8], g: usize) -> Option<Vec<([u8; 8], usize)>> {
    if seq.len() < g + 2 {
        return None;
    }
    let mut best_hash: Option<u64> = None;
    let mut best: Vec<([u8; 8], usize)> = Vec::new();
    let start = 1usize;
    let end = seq.len().saturating_sub(g + 1);
    for pos in start..=end {
        let slice = &seq[pos..pos + g];
        let h = rep_hash(slice)?;
        let bytes = encode_minimizer_rep(slice)?;
        match best_hash {
            None => {
                best_hash = Some(h);
                best.push((bytes, pos));
            }
            Some(min) if h < min => {
                best_hash = Some(h);
                best.clear();
                best.push((bytes, pos));
            }
            Some(min) if h == min => {
                best.push((bytes, pos));
            }
            _ => {}
        }
    }
    if best.is_empty() {
        None
    } else {
        Some(best)
    }
}

fn rep_hash(seq: &[u8]) -> Option<u64> {
    let mut h = 0u64;
    let mut ht = 0u64;
    let len = seq.len();
    for i in 0..len {
        let b = seq[i];
        if !matches!(b, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't') {
            return None;
        }
        h = h.rotate_left(1);
        ht = ht.rotate_left(1);
        let f = ((b & 6) >> 1) as usize;
        let rb = seq[len - 1 - i];
        let r = (((rb ^ 4) & 6) >> 1) as usize;
        h ^= REP_HASH_VALS[f];
        ht ^= REP_HASH_VALS[r];
    }
    let mut hashes = [h, ht];
    if hashes[1] < hashes[0] {
        hashes.swap(0, 1);
    }
    let mut buf = [0u8; 16];
    buf[..8].copy_from_slice(&hashes[0].to_le_bytes());
    buf[8..].copy_from_slice(&hashes[1].to_le_bytes());
    Some(wyhash(&buf, 0))
}

fn fill_revcomp(seq: &[u8], out: &mut Vec<u8>) -> bool {
    out.clear();
    out.reserve(seq.len());
    for &b in seq.iter().rev() {
        let comp = match b {
            b'A' | b'a' => b'T',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            b'T' | b't' => b'A',
            _ => return false,
        };
        out.push(comp);
    }
    true
}

fn format_positions_sample(positions: &[u64]) -> String {
    let mut out = String::new();
    let take = positions.len().min(3);
    for (idx, &pos_id) in positions.iter().take(take).enumerate() {
        let (unitig, rel_pos, is_km) = decode_pos_id(pos_id);
        if idx > 0 {
            out.push(',');
        }
        let kind = if is_km { 'k' } else { 'u' };
        out.push_str(&format!("{}:{}:{}", unitig, rel_pos, kind));
    }
    out
}

struct OnlistData {
    onlist: Option<Vec<bool>>,
    lengths: Vec<u32>,
    names: Vec<String>,
    shade_to_color: Vec<Option<u32>>,
    shade_sequences: Vec<bool>,
    use_shade: bool,
}

fn split_shade_name(name: &str) -> Option<(&str, &str)> {
    let needle = "_shade_";
    name.find(needle)
        .map(|idx| (&name[..idx], &name[idx + needle.len()..]))
}

fn read_onlist<R: Read>(reader: &mut R) -> Result<Option<Vec<bool>>> {
    read_onlist_with_lengths(reader).map(|data| data.onlist)
}

fn read_onlist_with_lengths<R: Read>(reader: &mut R) -> Result<OnlistData> {
    let num_trans = match read_i32_le(reader) {
        Ok(v) => v,
        Err(_) => {
            return Ok(OnlistData {
                onlist: None,
                lengths: Vec::new(),
                names: Vec::new(),
                shade_to_color: Vec::new(),
                shade_sequences: Vec::new(),
                use_shade: false,
            })
        }
    };
    if num_trans < 0 {
        return Err(Error::InvalidFormat("negative transcript count".into()));
    }
    let num_trans = num_trans as usize;

    let mut lengths = Vec::with_capacity(num_trans);
    for _ in 0..num_trans {
        let len = read_i32_le(reader)?;
        if len < 0 {
            return Err(Error::InvalidFormat("negative transcript length".into()));
        }
        lengths.push(len as u32);
    }
    let mut names = Vec::with_capacity(num_trans);
    let mut shade_to_color = vec![None; num_trans];
    let mut shade_sequences = vec![false; num_trans];
    let mut name_to_index: HashMap<String, usize> = HashMap::new();
    let mut use_shade = false;
    for idx in 0..num_trans {
        let name_len = read_u64_le(reader)? as usize;
        let mut buf = vec![0u8; name_len];
        reader.read_exact(&mut buf)?;
        let name = String::from_utf8_lossy(&buf).into_owned();
        if let Some((base, _variant)) = split_shade_name(&name) {
            use_shade = true;
            shade_sequences[idx] = true;
            if let Some(color_idx) = name_to_index.get(base) {
                shade_to_color[idx] = Some(*color_idx as u32);
            }
        } else {
            name_to_index.entry(name.clone()).or_insert(idx);
        }
        names.push(name);
    }

    let onlist_size = read_u64_le(reader)?;
    if onlist_size == 0 {
        return Ok(OnlistData {
            onlist: None,
            lengths,
            names,
            shade_to_color,
            shade_sequences,
            use_shade,
        });
    }
    let mut buf = vec![0u8; onlist_size as usize];
    reader.read_exact(&mut buf)?;
    let vals = unsafe { deserialize_roaring_to_vec(&buf) }
        .ok_or_else(|| Error::InvalidFormat("invalid onlist bitmap".into()))?;
    let mut onlist = vec![false; num_trans];
    for v in vals {
        let idx = v as usize;
        if idx < onlist.len() {
            onlist[idx] = true;
        }
    }
    Ok(OnlistData {
        onlist: Some(onlist),
        lengths,
        names,
        shade_to_color,
        shade_sequences,
        use_shade,
    })
}

fn filter_ec_by_fragment(
    index: &BifrostIndex,
    ec: &[u32],
    best_match: MatchInfo,
    fragment_length: i64,
) -> Vec<u32> {
    if fragment_length <= 0 {
        return ec.to_vec();
    }
    if index.transcript_lengths.is_empty() {
        return ec.to_vec();
    }
    if best_match.unitig_id >= index.ec_blocks.len() {
        return ec.to_vec();
    }
    let blocks = &index.ec_blocks[best_match.unitig_id];
    if blocks.is_empty() || blocks[0].positions.is_none() {
        return ec.to_vec();
    }
    let unitig_len = if best_match.unitig_id < index.unitigs.len() {
        index.unitigs[best_match.unitig_id].len()
    } else {
        index.k
    };
    let mut filtered = Vec::new();
    for &tr in ec {
        let (pos, forward) = match find_position_in_transcript(
            blocks,
            tr,
            best_match.unitig_pos,
            best_match.read_pos,
            best_match.used_revcomp,
            unitig_len,
            index.k,
        ) {
            Some((pos, forward)) => (pos, forward),
            None => {
                filtered.push(tr);
                continue;
            }
        };
        let len = index
            .transcript_lengths
            .get(tr as usize)
            .copied()
            .unwrap_or(0) as i64;
        if forward {
            if pos + fragment_length - 1 <= len {
                filtered.push(tr);
            }
        } else if pos - fragment_length >= 0 {
            filtered.push(tr);
        }
    }
    filtered
}

fn bias_hexamer_for_match(index: &BifrostIndex, best_match: MatchInfo) -> Option<usize> {
    let blocks = index.ec_blocks.get(best_match.unitig_id)?;
    let (lb, ub) = ec_block_at(blocks, best_match.unitig_pos as u32)?;
    let contig_start = lb as isize;
    let contig_len = (ub as isize) - contig_start;
    let pos = best_match.unitig_pos as isize - contig_start;
    let p = best_match.read_pos as isize;
    let pre = 2isize;
    let post = 4isize;
    let unitig = index.unitigs.get(best_match.unitig_id)?;
    if unitig.len() < 6 {
        return None;
    }
    let k = index.k as isize;
    let um_strand = !best_match.used_revcomp;

    if um_strand {
        if pos - p < pre {
            return None;
        }
        let start = contig_start + pos - p - pre;
        if start < 0 {
            return None;
        }
        let start = start as usize;
        if start + 6 > unitig.len() {
            return None;
        }
        hexamer_to_int(&unitig[start..start + 6], true)
    } else {
        if contig_len - 1 - pos - p < pre {
            return None;
        }
        let pos_ = (pos + p) + k - post;
        let start = contig_start + pos_;
        if start < 0 {
            return None;
        }
        let start = start as usize;
        if start + 6 > unitig.len() {
            return None;
        }
        hexamer_to_int(&unitig[start..start + 6], false)
    }
}

fn apply_strand_filter(
    index: &BifrostIndex,
    read_ec: &ReadEc,
    mode: StrandSpecific,
    is_first_read: bool,
    comprehensive: bool,
    report: &mut Option<&mut DebugReport>,
    header: &[u8],
) -> Option<Vec<u32>> {
    let ec = &read_ec.ec;
    let target = match mode {
        StrandSpecific::FR => is_first_read,
        StrandSpecific::RF => !is_first_read,
    };

    if comprehensive {
        let mut union = Vec::new();
        for hit in &read_ec.hits {
            if let Some(filtered) = filter_ec_for_hit(index, ec, hit, target) {
                if !filtered.is_empty() {
                    merge_sorted_unique_vec(&mut union, &filtered);
                }
            }
        }
        return Some(union);
    }

    let hit = read_ec.first_hit.as_ref()?;
    let filtered = filter_ec_for_hit(index, ec, hit, target)?;
    if filtered.len() < ec.len() {
        if let Some(r) = report.as_deref_mut() {
            r.record(
                header,
                DebugFailReason::Unknown,
                None,
                None,
                None,
                None,
                None,
                None,
                read_ec
                    .best_match
                    .as_ref()
                    .map(|m| m.used_revcomp)
                    .unwrap_or(false),
            );
        }
    }
    Some(filtered)
}

fn filter_ec_for_hit(
    index: &BifrostIndex,
    ec: &[u32],
    hit: &Hit,
    target: bool,
) -> Option<Vec<u32>> {
    let block = index.ec_blocks.get(hit.unitig_id)?;
    let block = block.get(hit.block_idx)?;
    let strands = block.strands.as_ref()?;
    let mut filtered = Vec::new();
    let um_strand = !hit.used_revcomp;
    for &tr in ec {
        let idx = match block.ec.binary_search(&tr) {
            Ok(v) => v,
            Err(_) => continue,
        };
        let sense = strands.get(idx).copied().unwrap_or(2);
        if sense == 2 || ((um_strand == (sense == 1)) == target) {
            filtered.push(tr);
        }
    }
    Some(filtered)
}

fn find_position_in_transcript(
    blocks: &[crate::index::EcBlock],
    tr: u32,
    unitig_pos: usize,
    read_pos: usize,
    used_revcomp: bool,
    unitig_len: usize,
    k: usize,
) -> Option<(i64, bool)> {
    let idx = unitig_pos as u32;
    let mc = ec_block_at(blocks, idx)?;
    let ecs = ec_blocks_leading_vals(blocks, idx);
    if ecs.is_empty() {
        return None;
    }
    let v_ec = ecs.last().copied()?;
    let rawpos = block_min_pos(v_ec, tr)?;
    let trpos = (rawpos & 0x7fff_ffff) as i64;
    let trsense = rawpos == (rawpos & 0x7fff_ffff);

    let csense = !used_revcomp;
    let um_dist = unitig_pos as i64;
    let um_size = unitig_len as i64;
    let p = read_pos as i64;
    let k = k as i64;
    let mc_first = mc.0 as i64;
    let mc_second = mc.1 as i64;
    let _um_dist_block = um_dist - mc_first;

    if trsense {
        if csense {
            let mut padding = 0i64;
            if trpos == 0 {
                let mut mc_cur = mc;
                for block in ecs.iter().rev().skip(1) {
                    if !block_contains(block, tr) {
                        padding = mc_cur.0 as i64;
                        break;
                    }
                    mc_cur = prev_block_at(blocks, mc_cur.0)?;
                }
            }
            let pos = trpos - p + um_dist + 1 - padding;
            Some((pos, csense))
        } else {
            let mut mc_cur = mc;
            let mut right_one = 0i64;
            let mut left_one = 0i64;
            let initial = mc_second;
            for (i, block) in ecs.iter().enumerate().rev() {
                if i == ecs.len() - 1 {
                    right_one = mc_cur.1 as i64;
                }
                if !block_contains(block, tr) {
                    left_one = mc_cur.1 as i64;
                    break;
                } else if i == 0 {
                    left_one = 0;
                }
                mc_cur = prev_block_at(blocks, mc_cur.0)?;
            }
            let padding = -(left_one + right_one - um_size + k - 1);
            let pos = trpos + p + k - (um_size - k - um_dist) + initial - 1 + padding;
            Some((pos, csense))
        }
    } else if csense {
        let mut left_one = 0i64;
        let mut right_one = 0i64;
        let mut unmapped_len = 0i64;
        let mut found_first_mapped = false;
        let ecs_all = ec_blocks_leading_vals(blocks, u32::MAX);
        let mut curr_mc = 0u32;
        for block in ecs_all {
            let mc_ = ec_block_at(blocks, curr_mc)?;
            if !block_contains(block, tr) && found_first_mapped {
                if unmapped_len == 0 {
                    left_one = mc_.0 as i64;
                }
                right_one = mc_.1 as i64;
                unmapped_len += mc_.1 as i64 - mc_.0 as i64;
            }
            if block_contains(block, tr) {
                found_first_mapped = true;
            }
            curr_mc = mc_.1;
        }
        let mut start = 0i64;
        start -= right_one - left_one;
        start += um_size - k;
        let pos = trpos + (-(um_dist - start)) + k + p;
        Some((pos, !csense))
    } else {
        let mut left_one = 0i64;
        let mut right_one = 0i64;
        let mut unmapped_len = 0i64;
        let mut found_first_mapped = false;
        let ecs_all = ec_blocks_leading_vals(blocks, u32::MAX);
        let mut curr_mc = 0u32;
        for block in ecs_all {
            let mc_ = ec_block_at(blocks, curr_mc)?;
            if !block_contains(block, tr) && found_first_mapped {
                if unmapped_len == 0 {
                    left_one = mc_.0 as i64;
                }
                right_one = mc_.1 as i64;
                unmapped_len += mc_.1 as i64 - mc_.0 as i64;
            }
            if block_contains(block, tr) {
                found_first_mapped = true;
            }
            curr_mc = mc_.1;
        }
        let unmapped_len = right_one - left_one;
        let padding = um_size - um_dist - unmapped_len - k + 1;
        let pos = trpos + padding - p;
        Some((pos, !csense))
    }
}

fn ec_blocks_leading_vals(
    blocks: &[crate::index::EcBlock],
    idx: u32,
) -> Vec<&crate::index::EcBlock> {
    if blocks.is_empty() {
        return Vec::new();
    }
    let mut lo = 0usize;
    let mut hi = blocks.len();
    while lo < hi {
        let mid = (lo + hi) / 2;
        if blocks[mid].lb <= idx {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    blocks[..lo].iter().collect()
}

fn block_index_for_position(blocks: &[crate::index::EcBlock], pos: usize) -> Option<usize> {
    if blocks.is_empty() {
        return None;
    }
    if blocks.len() == 1 {
        return Some(0);
    }
    let pos = pos as u32;
    let mut lo = 0usize;
    let mut hi = blocks.len();
    while lo < hi {
        let mid = (lo + hi) / 2;
        if blocks[mid].lb <= pos {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    if lo == 0 {
        return None;
    }
    Some(lo - 1)
}

fn ec_block_at(blocks: &[crate::index::EcBlock], idx: u32) -> Option<(u32, u32)> {
    if blocks.is_empty() {
        return None;
    }
    if blocks.len() == 1 {
        return Some((blocks[0].lb, blocks[0].ub));
    }
    let mut lo = 0usize;
    let mut hi = blocks.len();
    while lo < hi {
        let mid = (lo + hi) / 2;
        if blocks[mid].lb <= idx {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    if lo == 0 {
        return None;
    }
    let block = &blocks[lo - 1];
    Some((block.lb, block.ub))
}

fn prev_block_at(blocks: &[crate::index::EcBlock], lb: u32) -> Option<(u32, u32)> {
    if lb == 0 {
        ec_block_at(blocks, u32::MAX)
    } else {
        ec_block_at(blocks, lb - 1)
    }
}

fn block_contains(block: &crate::index::EcBlock, tr: u32) -> bool {
    block.ec.binary_search(&tr).is_ok()
}

fn block_min_pos(block: &crate::index::EcBlock, tr: u32) -> Option<u32> {
    let positions = block.positions.as_ref()?;
    let idx = block.ec.binary_search(&tr).ok()?;
    positions.get(idx).copied()
}

fn encode_kmer_pair(seq: &[u8]) -> Option<(u64, u64)> {
    let mut fwd = 0u64;
    for &b in seq {
        let v = match b {
            b'A' | b'a' => 0u64,
            b'C' | b'c' => 1u64,
            b'G' | b'g' => 2u64,
            b'T' | b't' => 3u64,
            _ => return None,
        };
        fwd = (fwd << 2) | v;
    }
    let mut rev = 0u64;
    for &b in seq.iter().rev() {
        let v = match b {
            b'A' | b'a' => 0u64,
            b'C' | b'c' => 1u64,
            b'G' | b'g' => 2u64,
            b'T' | b't' => 3u64,
            _ => return None,
        };
        let comp = 3 - v;
        rev = (rev << 2) | comp;
    }
    Some((fwd, rev))
}

fn lookup_ec(index: &KmerEcIndex, code: u64) -> Option<&[u32]> {
    let idx = index.mphf.try_hash(&code)? as usize;
    if index.keys_by_index.get(idx) == Some(&code) {
        return Some(&index.ecs[idx]);
    }
    None
}

fn read_unitig_sequences<R: Read + Seek>(
    reader: &mut R,
    dbg_start: u64,
    dbg_size: u64,
    k: u32,
    kmer_bytes: usize,
) -> Result<Vec<Vec<u8>>> {
    reader.seek(SeekFrom::Start(dbg_start))?;
    let _file_format_version = read_u64_le(reader)?;
    let _k = read_i32_le(reader)?;
    let _g = read_i32_le(reader)?;

    let v_unitigs_sz = read_u64_le(reader)?;
    let mut seqs: Vec<Vec<u8>> = Vec::with_capacity(v_unitigs_sz as usize);
    for _ in 0..v_unitigs_sz {
        let len = read_u64_le(reader)? as usize;
        let data_sz = len.div_ceil(4);
        let mut data = vec![0u8; data_sz];
        reader.read_exact(&mut data)?;
        seqs.push(decode_compressed_sequence(&data, len));
    }

    let km_unitigs_sz = read_u64_le(reader)?;
    for _ in 0..km_unitigs_sz {
        let mut buf = vec![0u8; kmer_bytes];
        reader.read_exact(&mut buf)?;
        seqs.push(decode_kmer(&buf, k as usize));
    }

    let h_kmers_ccov_sz = read_u64_le(reader)?;
    for _ in 0..h_kmers_ccov_sz {
        let mut buf = vec![0u8; kmer_bytes];
        reader.read_exact(&mut buf)?;
        seqs.push(decode_kmer(&buf, k as usize));
    }

    reader.seek(SeekFrom::Start(dbg_start + dbg_size))?;
    Ok(seqs)
}

fn decode_compressed_sequence(data: &[u8], len: usize) -> Vec<u8> {
    let mut out = Vec::with_capacity(len);
    for i in 0..len {
        let byte = data[i >> 2];
        let shift = (i & 0x3) << 1;
        let v = (byte >> shift) & 0x03;
        out.push(match v {
            0 => b'A',
            1 => b'C',
            2 => b'G',
            _ => b'T',
        });
    }
    out
}

fn decode_kmer(bytes: &[u8], k: usize) -> Vec<u8> {
    let mut out = Vec::with_capacity(k);
    let mut longs: Vec<u64> = Vec::with_capacity(bytes.len() / 8);
    for chunk in bytes.chunks_exact(8) {
        let mut arr = [0u8; 8];
        arr.copy_from_slice(chunk);
        longs.push(u64::from_le_bytes(arr));
    }
    for i in 0..k {
        let idx = i >> 5;
        let shift = 62 - ((i & 0x1F) << 1);
        let v = (longs[idx] >> shift) & 0x3;
        out.push(match v {
            0 => b'A',
            1 => b'C',
            2 => b'G',
            _ => b'T',
        });
    }
    out
}

fn read_u64_le<R: Read>(reader: &mut R) -> Result<u64> {
    let mut buf = [0u8; 8];
    reader.read_exact(&mut buf)?;
    Ok(u64::from_le_bytes(buf))
}

fn read_u32_le<R: Read>(reader: &mut R) -> Result<u32> {
    let mut buf = [0u8; 4];
    reader.read_exact(&mut buf)?;
    Ok(u32::from_le_bytes(buf))
}

fn read_i32_le<R: Read>(reader: &mut R) -> Result<i32> {
    let mut buf = [0u8; 4];
    reader.read_exact(&mut buf)?;
    Ok(i32::from_le_bytes(buf))
}
