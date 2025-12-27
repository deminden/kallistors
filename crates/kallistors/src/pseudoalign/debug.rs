use std::collections::HashMap;

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

impl DebugReport {
    pub(crate) fn new(max_traces: usize) -> Self {
        Self {
            counts: HashMap::new(),
            traces: Vec::new(),
            max_traces,
        }
    }

    #[allow(clippy::too_many_arguments)]
    pub(crate) fn record(
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
pub(crate) struct ReadDebugState {
    pub(crate) saw_valid_kmer: bool,
    pub(crate) saw_mphf_hit: bool,
    pub(crate) saw_positions: bool,
    pub(crate) saw_match: bool,
    pub(crate) saw_ec: bool,
    pub(crate) intersection_empty: bool,
    pub(crate) first_mphf_miss: Option<(usize, usize)>,
    pub(crate) first_no_positions: Option<(usize, usize)>,
    pub(crate) first_no_match: Option<(usize, usize)>,
    pub(crate) first_empty_ec: Option<(usize, usize)>,
    pub(crate) first_no_match_positions: Option<(usize, usize, usize, String)>,
    pub(crate) used_revcomp: bool,
}

pub(crate) fn debug_reason(
    state: &ReadDebugState,
) -> (DebugFailReason, Option<usize>, Option<usize>) {
    if !state.saw_valid_kmer {
        return (DebugFailReason::NoValidKmer, None, None);
    }
    if !state.saw_mphf_hit {
        return (
            DebugFailReason::NoMinimizerHit,
            state.first_mphf_miss.map(|v| v.0),
            state.first_mphf_miss.map(|v| v.1),
        );
    }
    if !state.saw_positions {
        return (
            DebugFailReason::NoPositions,
            state.first_no_positions.map(|v| v.0),
            state.first_no_positions.map(|v| v.1),
        );
    }
    if !state.saw_match {
        return (
            DebugFailReason::NoKmerMatch,
            state.first_no_match.map(|v| v.0),
            state.first_no_match.map(|v| v.1),
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

pub(crate) fn format_positions_sample(positions: &[u64]) -> String {
    let mut out = String::new();
    let take = positions.len().min(3);
    for (idx, &pos_id) in positions.iter().take(take).enumerate() {
        let (unitig, rel_pos, is_km) = crate::index::bifrost::decode_pos_id(pos_id);
        if idx > 0 {
            out.push(',');
        }
        let kind = if is_km { 'k' } else { 'u' };
        out.push_str(&format!("{}:{}:{}", unitig, rel_pos, kind));
    }
    out
}
