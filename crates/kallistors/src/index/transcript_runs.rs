//! Exact transcript-run scratch and stable per-unitig assembly.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub(super) struct TRInfo {
    pub(super) trid: u32,
    // Transcript offset with reverse orientation in the high bit. start/stop
    // delimit the covered unitig k-mer starts, using half-open coordinates.
    pub(super) pos: u32,
    pub(super) start: u32,
    pub(super) stop: u32,
}

pub(super) struct TrInfoShard {
    words: Vec<u32>,
    compact: bool,
}
impl TrInfoShard {
    pub(super) fn new(compact: bool) -> Self {
        Self {
            words: Vec::new(),
            compact,
        }
    }
    #[inline]
    pub(super) fn push(&mut self, unitig: usize, info: TRInfo, kmer_count: u32) {
        debug_assert!(u32::try_from(unitig).is_ok());
        let unitig = unitig as u32;
        if self.compact && info.start == 0 && info.stop == kmer_count {
            // Whole-unitig bounds are exact and already stored with the unitig.
            self.words
                .extend_from_slice(&[unitig | (1 << 31), info.trid, info.pos]);
        } else {
            self.words
                .extend_from_slice(&[unitig, info.trid, info.pos, info.start, info.stop]);
        }
    }
    #[inline]
    pub(super) fn iter<'a>(
        &'a self,
        counts: &'a [u32],
    ) -> impl Iterator<Item = (usize, TRInfo)> + 'a {
        let mut words = self.words.as_slice();
        std::iter::from_fn(move || {
            let &first = words.first()?;
            let whole = self.compact && first & (1 << 31) != 0;
            let unitig = (if whole { first & !(1 << 31) } else { first }) as usize;
            let info = TRInfo {
                trid: words[1],
                pos: words[2],
                start: if whole { 0 } else { words[3] },
                stop: if whole { counts[unitig] } else { words[4] },
            };
            words = &words[if whole { 3 } else { 5 }..];
            Some((unitig, info))
        })
    }
}

pub(super) struct TranscriptLocations {
    pub(super) offsets: Vec<usize>,
    pub(super) records: Vec<TRInfo>,
}
impl TranscriptLocations {
    pub(super) fn from_shards(shards: Vec<TrInfoShard>, counts: &[u32]) -> Self {
        let unitigs = counts.len();
        let mut offsets = vec![0usize; unitigs + 1];
        for shard in &shards {
            for (unitig, _) in shard.iter(counts) {
                offsets[unitig + 1] += 1;
            }
        }
        for i in 1..offsets.len() {
            offsets[i] += offsets[i - 1];
        }
        let mut cursor = offsets[..unitigs].to_vec();
        let mut records = vec![TRInfo::default(); offsets[unitigs]];
        for shard in shards {
            for (unitig, info) in shard.iter(counts) {
                records[cursor[unitig]] = info;
                cursor[unitig] += 1;
            }
        }
        Self { offsets, records }
    }
}
