use crate::bias::hexamer_to_int;

use super::{BifrostIndex, MatchInfo};

pub(super) fn bias_hexamer_for_match(index: &BifrostIndex, best_match: MatchInfo) -> Option<usize> {
    let block_idx = index
        .flat_ec
        .block_index_for_position(best_match.unitig_id, best_match.unitig_pos)?;
    let (lb, ub) = index
        .flat_ec
        .block_bounds(best_match.unitig_id, block_idx)?;
    let contig_start = lb as isize;
    let contig_len = (ub as isize) - contig_start;
    let pos = best_match.unitig_pos as isize - contig_start;
    let p = best_match.read_pos as isize;
    // The bias window is anchored at the inferred read start, two bases before
    // it and four from it onward. Reverse matches require the opposite window
    // orientation; the matched k-mer's start alone is not the sampling position.
    let pre = 2isize;
    let post = 4isize;
    let unitig = index.unitig_sequence(best_match.unitig_id)?;
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
