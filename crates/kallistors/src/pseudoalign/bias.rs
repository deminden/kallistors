use crate::bias::hexamer_to_int;

use super::{BifrostIndex, MatchInfo, ec_block_at};

pub(super) fn bias_hexamer_for_match(index: &BifrostIndex, best_match: MatchInfo) -> Option<usize> {
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
