use super::KmerEcIndex;

pub(crate) fn ec_has_offlist(ec: &[u32], onlist: &[bool]) -> bool {
    for &t in ec {
        let idx = t as usize;
        if idx >= onlist.len() || !onlist[idx] {
            return true;
        }
    }
    false
}

pub(crate) fn filter_shades(ec: &[u32], shade_sequences: &[bool], out: &mut Vec<u32>) {
    out.clear();
    for &tr in ec {
        let idx = tr as usize;
        if idx < shade_sequences.len() && shade_sequences[idx] {
            continue;
        }
        out.push(tr);
    }
}

pub(crate) fn filter_shades_in_place(ec: &mut Vec<u32>, shade_sequences: &[bool]) {
    ec.retain(|&tr| {
        let idx = tr as usize;
        idx >= shade_sequences.len() || !shade_sequences[idx]
    });
}

pub(crate) fn merge_sorted_unique_vec(base: &mut Vec<u32>, extra: &[u32]) {
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

pub(crate) fn intersect_sorted(a: &[u32], b: &[u32], out: &mut Vec<u32>) {
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

pub(crate) fn merge_sorted_unique(a: &mut Vec<u32>, b: &[u32]) {
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

pub(crate) fn encode_kmer_pair(seq: &[u8]) -> Option<(u64, u64)> {
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

pub(crate) fn lookup_ec(index: &KmerEcIndex, code: u64) -> Option<&[u32]> {
    let idx = index.mphf.try_hash(&code)? as usize;
    if index.keys_by_index.get(idx) == Some(&code) {
        return Some(&index.ecs[idx]);
    }
    None
}
