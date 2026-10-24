//! Native BooPHF construction and serialization, with stable level order.

use super::{write_f64, write_i32, write_u64};
use crate::index::bifrost::{ReducedDomain, minimizer_hash};
use crate::index::hashing::xorshift_next;
use crate::index::workers::IndexWorkers;
use crate::{Error, Result};
use rayon::prelude::*;
use std::io::Write;

const BOO_GAMMA: f64 = 2.0;
const BOO_LEVELS: usize = 25;

struct BooBitVector {
    size: u64,
    words: Vec<u64>,
    ranks: Vec<u64>,
}

impl BooBitVector {
    fn new(size: u64) -> Result<Self> {
        let nchar = size
            .checked_div(64)
            .and_then(|v| v.checked_add(1))
            .ok_or_else(|| Error::InvalidFormat("MPHF bitvector too large".into()))?;
        let nchar = usize::try_from(nchar)
            .map_err(|_| Error::InvalidFormat("MPHF bitvector too large".into()))?;
        Ok(Self {
            size,
            words: vec![0; nchar],
            ranks: Vec::new(),
        })
    }

    fn set(&mut self, pos: u64) {
        self.words[(pos >> 6) as usize] |= 1u64 << (pos & 63);
    }

    fn build_ranks(&mut self, offset: u64) -> u64 {
        self.ranks.clear();
        let mut rank = offset;
        for (idx, &word) in self.words.iter().enumerate() {
            if (idx as u64 * 64).is_multiple_of(512) {
                self.ranks.push(rank);
            }
            rank += u64::from(word.count_ones());
        }
        rank
    }

    fn write<W: Write>(&self, writer: &mut W) -> Result<()> {
        write_u64(writer, self.size)?;
        write_u64(writer, self.words.len() as u64)?;
        for &word in &self.words {
            write_u64(writer, word)?;
        }
        write_u64(writer, self.ranks.len() as u64)?;
        for &rank in &self.ranks {
            write_u64(writer, rank)?;
        }
        Ok(())
    }
}

pub(super) trait MphfIndex: Copy + Send + Sync {
    fn from_index(index: usize) -> Self;
    fn index(self) -> usize;
}
impl MphfIndex for u32 {
    fn from_index(index: usize) -> Self {
        index as u32
    }
    fn index(self) -> usize {
        self as usize
    }
}
impl MphfIndex for usize {
    fn from_index(index: usize) -> Self {
        index
    }
    fn index(self) -> usize {
        self
    }
}

pub(super) trait MphfBucket: Copy + Send + Sync {
    fn from_bucket(bucket: u64) -> Self;
    fn bucket(self) -> u64;
}
impl MphfBucket for u32 {
    fn from_bucket(bucket: u64) -> Self {
        bucket as u32
    }
    fn bucket(self) -> u64 {
        self as u64
    }
}
impl MphfBucket for u64 {
    fn from_bucket(bucket: u64) -> Self {
        bucket
    }
    fn bucket(self) -> u64 {
        self
    }
}

pub(super) fn write_boo_mphf(out: &mut Vec<u8>, keys: &[u64], threads: usize) -> Result<()> {
    if boo_level_domains(keys.len())?
        .iter()
        .all(|&domain| domain <= u32::MAX as u64)
    {
        write_boo_mphf_ids::<u32, u32>(out, keys, threads)
    } else if keys.len() <= u32::MAX as usize {
        write_boo_mphf_ids::<u32, u64>(out, keys, threads)
    } else {
        write_boo_mphf_ids::<usize, u64>(out, keys, threads)
    }
}

#[derive(Clone, Copy, Default)]
struct MphfMarks {
    seen: u64,
    collisions: u64,
}

pub(super) fn write_boo_mphf_ids<Id: MphfIndex, Bucket: MphfBucket>(
    out: &mut Vec<u8>,
    keys: &[u64],
    threads: usize,
) -> Result<()> {
    write_boo_mphf_ids_first::<Id, Bucket>(out, keys, threads, keys.len() > (1 << 28))
}

pub(super) fn write_boo_mphf_ids_first<Id: MphfIndex, Bucket: MphfBucket>(
    out: &mut Vec<u8>,
    keys: &[u64],
    threads: usize,
    implicit_first: bool,
) -> Result<()> {
    let pool = if threads > 1 && keys.len() >= 1024 {
        Some(IndexWorkers::new(threads)?)
    } else {
        None
    };
    let domains = boo_level_domains(keys.len())?;
    // The first level visits every key in numeric order. Generate those IDs
    // instead of allocating and scanning a graph-sized identity vector.
    // Subsequent levels retain their exact stable survivor order.
    let mut survivors: Vec<Id> = if implicit_first {
        Vec::new()
    } else {
        (0..keys.len()).map(Id::from_index).collect()
    };
    let mut bitvectors = Vec::with_capacity(BOO_LEVELS);
    let mut offset = 0u64;

    for (level, &domain) in domains.iter().enumerate() {
        let mut bitvector = BooBitVector::new(domain)?;
        let all_keys = implicit_first && level == 0;
        let count = if all_keys {
            keys.len()
        } else {
            survivors.len()
        };
        let id_at = |index: usize| {
            if all_keys {
                Id::from_index(index)
            } else {
                survivors[index]
            }
        };
        if level + 1 < BOO_LEVELS && count != 0 {
            if let Some(pool) = pool.as_ref().filter(|_| count >= 1024) {
                let chunk = count.div_ceil(threads.min(count.div_ceil(256)).max(1));
                let reduced = ReducedDomain::new(domain);
                let shards = pool.install(|| {
                    (0..count.div_ceil(chunk))
                        .into_par_iter()
                        .map(|batch| {
                            let first = batch * chunk;
                            let end = (first + chunk).min(count);
                            // Both fields are touched for each hashed bit. Keep
                            // them on the same cache line without atomic updates.
                            let mut marks = vec![MphfMarks::default(); bitvector.words.len()];
                            let buckets = (first..end)
                                .map(|index| {
                                    let key_idx = id_at(index);
                                    let bucket = reduced.remainder(boo_level_hash(
                                        &keys[key_idx.index()].to_be_bytes(),
                                        level,
                                    ));
                                    let word = (bucket >> 6) as usize;
                                    let mask = 1 << (bucket & 63);
                                    let mark = &mut marks[word];
                                    mark.collisions |= mark.seen & mask;
                                    mark.seen |= mask;
                                    Bucket::from_bucket(bucket)
                                })
                                .collect::<Vec<_>>();
                            (marks, buckets)
                        })
                        .collect::<Vec<_>>()
                });
                pool.install(|| {
                    bitvector
                        .words
                        .par_iter_mut()
                        .enumerate()
                        .for_each(|(word, dest)| {
                            let mut seen = 0;
                            let mut collisions = 0;
                            for (marks, _) in &shards {
                                let mark = marks[word];
                                collisions |= mark.collisions | (seen & mark.seen);
                                seen |= mark.seen;
                            }
                            *dest = seen & !collisions;
                        })
                });
                // The merged bitvector owns every collision decision. Marks
                // are dead now; release them before allocating survivor shards.
                let shards = pool.install(|| {
                    shards
                        .into_par_iter()
                        .map(|(marks, buckets)| {
                            crate::index::scratch_pages::retire(marks);
                            buckets
                        })
                        .collect::<Vec<_>>()
                });
                if !implicit_first || all_keys {
                    let next = pool.install(|| {
                        (0..count.div_ceil(chunk))
                            .into_par_iter()
                            .zip(shards.par_iter())
                            .map(|(batch, buckets)| {
                                buckets
                                    .iter()
                                    .enumerate()
                                    .filter_map(|(local, &bucket)| {
                                        let key_idx = id_at(batch * chunk + local);
                                        let bucket = bucket.bucket();
                                        (bitvector.words[(bucket >> 6) as usize]
                                            & (1 << (bucket & 63))
                                            == 0)
                                            .then_some(key_idx)
                                    })
                                    .collect::<Vec<_>>()
                            })
                            .collect::<Vec<_>>()
                    });
                    // Every next ID has been copied. Neither original IDs nor
                    // per-key buckets need to overlap the final concatenation.
                    pool.install(|| {
                        shards
                            .into_par_iter()
                            .for_each(crate::index::scratch_pages::retire)
                    });
                    crate::index::scratch_pages::retire(std::mem::take(&mut survivors));
                    let mut combined = Vec::with_capacity(next.iter().map(Vec::len).sum());
                    for mut shard in next {
                        let written = shard.len();
                        combined.append(&mut shard);
                        crate::index::scratch_pages::discard_unused(&mut shard, written);
                    }
                    survivors = combined;
                } else {
                    // Each disjoint source chunk compacts forward, never
                    // overwriting an unread ID. Join before concatenating the
                    // live prefixes in original chunk order within this same
                    // allocation; both moves preserve stable survivor order.
                    let counts = pool.install(|| {
                        survivors
                            .par_chunks_mut(chunk)
                            .zip(shards.par_iter())
                            .map(|(ids, buckets)| {
                                let mut written = 0;
                                for (read, &bucket) in buckets.iter().enumerate() {
                                    let bucket = bucket.bucket();
                                    if bitvector.words[(bucket >> 6) as usize]
                                        & (1 << (bucket & 63))
                                        == 0
                                    {
                                        ids[written] = ids[read];
                                        written += 1;
                                    }
                                }
                                written
                            })
                            .collect::<Vec<_>>()
                    });
                    pool.install(|| {
                        shards
                            .into_par_iter()
                            .for_each(crate::index::scratch_pages::retire)
                    });
                    let old_len = survivors.len();
                    let mut written = 0;
                    for (batch, count) in counts.into_iter().enumerate() {
                        let first = batch * chunk;
                        survivors.copy_within(first..first + count, written);
                        written += count;
                    }
                    survivors.truncate(written);
                    crate::index::scratch_pages::discard_unused(&mut survivors, old_len);
                }
            } else {
                let domain_usize = usize::try_from(domain)
                    .map_err(|_| Error::InvalidFormat("MPHF domain too large".into()))?;
                let mut counts = vec![0u8; domain_usize];
                let mut buckets = Vec::with_capacity(count);
                let reduced_domain = ReducedDomain::new(domain);
                for index in 0..count {
                    let key_idx = id_at(index);
                    let bucket = reduced_domain
                        .remainder(boo_level_hash(&keys[key_idx.index()].to_be_bytes(), level));
                    let count = &mut counts[bucket as usize];
                    *count = count.saturating_add(1).min(2);
                    buckets.push(Bucket::from_bucket(bucket));
                }

                let mut next_survivors = Vec::new();
                for (index, bucket) in buckets.into_iter().enumerate() {
                    let key_idx = id_at(index);
                    let bucket = bucket.bucket();
                    if counts[bucket as usize] == 1 {
                        bitvector.set(bucket);
                    } else {
                        next_survivors.push(key_idx);
                    }
                }
                survivors = next_survivors;
            }
        }
        offset = bitvector.build_ranks(offset);
        bitvectors.push(bitvector);
    }

    write_f64(out, BOO_GAMMA)?;
    write_i32(out, BOO_LEVELS as i32)?;
    write_u64(out, offset)?;
    write_u64(out, keys.len() as u64)?;
    for bitvector in &bitvectors {
        bitvector.write(out)?;
    }

    write_u64(out, survivors.len() as u64)?;
    for (idx, &key_idx) in survivors.iter().enumerate() {
        out.write_all(&keys[key_idx.index()].to_be_bytes())?;
        write_u64(out, idx as u64)?;
    }
    Ok(())
}

fn boo_level_domains(n: usize) -> Result<Vec<u64>> {
    let nelem = n as f64;
    let hash_domain = (nelem * BOO_GAMMA).ceil();
    let proba_collision = if n == 0 {
        1.0
    } else {
        1.0 - (((BOO_GAMMA * nelem - 1.0) / (BOO_GAMMA * nelem)).powf(nelem - 1.0))
    };
    let mut out = Vec::with_capacity(BOO_LEVELS);
    for level in 0..BOO_LEVELS {
        let raw = (hash_domain * proba_collision.powi(level as i32)) as u64;
        let domain = raw
            .checked_add(63)
            .map(|v| (v / 64) * 64)
            .ok_or_else(|| Error::InvalidFormat("MPHF domain too large".into()))?;
        out.push(domain.max(64));
    }
    Ok(out)
}

fn boo_level_hash(key: &[u8; 8], level: usize) -> u64 {
    if level == 0 {
        return minimizer_hash(key, 0xAAAA_AAAA_5555_5555);
    }
    let mut s0 = minimizer_hash(key, 0xAAAA_AAAA_5555_5555);
    let mut s1 = minimizer_hash(key, 0x3333_3333_CCCC_CCCC);
    if level == 1 {
        return s1;
    }
    let mut out = 0;
    for _ in 2..=level {
        out = xorshift_next(&mut s0, &mut s1);
    }
    out
}
