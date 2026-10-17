//! Build aligned minimizer descriptors while looking up bitmap-selected positions.
use super::{
    BitmapPositions, BooPhf, EncodedUnitigs, LOAD_BATCH_POSITIONS, LOAD_MIN_PAR_CHUNK, MinzOffsets,
    MinzPositionIndex, UnitigBlockState, advance_unitig_state,
};
use crate::pseudoalign::fingerprints::{self, DISABLED};
use crate::timing::{self, Stage};
use rayon::prelude::*;
use rdst::RadixSort;
use std::sync::atomic::{AtomicU16, Ordering};

const ORD_BITS: u32 = 26;
const ORD_LIMIT: usize = 1 << ORD_BITS;
const ORD_MASK: u64 = (1 << ORD_BITS) - 1;
const UID_LIMIT: usize = (1 << 21) - 1;
const UID_MASK: u64 = UID_LIMIT as u64;
const MULTI_TAG: u64 = UID_MASK << 42;
const LOCAL_BITS: u32 = 20;
const LOCAL_MASK: u64 = (1 << LOCAL_BITS) - 1;
const LOOKUP_KEY_MASK: u64 = (1 << 46) - 1;

pub(super) struct Context<'a> {
    pub(super) encoded: &'a EncodedUnitigs,
    pub(super) mphf: &'a BooPhf,
    pub(super) pool: &'a rayon::ThreadPool,
    pub(super) buckets: usize,
    pub(super) k: usize,
}
pub(super) fn eligible(
    ctx: &Context<'_>,
    ordinary: &[BitmapPositions],
    km: &[BitmapPositions],
    special: usize,
) -> bool {
    usize::BITS == 64
        && ctx.pool.current_num_threads() > 1
        && ctx.k == 31
        && ctx.encoded.g == 23
        && ordinary.len() <= 1
        && km.len() <= 1
        && ctx.buckets < ORD_LIMIT
        && ctx.encoded.unitig_bases.metadata.len() < UID_LIMIT
        && ctx.encoded.km_bases.metadata.len() < UID_LIMIT
        && (0..ctx.encoded.unitig_bases.metadata.len()).all(|i| {
            ctx.encoded
                .unitig_bases
                .len(i)
                .is_some_and(|n| n < (1 << 16))
        })
        && ordinary
            .iter()
            .chain(km)
            .map(BitmapPositions::len)
            .chain(std::iter::once(special))
            .try_fold(0usize, usize::checked_add)
            .is_some_and(|n| n < ORD_LIMIT)
}

// Packed only during parallel lookup, never in the runtime descriptor table.
// 46 key bits (bucket26/local20), position38 and fingerprint12 fit in 96 bits.
#[derive(Clone, Copy)]
#[repr(C, packed(4))]
struct LookupRecord {
    low: u64,
    high: u32,
}
const _: () = assert!(std::mem::size_of::<LookupRecord>() == 12);
impl LookupRecord {
    fn new(bucket: u32, local: u32, position: u64, fingerprint: u16) -> Self {
        debug_assert!(fingerprint <= DISABLED);
        let key = (u64::from(bucket) << LOCAL_BITS) | u64::from(local);
        Self {
            low: key | (position << 46),
            high: ((position >> 18) as u32) | (u32::from(fingerprint) << 20),
        }
    }
    fn key(self) -> u64 {
        self.low & LOOKUP_KEY_MASK
    }
    fn bucket(self) -> usize {
        (self.key() >> LOCAL_BITS) as usize
    }
    fn local(self) -> u32 {
        (self.key() & LOCAL_MASK) as u32
    }
    fn position(self) -> u64 {
        (self.low >> 46) | (u64::from(self.high & ((1 << 20) - 1)) << 18)
    }
    fn fingerprint(self) -> u16 {
        (self.high >> 20) as u16
    }
}
impl rdst::RadixKey for LookupRecord {
    // Fill the leading radix byte: all46 key bits retain their order.
    const LEVELS: usize = 6;
    fn get_level(&self, level: usize) -> u8 {
        ((self.key() << 2) >> (level * 8)) as u8
    }
}
#[derive(Clone, Copy)]
#[repr(C)]
struct Event {
    key: u64,
    position: u64,
}
impl Event {
    fn new(bucket: usize, ordinal: u32, position: u64) -> Self {
        Self {
            key: ((bucket as u64) << ORD_BITS) | u64::from(ordinal),
            position,
        }
    }
    fn bucket(self) -> usize {
        (self.key >> ORD_BITS) as usize
    }
}
impl rdst::RadixKey for Event {
    // Bucket26 + ordinal26 use52 bits; pad to seven complete radix bytes.
    const LEVELS: usize = 7;
    fn get_level(&self, level: usize) -> u8 {
        ((self.key << 4) >> (level * 8)) as u8
    }
}
fn position_code(position: u64) -> Option<u64> {
    let uid = position >> 32;
    let rel = position & 0x7fff_ffff;
    (uid < UID_MASK && rel < (1 << 16))
        .then_some((uid << 16) | rel | (((position >> 31) & 1) << 37))
}
fn decode_position(code: u64) -> u64 {
    ((code >> 16 & UID_MASK) << 32) | (code & 0xffff) | ((code >> 37 & 1) << 31)
}

struct Builder {
    states: Vec<u64>,
    events: Vec<Event>,
}
impl Builder {
    fn new(buckets: usize, input: usize) -> Self {
        Self {
            states: vec![0; buckets],
            events: Vec::with_capacity(input.saturating_sub(buckets).saturating_mul(2).min(input)),
        }
    }
    #[inline]
    fn record_code(&mut self, bucket: usize, code: u64, ordinal: u32) {
        record_code_state(
            &mut self.states[bucket],
            &mut self.events,
            bucket,
            code,
            ordinal,
        );
    }
    fn record_existing(&mut self, bucket: usize, old: u64, current: Event) {
        record_existing_state(
            &mut self.states[bucket],
            &mut self.events,
            bucket,
            old,
            current,
        );
    }
    fn record_raw(&mut self, bucket: usize, position: u64, ordinal: u32) {
        if let Some(code) = position_code(position) {
            self.record_code(bucket, code, ordinal);
            return;
        }
        let old = self.states[bucket];
        let current = Event::new(bucket, ordinal, position);
        if old == 0 {
            self.states[bucket] = MULTI_TAG | 1;
            self.events.push(current);
        } else {
            self.record_existing(bucket, old, current);
        }
    }
    fn finish(mut self, pool: &rayon::ThreadPool) -> MinzPositionIndex {
        if self.events.len() > 1 {
            pool.install(|| self.events.radix_sort_builder().with_parallel(true).sort());
        }
        pool.install(|| {
            self.states.par_iter_mut().for_each(|state| {
                if *state == 0 {
                    *state = 1 << 63;
                } else if *state & MULTI_TAG != MULTI_TAG {
                    *state = decode_position(*state >> ORD_BITS);
                }
            })
        });
        let mut at = 0;
        while at < self.events.len() {
            let bucket = self.events[at].bucket();
            let mut end = at + 1;
            while end < self.events.len() && self.events[end].bucket() == bucket {
                end += 1;
            }
            debug_assert_eq!(self.states[bucket] & ORD_MASK, (end - at) as u64);
            let crowded = self.events[at..end].iter().any(|event| {
                event.position >> 32 == u32::MAX as u64 && event.position & 0x8000_0000 != 0
            });
            self.states[bucket] = (1 << 63)
                | (((end - at) as u64) << 32)
                | at as u64
                | if crowded { 1 << 62 } else { 0 };
            at = end;
        }
        let values = self.events.iter().map(|event| event.position).collect();
        if timing::is_enabled() {
            eprintln!("[index-one-pass] overflow_positions={}", self.events.len());
        }
        crate::index::scratch_pages::retire(self.events);
        MinzPositionIndex {
            offsets: MinzOffsets::Compact(Vec::new()),
            values,
            inline: Some(self.states),
            inline_building: false,
        }
    }
}
#[inline]
fn record_code_state(
    state: &mut u64,
    events: &mut Vec<Event>,
    bucket: usize,
    code: u64,
    ordinal: u32,
) {
    let old = *state;
    if old == 0 {
        *state = (code << ORD_BITS) | (u64::from(ordinal) + 1);
        return;
    }
    record_existing_state(
        state,
        events,
        bucket,
        old,
        Event::new(bucket, ordinal, decode_position(code)),
    );
}
fn record_existing_state(
    state: &mut u64,
    events: &mut Vec<Event>,
    bucket: usize,
    old: u64,
    current: Event,
) {
    if old & MULTI_TAG == MULTI_TAG {
        *state = old + 1;
    } else {
        events.push(Event::new(
            bucket,
            ((old & ORD_MASK) - 1) as u32,
            decode_position(old >> ORD_BITS),
        ));
        *state = MULTI_TAG | 2;
    }
    events.push(current);
}

fn merge_ordinary_chunks(
    pool: &rayon::ThreadPool,
    builder: &mut Builder,
    fingerprints: &[AtomicU16],
    chunks: &[Vec<LookupRecord>],
    base: u32,
) {
    let width = builder
        .states
        .len()
        .div_ceil(pool.current_num_threads())
        .max(1);
    let events: Vec<Vec<Event>> = pool.install(|| {
        builder
            .states
            .par_chunks_mut(width)
            .enumerate()
            .map(|(owner, states)| {
                let first = owner * width;
                let end = first + states.len();
                let mut events = Vec::new();
                // Every source chunk is sorted by bucket, then original ordinal.
                // Slice out this owner's range without rescanning other buckets.
                for chunk in chunks {
                    let lo = chunk.partition_point(|record| record.bucket() < first);
                    let hi = lo + chunk[lo..].partition_point(|record| record.bucket() < end);
                    for &record in &chunk[lo..hi] {
                        let bucket = record.bucket();
                        fingerprints::observe_serial(&fingerprints[bucket], record.fingerprint());
                        record_code_state(
                            &mut states[bucket - first],
                            &mut events,
                            bucket,
                            record.position(),
                            base + record.local(),
                        );
                    }
                }
                events
            })
            .collect()
    });
    for mut shard in events {
        builder.events.append(&mut shard);
    }
}

fn ordinary_batch(
    ctx: &Context<'_>,
    builder: &mut Builder,
    fingerprints: &[AtomicU16],
    positions: &[u64],
    base: u32,
) {
    let chunk_len =
        (positions.len() / (ctx.pool.current_num_threads() * 4).max(1)).max(LOAD_MIN_PAR_CHUNK);
    let chunks: Vec<Vec<LookupRecord>> = ctx.pool.install(|| {
        positions
            .par_chunks(chunk_len)
            .enumerate()
            .map(|(chunk, positions)| {
                let mut records = Vec::with_capacity(positions.len());
                for (i, &position) in positions.iter().enumerate() {
                    if let Some(bytes) = ctx.encoded.unitig_bases.minimizer_at(
                        (position >> 32) as usize,
                        position as u32 as usize,
                        23,
                    ) && let Some(bucket) = ctx.mphf.lookup(&bytes)
                    {
                        let code = ((position >> 32) << 16) | (position & 0xffff);
                        records.push(LookupRecord::new(
                            bucket as u32,
                            (chunk * chunk_len + i) as u32,
                            code,
                            fingerprints::word_fingerprint(u64::from_le_bytes(bytes)),
                        ));
                    }
                }
                records
                    .radix_sort_builder()
                    .with_parallel(false)
                    .with_single_threaded_tuner()
                    .sort();
                records
            })
            .collect()
    });
    if ctx.pool.current_num_threads() > 1
        && positions.len() >= LOAD_MIN_PAR_CHUNK * ctx.pool.current_num_threads()
    {
        merge_ordinary_chunks(ctx.pool, builder, fingerprints, &chunks, base);
    } else {
        for chunk in chunks {
            for record in chunk {
                let bucket = record.bucket();
                fingerprints::observe_serial(&fingerprints[bucket], record.fingerprint());
                builder.record_code(bucket, record.position(), base + record.local());
            }
        }
    }
}
fn km_batch(
    ctx: &Context<'_>,
    builder: &mut Builder,
    fingerprints: &[AtomicU16],
    positions: &[u32],
    base: u32,
) {
    let chunk_len =
        (positions.len() / (ctx.pool.current_num_threads() * 4).max(1)).max(LOAD_MIN_PAR_CHUNK);
    let chunks: Vec<Vec<LookupRecord>> = ctx.pool.install(|| {
        positions
            .par_chunks(chunk_len)
            .enumerate()
            .map(|(chunk, positions)| {
                let mut records = Vec::with_capacity(positions.len());
                for (i, &position) in positions.iter().enumerate() {
                    let uid = position as usize / 9;
                    let rel = position as usize % 9;
                    if let Some(bytes) = ctx.encoded.km_bases.minimizer_at(uid, rel, 23)
                        && let Some(bucket) = ctx.mphf.lookup(&bytes)
                    {
                        // Eligibility bounds UID/position widths. The shared
                        // record retains the KM bit in its spare high padding.
                        let code = ((uid as u64) << 16) | (1 << 37) | rel as u64;
                        records.push(LookupRecord::new(
                            bucket as u32,
                            (chunk * chunk_len + i) as u32,
                            code,
                            DISABLED,
                        ));
                    }
                }
                records
                    .radix_sort_builder()
                    .with_parallel(false)
                    .with_single_threaded_tuner()
                    .sort();
                records
            })
            .collect()
    });
    if ctx.pool.current_num_threads() > 1
        && positions.len() >= LOAD_MIN_PAR_CHUNK * ctx.pool.current_num_threads()
    {
        merge_ordinary_chunks(ctx.pool, builder, fingerprints, &chunks, base);
    } else {
        for chunk in chunks {
            for record in chunk {
                let bucket = record.bucket();
                fingerprints[bucket].store(DISABLED, Ordering::Relaxed);
                builder.record_code(bucket, record.position(), base + record.local());
            }
        }
    }
}
pub(super) fn build(
    ctx: Context<'_>,
    ordinary: Vec<BitmapPositions>,
    km: Vec<BitmapPositions>,
    specials: Vec<(u64, u64)>,
) -> (MinzPositionIndex, Vec<AtomicU16>) {
    debug_assert!(eligible(&ctx, &ordinary, &km, specials.len()));
    let _count = timing::scoped(Stage::MinimizerCountPass);
    let input: usize = ordinary
        .iter()
        .chain(&km)
        .map(BitmapPositions::len)
        .sum::<usize>()
        + specials.len();
    let mut builder = Builder::new(ctx.buckets, input);
    let fingerprints: Vec<_> = (0..ctx.buckets).map(|_| AtomicU16::new(0)).collect();
    let mut ordinal = 0u32;
    let mut batch = Vec::with_capacity(LOAD_BATCH_POSITIONS);
    let mut base = ordinal;
    for bitmap in ordinary {
        if ctx.encoded.unitig_bases.metadata.is_empty() {
            ordinal += bitmap.len() as u32;
            bitmap.retire();
            continue;
        }
        let mut state = UnitigBlockState {
            unitig_id: 0,
            total_unitig_len: 0,
            current_unitig_len: ctx.encoded.unitig_minimizer_count(0) as u64,
        };
        for position in bitmap.iter() {
            advance_unitig_state(&mut state, ctx.encoded, u64::from(position));
            let rel = u64::from(position) - state.total_unitig_len;
            if batch.is_empty() {
                base = ordinal;
            }
            batch.push(((state.unitig_id as u64) << 32) | rel);
            ordinal += 1;
            if batch.len() == LOAD_BATCH_POSITIONS {
                ordinary_batch(&ctx, &mut builder, &fingerprints, &batch, base);
                batch.clear();
            }
        }
        if !batch.is_empty() {
            ordinary_batch(&ctx, &mut builder, &fingerprints, &batch, base);
            batch.clear();
        }
        bitmap.retire();
    }
    crate::index::scratch_pages::retire(batch);
    let mut batch = Vec::with_capacity(LOAD_BATCH_POSITIONS);
    for bitmap in km {
        for position in bitmap.iter() {
            if batch.is_empty() {
                base = ordinal;
            }
            batch.push(position);
            ordinal += 1;
            if batch.len() == LOAD_BATCH_POSITIONS {
                km_batch(&ctx, &mut builder, &fingerprints, &batch, base);
                batch.clear();
            }
        }
        if !batch.is_empty() {
            km_batch(&ctx, &mut builder, &fingerprints, &batch, base);
            batch.clear();
        }
        bitmap.retire();
    }
    crate::index::scratch_pages::retire(batch);
    for &(minimizer, position) in &specials {
        if let Some(bucket) = ctx.mphf.lookup(&minimizer.to_le_bytes()) {
            fingerprints[bucket as usize].store(DISABLED, Ordering::Relaxed);
            builder.record_raw(bucket as usize, position, ordinal);
        }
        ordinal += 1;
    }
    crate::index::scratch_pages::retire(specials);
    debug_assert_eq!(ordinal as usize, input);
    drop(_count);
    let _fill = timing::scoped(Stage::MinimizerFillPass);
    (builder.finish(ctx.pool), fingerprints)
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn records_and_finalization_preserve_zero_positions_full_specials_and_encounter_order() {
        let mut radix_records = Vec::new();
        let mut radix_events = Vec::new();
        for bucket in [0, 1, 255, 256, 1 << 25, ORD_LIMIT - 1] {
            for ordinal in [0, 1, LOCAL_MASK as u32] {
                radix_records.push(LookupRecord::new(bucket as u32, ordinal, 0, 1));
                radix_events.push(Event::new(bucket, ordinal, 0));
            }
        }
        radix_records.reverse();
        radix_events.reverse();
        let mut expected: Vec<_> = radix_records.iter().map(|record| record.key()).collect();
        expected.sort_unstable();
        radix_records
            .radix_sort_builder()
            .with_parallel(false)
            .with_single_threaded_tuner()
            .sort();
        assert_eq!(
            radix_records
                .iter()
                .map(|record| record.key())
                .collect::<Vec<_>>(),
            expected
        );
        let mut expected: Vec<_> = radix_events.iter().map(|event| event.key).collect();
        expected.sort_unstable();
        radix_events
            .radix_sort_builder()
            .with_parallel(false)
            .with_single_threaded_tuner()
            .sort();
        assert_eq!(
            radix_events
                .iter()
                .map(|event| event.key)
                .collect::<Vec<_>>(),
            expected
        );
        for fingerprint in [0, 1, 254, 255, 256, 4094, DISABLED] {
            let record = LookupRecord::new(
                (ORD_LIMIT - 1) as u32,
                LOCAL_MASK as u32,
                ((UID_LIMIT - 1) as u64) << 16 | 65535 | (1 << 37),
                fingerprint,
            );
            assert_eq!(record.bucket(), ORD_LIMIT - 1);
            assert_eq!(record.local(), LOCAL_MASK as u32);
            assert_eq!(
                record.position(),
                ((UID_LIMIT - 1) as u64) << 16 | 65535 | (1 << 37)
            );
            assert_eq!(record.fingerprint(), fingerprint);
        }
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(2)
            .build()
            .unwrap();
        let cases = [
            (1, 0),
            (3, (2 << 32) | 3),
            (1, (4 << 32) | 9),
            (1, 0x8000_0008),
            (3, 0xffff_ffff_0000_0001),
            (1, u64::MAX),
            (2, (7 << 32) | 0x8000_0001),
            (4, (1 << 32) | 1),
            (4, (2 << 32) | 2),
        ];
        for len in [0, 1, cases.len()] {
            let mut builder = Builder::new(5, len);
            let mut expected = vec![Vec::new(); 5];
            for (ordinal, &(bucket, position)) in cases[..len].iter().enumerate().rev() {
                builder.record_raw(bucket, position, ordinal as u32);
            }
            for &(bucket, position) in &cases[..len] {
                expected[bucket].push(position);
            }
            let index = builder.finish(&pool);
            for (bucket, expected) in expected.iter().enumerate() {
                assert_eq!(index.get(bucket), expected);
                assert_eq!(
                    index.get_with_overcrowded(bucket).1,
                    expected
                        .iter()
                        .any(|&p| p >> 32 == u32::MAX as u64 && p & 0x8000_0000 != 0)
                );
            }
        }
        for workers in [1, 2, 8] {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(workers)
                .build()
                .unwrap();
            let mut serial = Builder::new(17, 500);
            let mut parallel = Builder::new(17, 500);
            let serial_fp: Vec<_> = (0..17).map(|_| AtomicU16::new(0)).collect();
            let parallel_fp: Vec<_> = (0..17).map(|_| AtomicU16::new(0)).collect();
            serial.record_raw(4, u64::MAX, 0);
            parallel.record_raw(4, u64::MAX, 0);
            for fp in [&serial_fp, &parallel_fp] {
                fp[4].store(DISABLED, Ordering::Relaxed);
            }
            let chunks: Vec<Vec<LookupRecord>> = (0..3)
                .map(|chunk| {
                    let mut records: Vec<_> = (0..128)
                        .map(|i| {
                            LookupRecord::new(
                                (i * 7 % 17) as u32,
                                (chunk * 128 + i) as u32,
                                ((i as u64) << 16) | (chunk as u64) | (u64::from(i % 3 == 0) << 37),
                                ((i * 31 + chunk * 1365) % 4094 + 1) as u16,
                            )
                        })
                        .collect();
                    records.sort_unstable_by_key(|record| record.key());
                    records
                })
                .collect();
            for base in [37, 1000] {
                for chunk in &chunks {
                    for &record in chunk {
                        fingerprints::observe_serial(
                            &serial_fp[record.bucket()],
                            record.fingerprint(),
                        );
                        serial.record_code(
                            record.bucket(),
                            record.position(),
                            base + record.local(),
                        );
                    }
                }
                merge_ordinary_chunks(&pool, &mut parallel, &parallel_fp, &chunks, base);
            }
            let serial = serial.finish(&pool);
            let parallel = parallel.finish(&pool);
            for bucket in 0..17 {
                assert_eq!(
                    parallel.get_with_overcrowded(bucket),
                    serial.get_with_overcrowded(bucket)
                );
                assert_eq!(
                    parallel_fp[bucket].load(Ordering::Relaxed),
                    serial_fp[bucket].load(Ordering::Relaxed)
                );
            }
        }
    }
}
