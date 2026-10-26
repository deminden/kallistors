//! Bounded source-ordered construction for large ordinary-only minimizer maps.
//! Unsupported or unexpectedly repetitive inputs retain the established loader.
use super::loader_probe::CachedUnitigMinimizers;
use super::one_pass::Context;
use super::*;
use crate::index::bifrost::MphfBatchKeys;
use crate::pseudoalign::fingerprints;

const MULTI: u64 = 1 << 63;
const UID_LIMIT: usize = 1 << 30;
const INPUT_LIMIT: usize = 1 << 29;
const EVENT_BYTES: usize = 800 * 1024 * 1024;

// Original encounter order is (unitig, relative position). No global ordinal
// is needed: every ordinary bitmap is sorted, and only one is supported here.
// Eligibility rejects every unitig of length65536 or greater. With
// g=23 and a valid global bitmap position, advancing to its containing
// unitig proves a relative minimizer offset below65536. Keep that proof
// visible to inlined checked extraction; no unchecked sequence access.
fn unitig_length_eligible(length: usize) -> bool {
    length < (1 << 16)
}
#[inline]
fn relative_offset(raw: u64, unitig_start: u64) -> usize {
    let relative = raw - unitig_start;
    debug_assert!(relative < (1 << 16));
    usize::from(relative as u16)
}

fn source_code(unitig: usize, relative: usize) -> u64 {
    debug_assert!(unitig < UID_LIMIT && relative < (1 << 16));
    ((unitig as u64) << 16) | relative as u64
}
fn position(code: u64) -> u64 {
    ((code >> 16) << 32) | (code & 0xffff)
}

#[derive(Clone, Copy)]
#[repr(C, packed(2))]
struct Event {
    source_low: u32,
    source_high: u16,
    bucket: u32,
}
impl Event {
    fn new(bucket: u32, source: u64) -> Self {
        Self {
            source_low: source as u32,
            source_high: (source >> 32) as u16,
            bucket,
        }
    }
    fn source(self) -> u64 {
        u64::from(self.source_low) | (u64::from(self.source_high) << 32)
    }
}
impl rdst::RadixKey for Event {
    const LEVELS: usize = 10;
    fn get_level(&self, level: usize) -> u8 {
        if level < 6 {
            (self.source() >> (level * 8)) as u8
        } else {
            (self.bucket >> ((level - 6) * 8)) as u8
        }
    }
}

#[derive(Clone, Copy)]
#[repr(C)]
struct Record {
    bucket: u32,
    source_low: u32,
    source_high_and_fingerprint: u32,
}
impl Record {
    fn new(bucket: u32, source: u64, fingerprint: u16) -> Self {
        debug_assert!(source < (1 << 48) && fingerprint <= fingerprints::DISABLED);
        Self {
            bucket,
            source_low: source as u32,
            source_high_and_fingerprint: ((source >> 32) as u32) | (u32::from(fingerprint) << 16),
        }
    }
    fn source(self) -> u64 {
        u64::from(self.source_low) | (u64::from(self.source_high_and_fingerprint & 0xffff) << 32)
    }
    fn fingerprint(self) -> u16 {
        (self.source_high_and_fingerprint >> 16) as u16
    }
}
impl rdst::RadixKey for Record {
    const LEVELS: usize = 4;
    fn get_level(&self, level: usize) -> u8 {
        (self.bucket >> (level * 8)) as u8
    }
}
const _: () = assert!(std::mem::size_of::<Event>() == 10);
const _: () = assert!(std::mem::size_of::<Record>() == 12);

fn event_limit(input: usize) -> usize {
    (EVENT_BYTES / std::mem::size_of::<Event>()).min(input / 3)
}
fn size_eligible(buckets: usize, input: usize) -> bool {
    usize::BITS == 64
        && (1 << 28..u32::MAX as usize).contains(&buckets)
        && (buckets..INPUT_LIMIT).contains(&input)
        && input
            .checked_sub(buckets)
            .and_then(|excess| excess.checked_mul(2))
            .is_some_and(|bound| bound <= event_limit(input))
}
pub(super) fn eligible(
    ctx: &Context<'_>,
    ordinary: &[BitmapPositions],
    km: &[BitmapPositions],
    special: usize,
) -> bool {
    if ctx.k != 31
        || ctx.encoded.g != 23
        || ordinary.len() != 1
        || km.iter().any(|bitmap| bitmap.len() != 0)
        || special != 0
        || !ctx.encoded.km_bases.metadata.is_empty()
        || ctx.encoded.unitig_bases.metadata.is_empty()
        || ctx.encoded.unitig_bases.metadata.len() >= UID_LIMIT
        || ctx.mphf.size() != ctx.buckets as u64
        || !size_eligible(ctx.buckets, ordinary[0].len())
    {
        return false;
    }
    let mut windows = 0u64;
    for unitig in 0..ctx.encoded.unitig_bases.metadata.len() {
        let Some(length) = ctx.encoded.unitig_bases.len(unitig) else {
            return false;
        };
        if !unitig_length_eligible(length) {
            return false;
        }
        windows += ctx.encoded.unitig_minimizer_count(unitig) as u64;
    }
    ordinary[0]
        .last()
        .is_some_and(|last| u64::from(last) < windows)
}

struct Builder {
    states: Vec<u64>,
    events: Vec<Event>,
    limit: usize,
}
impl Builder {
    fn new(buckets: usize, input: usize, limit: usize) -> Self {
        Self {
            states: vec![0; buckets],
            events: Vec::with_capacity(limit.min(input)),
            limit,
        }
    }
    fn retire(self) {
        crate::index::scratch_pages::retire(self.states);
        crate::index::scratch_pages::retire(self.events);
    }
    fn finish(mut self, pool: &ThreadPool) -> (MinzPositionIndex, Vec<AtomicU16>) {
        if self.events.len() > 1 {
            pool.install(|| self.events.radix_sort_builder().with_parallel(true).sort());
        }
        // One ordered sequential-address pass materializes the runtime
        // filter while completing the existing descriptor conversion. No
        // second random-access table is needed during bucket construction.
        let fingerprints = pool.install(|| {
            self.states
                .par_iter_mut()
                .map(|state| {
                    let fingerprint = ((*state & FINGERPRINT_MASK) >> FINGERPRINT_SHIFT) as u16;
                    *state &= !FINGERPRINT_MASK;
                    if *state == 0 {
                        *state = MULTI;
                    } else if *state & MULTI == 0 {
                        *state = position(*state - 1);
                    }
                    AtomicU16::new(fingerprint)
                })
                .collect::<Vec<_>>()
        });
        let mut at = 0;
        while at < self.events.len() {
            let bucket = self.events[at].bucket as usize;
            let mut end = at + 1;
            while end < self.events.len() && self.events[end].bucket as usize == bucket {
                end += 1;
            }
            debug_assert_eq!(self.states[bucket] & !MULTI, (end - at) as u64);
            // Total input<2^29 proves both offset and count fit the existing
            // ordinary-only descriptor without touching marker flags61/62.
            self.states[bucket] = MULTI | (((end - at) as u64) << 32) | at as u64;
            at = end;
        }
        if timing::is_enabled() {
            eprintln!(
                "[index-large-one-pass] overflow_positions={} event_limit={}",
                self.events.len(),
                self.limit
            );
        }
        let values = if pool.current_num_threads() > 1 {
            pool.install(|| {
                self.events
                    .par_iter()
                    .map(|event| position(event.source()))
                    .collect::<Vec<_>>()
            })
        } else {
            self.events
                .iter()
                .map(|event| position(event.source()))
                .collect::<Vec<_>>()
        };
        crate::index::scratch_pages::retire(self.events);
        (
            MinzPositionIndex {
                offsets: MinzOffsets::Compact(Vec::new()),
                values,
                inline: Some(self.states),
                inline_building: false,
            },
            fingerprints,
        )
    }
}
// Source+1 uses at most47 bits (even the maximum46-bit source+1 fits).
// Multi-row counts are bounded below2^29, and bit63 retains its existing flag.
// Bits47..58 can carry fingerprint consensus until finish clears them BEFORE
// converting singleton sources to runtime positions (which can use bit61).
const FINGERPRINT_SHIFT: u32 = 47;
const FINGERPRINT_MASK: u64 = (fingerprints::DISABLED as u64) << FINGERPRINT_SHIFT;

#[inline]
fn record_fingerprinted(
    state: &mut u64,
    events: &mut Vec<Event>,
    bucket: u32,
    source: u64,
    fingerprint: u16,
) {
    debug_assert!(fingerprint <= fingerprints::DISABLED);
    if *state == 0 {
        // Both descriptor and consensus are unseen. This case emits no
        // overflow events; source+1 is nonzero and cannot overlap FP bits.
        let observed = if fingerprint == 0 {
            fingerprints::DISABLED
        } else {
            fingerprint
        };
        debug_assert_eq!((source + 1) & FINGERPRINT_MASK, 0);
        *state = (source + 1) | (u64::from(observed) << FINGERPRINT_SHIFT);
        return;
    }
    let old = ((*state & FINGERPRINT_MASK) >> FINGERPRINT_SHIFT) as u16;
    let observed = if fingerprint == 0
        || fingerprint == fingerprints::DISABLED
        || (old != 0 && old != fingerprint)
    {
        fingerprints::DISABLED
    } else {
        fingerprint
    };
    let mut descriptor = *state & !FINGERPRINT_MASK;
    record(&mut descriptor, events, bucket, source);
    debug_assert_eq!(descriptor & FINGERPRINT_MASK, 0);
    *state = descriptor | (u64::from(observed) << FINGERPRINT_SHIFT);
}

fn record(state: &mut u64, events: &mut Vec<Event>, bucket: u32, source: u64) {
    let old = *state;
    if old == 0 {
        *state = source + 1;
    } else {
        if old & MULTI == 0 {
            events.push(Event::new(bucket, old - 1));
            *state = MULTI | 2;
        } else {
            *state = old + 1;
        }
        events.push(Event::new(bucket, source));
    }
}

// Called only with a live in-bounds destination; the vector never moves
// during a flush. Hardware hints do not read/write any Rust value, retain a
// pointer, change publication order, or require a separate allocation.
#[inline]
fn prefetch_owned_descriptor(slot: &u64) {
    #[cfg(target_arch = "aarch64")]
    // SAFETY: an immutable in-bounds address supplies a cache hint only.
    unsafe {
        core::arch::asm!("prfm pstl1keep, [{address}]", address=in(reg) slot as *const u64,
            options(nostack, readonly, preserves_flags));
    }
    #[cfg(target_arch = "x86_64")]
    // SAFETY: SSE is in the x86-64 baseline and the supplied address is live.
    unsafe {
        core::arch::x86_64::_mm_prefetch(
            (slot as *const u64).cast::<i8>(),
            core::arch::x86_64::_MM_HINT_T0,
        );
    }
    #[cfg(not(any(target_arch = "aarch64", target_arch = "x86_64")))]
    let _ = slot;
}

// A single owner already has exclusive access to every bucket. Publish
// original encounter order directly instead of sorting owner-routing records.
// The common caller checks the event bound before any mutation. Final multi
// rows use the same full source-order sort and descriptor construction.
fn flush_serial<const BATCH: usize>(
    ctx: &Context<'_>,
    builder: &mut Builder,
    positions: &[u32],
    boundary: &mut UnitigBlockState,
) {
    let mut probe = CachedUnitigMinimizers::new(ctx.encoded);
    for queries in positions.chunks(BATCH) {
        let mut keys = MphfBatchKeys::<BATCH>::new();
        let mut sources = [0; BATCH];
        for (lane, &raw) in queries.iter().enumerate() {
            let raw = u64::from(raw);
            advance_unitig_state(boundary, ctx.encoded, raw);
            let relative = relative_offset(raw, boundary.total_unitig_len);
            sources[lane] = source_code(boundary.unitig_id, relative);
            keys.set(lane, probe.at23(boundary.unitig_id, relative));
        }
        let buckets = ctx.mphf.lookup_batch(&keys);
        for bucket in buckets.found_ranks() {
            if let Some(destination) = builder.states.get(bucket as usize) {
                prefetch_owned_descriptor(destination);
            }
        }
        for ((bucket, key), &source) in buckets
            .into_iter()
            .zip(keys.iter())
            .zip(&sources)
            .take(queries.len())
        {
            if let (Some(bucket), Some(key)) = (bucket, key) {
                let bucket = bucket as usize;
                record_fingerprinted(
                    &mut builder.states[bucket],
                    &mut builder.events,
                    bucket as u32,
                    source,
                    fingerprints::word_fingerprint(u64::from_le_bytes(key)),
                );
            }
        }
    }
    debug_assert!(builder.events.len() <= builder.limit);
}

// Only the owner-range boundary predicates need to be monotone.
// Stable interval routing preserves encounter order inside every bucket;
// full bucket sorting is not required for the owned descriptor updates.
// Out-of-domain records keep a final sentinel group, so even malformed
// lookup output remains outside each owner's [first,end) slice as before.
fn route_records(records: Vec<Record>, width: usize, buckets: usize) -> Vec<Record> {
    if records.len() < 2 {
        return records;
    }
    let owners = buckets.div_ceil(width);
    let owner = |bucket: u32| {
        let bucket = bucket as usize;
        if bucket < buckets {
            bucket / width
        } else {
            owners
        }
    };
    let mut next = vec![0usize; owners + 1];
    for record in &records {
        next[owner(record.bucket)] += 1;
    }
    let mut total = 0;
    for count in &mut next {
        let length = *count;
        *count = total;
        total += length;
    }
    debug_assert_eq!(total, records.len());
    let mut routed = vec![Record::new(0, 0, 0); records.len()];
    for record in &records {
        let destination = &mut next[owner(record.bucket)];
        routed[*destination] = *record;
        *destination += 1;
    }
    crate::index::scratch_pages::retire(records);
    routed
}

fn flush<const BATCH: usize, const MIN_CHUNK: usize>(
    ctx: &Context<'_>,
    builder: &mut Builder,
    positions: &[u32],
    boundary: &mut UnitigBlockState,
    seen: usize,
) -> bool {
    // At most two events per new position, and never more events than all
    // positions seen. Refuse before mutation/allocation if the bound is too big.
    let maximum = builder
        .events
        .len()
        .saturating_add(positions.len().saturating_mul(2))
        .min(seen.saturating_add(positions.len()));
    if maximum > builder.limit {
        return false;
    }
    if ctx.pool.current_num_threads() == 1 {
        flush_serial::<BATCH>(ctx, builder, positions, boundary);
        return true;
    }
    let chunk = (positions.len() / ctx.pool.current_num_threads().saturating_mul(4).max(1))
        .max(MIN_CHUNK)
        .next_power_of_two();
    let mut initial = Vec::with_capacity(positions.len().div_ceil(chunk));
    for entries in positions.chunks(chunk) {
        advance_unitig_state(boundary, ctx.encoded, u64::from(entries[0]));
        initial.push(*boundary);
    }
    let width = ctx.buckets.div_ceil(ctx.pool.current_num_threads()).max(1);
    let chunks: Vec<Vec<Record>> = ctx.pool.install(|| {
        positions
            .par_chunks(chunk)
            .zip(initial.par_iter())
            .map(|(entries, &start)| {
                let mut state = start;
                let mut probe = CachedUnitigMinimizers::new(ctx.encoded);
                let mut records = Vec::with_capacity(entries.len());
                for queries in entries.chunks(BATCH) {
                    let mut keys = MphfBatchKeys::<BATCH>::new();
                    let mut sources = [0; BATCH];
                    for (lane, &raw) in queries.iter().enumerate() {
                        let raw = u64::from(raw);
                        advance_unitig_state(&mut state, ctx.encoded, raw);
                        let relative = relative_offset(raw, state.total_unitig_len);
                        sources[lane] = source_code(state.unitig_id, relative);
                        keys.set(lane, probe.at23(state.unitig_id, relative));
                    }
                    let buckets = ctx.mphf.lookup_batch(&keys);
                    for ((bucket, key), &source) in buckets
                        .into_iter()
                        .zip(keys.iter())
                        .zip(&sources)
                        .take(queries.len())
                    {
                        if let (Some(bucket), Some(key)) = (bucket, key) {
                            records.push(Record::new(
                                bucket as u32,
                                source,
                                fingerprints::word_fingerprint(u64::from_le_bytes(key)),
                            ));
                        }
                    }
                }
                route_records(records, width, ctx.buckets)
            })
            .collect()
    });
    let shards: Vec<Vec<Event>> = ctx.pool.install(|| {
        builder
            .states
            .par_chunks_mut(width)
            .enumerate()
            .map(|(owner, states)| {
                let first = owner * width;
                let end = first + states.len();
                let mut events = Vec::new();
                for chunk in &chunks {
                    let lo = chunk.partition_point(|record| (record.bucket as usize) < first);
                    let hi =
                        lo + chunk[lo..].partition_point(|record| (record.bucket as usize) < end);
                    let records = &chunk[lo..hi];
                    for (index, &item) in records.iter().enumerate() {
                        if let Some(future) = records.get(index + 32) {
                            // This owner-grouped slice belongs to the same exclusive
                            // bucket range. Borrow the future destination only
                            // for a hint, before updating the current record.
                            if let Some(destination) = states.get(future.bucket as usize - first) {
                                prefetch_owned_descriptor(destination);
                            }
                        }
                        let bucket = item.bucket as usize;
                        record_fingerprinted(
                            &mut states[bucket - first],
                            &mut events,
                            item.bucket,
                            item.source(),
                            item.fingerprint(),
                        );
                    }
                }
                events
            })
            .collect()
    });
    for mut shard in shards {
        builder.events.append(&mut shard);
    }
    debug_assert!(builder.events.len() <= builder.limit);
    true
}

pub(super) fn try_build(
    ctx: Context<'_>,
    ordinary: &BitmapPositions,
) -> Option<(MinzPositionIndex, Vec<AtomicU16>)> {
    try_build_sized::<64, LOAD_BATCH_POSITIONS, LOAD_MIN_PAR_CHUNK>(
        &ctx,
        ordinary,
        event_limit(ordinary.len()),
    )
}
fn try_build_sized<const BATCH: usize, const POSITIONS: usize, const MIN_CHUNK: usize>(
    ctx: &Context<'_>,
    ordinary: &BitmapPositions,
    limit: usize,
) -> Option<(MinzPositionIndex, Vec<AtomicU16>)> {
    assert!(POSITIONS > 0 && MIN_CHUNK > 0);
    let count_time = timing::scoped(Stage::MinimizerCountPass);
    let mut builder = Builder::new(ctx.buckets, ordinary.len(), limit);
    let mut boundary = UnitigBlockState {
        unitig_id: 0,
        total_unitig_len: 0,
        current_unitig_len: ctx.encoded.unitig_minimizer_count(0) as u64,
    };
    let mut batch = Vec::with_capacity(POSITIONS);
    let mut seen = 0;
    for raw in ordinary.iter() {
        batch.push(raw);
        if batch.len() == POSITIONS {
            if !flush::<BATCH, MIN_CHUNK>(ctx, &mut builder, &batch, &mut boundary, seen) {
                builder.retire();
                return None;
            }
            seen += batch.len();
            batch.clear();
        }
    }
    if !batch.is_empty()
        && !flush::<BATCH, MIN_CHUNK>(ctx, &mut builder, &batch, &mut boundary, seen)
    {
        builder.retire();
        return None;
    }
    crate::index::scratch_pages::retire(batch);
    drop(count_time);
    let _finish_time = timing::scoped(Stage::MinimizerFillPass);
    Some(builder.finish(ctx.pool))
}

#[cfg(test)]
mod tests;
