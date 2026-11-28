use super::utils::EcCountMap;
use std::ops::Range;
use std::sync::Arc;
use std::thread;
use std::time::Instant;

use crossbeam_channel::{Sender, bounded};

use crate::io::{PackedSeqBatch, PackedSeqBatchSource, ReadSource};
use crate::timing::{self, Stage};
use crate::{Error, Result};

use super::{BifrostIndex, EcCounts, FragmentFilter, PseudoalignOptions, Strand};

// Original kallisto's paired FASTQ buffer holds both sequences and one trailing
// byte per mate. Its sampling quota belongs to that buffer, not a parser chunk.
const PAIRED_SAMPLING_BUFFER_BYTES: usize = 1 << 23;

struct PairedChunk {
    left: Arc<PackedSeqBatch>,
    right: Arc<PackedSeqBatch>,
    reads: Range<usize>,
}

#[derive(Default)]
struct PairedBatchBuilder {
    chunks: Vec<PairedChunk>,
    bytes: usize,
}

impl PairedBatchBuilder {
    fn finish(&mut self) -> Option<Vec<PairedChunk>> {
        self.bytes = 0;
        (!self.chunks.is_empty()).then(|| std::mem::take(&mut self.chunks))
    }

    fn push(
        &mut self,
        left: PackedSeqBatch,
        right: PackedSeqBatch,
        byte_limit: usize,
    ) -> Result<Vec<Vec<PairedChunk>>> {
        if left.len() != right.len() {
            return Err(Error::InvalidFormat("paired FASTQ length mismatch".into()));
        }
        let left = Arc::new(left);
        let right = Arc::new(right);
        let mut ready = Vec::new();
        if byte_limit == 0 {
            if let Some(pending) = self.finish() {
                ready.push(pending);
            }
            if !left.is_empty() {
                ready.push(vec![PairedChunk {
                    reads: 0..left.len(),
                    left,
                    right,
                }]);
            }
            return Ok(ready);
        }
        let mut start = 0;
        for i in 0..left.len() {
            let pair_bytes = left.seq(i).len() + right.seq(i).len() + 2;
            // Upstream uses strict '< limit', retaining a boundary pair for
            // the next buffer. Keep parser batches shared when it splits one.
            if self.bytes + pair_bytes >= byte_limit {
                if start < i {
                    self.chunks.push(PairedChunk {
                        left: left.clone(),
                        right: right.clone(),
                        reads: start..i,
                    });
                }
                if let Some(pending) = self.finish() {
                    ready.push(pending);
                }
                start = i;
            }
            self.bytes += pair_bytes;
            // An oversized pair must still make progress rather than leaving
            // an empty buffer loop. It becomes a single-pair batch.
            if pair_bytes >= byte_limit {
                self.chunks.push(PairedChunk {
                    left: left.clone(),
                    right: right.clone(),
                    reads: i..i + 1,
                });
                ready.push(self.finish().expect("oversized pair is nonempty"));
                start = i + 1;
            }
        }
        if start < left.len() {
            self.chunks.push(PairedChunk {
                reads: start..left.len(),
                left,
                right,
            });
        }
        Ok(ready)
    }
}

fn produce_batches<R: PackedSeqBatchSource>(
    reader: &mut R,
    tx: &Sender<PackedSeqBatch>,
    batch_size: usize,
) -> Result<()> {
    while let Some(batch) = reader.next_packed_seq_batch(batch_size)? {
        if tx.send(batch).is_err() {
            return Ok(());
        }
    }
    Ok(())
}

/// One exact-match table for all workers, with the same total budget as the
/// private tables it replaces, so a k-mer resolved by one worker hits for all.
fn shared_match_table(
    entries_per_worker: usize,
    threads: usize,
) -> Result<Option<Arc<super::core::SharedMatchTable>>> {
    if std::env::var_os("KALLISTORS_PRIVATE_MATCH_CACHE").is_some() {
        return Ok(None);
    }
    super::core::SharedMatchTable::new(entries_per_worker.saturating_mul(threads)).map(Some)
}

/// Fragment-length sampling budget for numbered paired buffers, independent of
/// thread timing.
///
/// Workers take buffers in input order. Buffer `i` (1-based) may record up to
/// `goal - observations(buffers 1..=i - lag)` observations. With `lag` equal
/// to the worker count this is the budget a worker sees when the `lag - 1`
/// buffers handed out just before its own are still running, which is how a
/// shared, submit-on-finish counter behaves when workers run at the same
/// speed. Mapping uses `lag = 1`, which makes the sample independent of the
/// worker count as well as of timing. A worker only ever waits for strictly
/// earlier buffers, so the wait always ends.
struct OrderedSampler {
    goal: usize,
    lag: u64,
    state: std::sync::Mutex<OrderedSamplerState>,
    ready: std::sync::Condvar,
    /// First buffer known to have a zero budget; every later one has too.
    finished_from: std::sync::atomic::AtomicU64,
}

struct OrderedSamplerState {
    /// `prefix[n]`: observations recorded by buffers `1..=n`, for the
    /// contiguous run of completed buffers.
    prefix: Vec<usize>,
    /// Completed buffers beyond that run.
    pending: std::collections::BTreeMap<u64, usize>,
    aborted: bool,
}

impl OrderedSampler {
    fn new(goal: usize, lag: usize) -> Self {
        Self {
            goal,
            lag: lag.max(1) as u64,
            state: std::sync::Mutex::new(OrderedSamplerState {
                prefix: vec![0],
                pending: std::collections::BTreeMap::new(),
                aborted: false,
            }),
            ready: std::sync::Condvar::new(),
            finished_from: std::sync::atomic::AtomicU64::new(if goal == 0 { 0 } else { u64::MAX }),
        }
    }

    /// Whether every buffer handed out from now on has a zero budget.
    fn finished(&self) -> bool {
        self.finished_from
            .load(std::sync::atomic::Ordering::Acquire)
            != u64::MAX
    }

    fn budget(&self, sequence: u64) -> usize {
        // Budgets only shrink with the sequence number, but an earlier buffer
        // can still be finishing after a later one found its budget empty.
        if sequence
            >= self
                .finished_from
                .load(std::sync::atomic::Ordering::Acquire)
        {
            return 0;
        }
        if sequence <= self.lag {
            return self.goal;
        }
        let needed = (sequence - self.lag) as usize;
        let mut state = self.state.lock().expect("fragment sampler");
        while state.prefix.len() <= needed && !state.aborted {
            state = self.ready.wait(state).expect("fragment sampler");
        }
        if state.aborted {
            return 0;
        }
        let budget = self.goal.saturating_sub(state.prefix[needed]);
        if budget == 0 {
            // Prefix sums only grow, so no later buffer can sample either.
            self.finished_from
                .fetch_min(sequence, std::sync::atomic::Ordering::AcqRel);
        }
        budget
    }

    fn complete(&self, sequence: u64, observations: usize) {
        let mut state = self.state.lock().expect("fragment sampler");
        state.pending.insert(sequence, observations);
        loop {
            let next = state.prefix.len() as u64;
            let Some(observed) = state.pending.remove(&next) else {
                break;
            };
            let total = state.prefix[state.prefix.len() - 1] + observed;
            state.prefix.push(total);
        }
        drop(state);
        self.ready.notify_all();
    }

    fn abort(&self) {
        self.state.lock().expect("fragment sampler").aborted = true;
        self.ready.notify_all();
    }
}

struct AbortOnDrop<'a>(&'a OrderedSampler);

impl Drop for AbortOnDrop<'_> {
    fn drop(&mut self) {
        self.0.abort();
    }
}

/// Graph k-mer filter that lets workers reject absent k-mers in one load.
fn shared_kmer_filter(
    index: &BifrostIndex,
    threads: usize,
) -> Option<Arc<super::core::KmerFilter>> {
    if std::env::var_os("KALLISTORS_NO_KMER_FILTER").is_some() {
        return None;
    }
    if let Some(filter) = &index.kmer_filter {
        return Some(Arc::clone(filter));
    }
    let start = Instant::now();
    let filter = super::core::KmerFilter::build(index, threads).map(Arc::new);
    if std::env::var_os("KALLISTORS_TIMINGS").is_some() {
        eprintln!(
            "[kmer-filter] built={} in {:.3}s",
            filter.is_some(),
            start.elapsed().as_secs_f64()
        );
    }
    filter
}

/// Runs `produce` on its own thread, feeding batches through a one-slot queue
/// to `consume` on this thread, so input decoding overlaps mapping.
///
/// A producer error is reported ahead of the consumer's: when the producer
/// fails, its queue closes and the consumer sees a short (or mismatched)
/// stream, which is a consequence rather than the cause.
pub(crate) fn with_producer<T, P, C, O>(produce: P, consume: C) -> Result<O>
where
    T: Send,
    P: FnOnce(&Sender<T>) -> Result<()> + Send,
    C: FnOnce(&mut dyn FnMut() -> Result<Option<T>>) -> Result<O>,
{
    thread::scope(|scope| {
        let (tx, rx) = bounded::<T>(1);
        let read_start = Instant::now();
        let producer = scope.spawn(move || produce(&tx));
        let consumed = consume(&mut || Ok(rx.recv().ok()));
        // Disconnect before joining: an early consumer exit can leave the
        // producer blocked in send even though no more batches are needed.
        drop(rx);
        let produced = producer
            .join()
            .unwrap_or_else(|_| Err(Error::InvalidFormat("producer thread panicked".into())));
        timing::add_duration(Stage::FastqReadDecompress, read_start.elapsed());
        produced?;
        consumed
    })
}

pub fn pseudoalign_single_end_bifrost_with_options_threaded<
    R: ReadSource + PackedSeqBatchSource + Send,
>(
    index: &BifrostIndex,
    reader: &mut R,
    strand: Strand,
    filter: Option<FragmentFilter>,
    options: PseudoalignOptions,
    threads: usize,
) -> Result<EcCounts> {
    if threads <= 1 {
        return single_end_pipeline(
            index,
            &mut || reader.next_packed_seq_batch(4_096),
            strand,
            filter,
            options,
            threads,
        );
    }
    with_producer(
        |tx| produce_batches(reader, tx, 4_096),
        |next| single_end_pipeline(index, next, strand, filter, options, threads),
    )
}

/// Single-end mapping over sequence batches in input order, from any source.
pub(crate) fn single_end_pipeline(
    index: &BifrostIndex,
    next: &mut dyn FnMut() -> Result<Option<PackedSeqBatch>>,
    strand: Strand,
    filter: Option<FragmentFilter>,
    options: PseudoalignOptions,
    threads: usize,
) -> Result<EcCounts> {
    let cache_entries = options
        .match_cache_size
        .allocation(threads)?
        .entries_per_worker();
    if threads <= 1 {
        let _cache_budget = super::core::MatchCacheBudgetScope::enter(cache_entries)?;
        // Synchronous counting uses the same bias-aware sequence batches
        // and cap semantics without spawning a producer or worker.
        let mut pairs = super::core::ScopedEcPairs::new(index, threads.max(1));
        let mut counts = EcCounts::new(options.bias);
        let mut map = EcCountMap::default();
        while let Some(batch) = next()? {
            super::single::pseudoalign_single_end_bifrost_batch_into(
                index,
                batch,
                strand,
                filter,
                options,
                &mut counts,
                &mut map,
                &mut pairs,
            );
        }
        return Ok(counts);
    }

    let threads = threads.max(1);
    let mut merged = EcCounts::new(options.bias);
    let mut merged_map: EcCountMap = EcCountMap::default();
    let mut final_error: Option<Error> = None;
    let key_pool = Arc::new(super::utils::SharedEcPool::default());

    let shared_matches = shared_match_table(cache_entries, threads)?;
    let shared_matches = &shared_matches;
    let kmer_filter = shared_kmer_filter(index, threads);
    let kmer_filter = &kmer_filter;
    thread::scope(|scope| {
        // Rendezvous dispatch bounds read storage independently of worker
        // count without dropping or reordering input.
        let (work_tx, work_rx) = bounded::<PackedSeqBatch>(0);
        let mut handles = Vec::with_capacity(threads);

        let pseudoalign_start = Instant::now();
        for _ in 0..threads {
            let rx = work_rx.clone();
            let key_pool = Arc::clone(&key_pool);
            handles.push(
                scope.spawn(move || -> Result<(EcCounts, super::utils::SharedEcKeys)> {
                    let _cache_budget = match &shared_matches {
                        Some(table) => super::core::MatchCacheBudgetScope::enter_shared(
                            Arc::clone(table),
                            kmer_filter.clone(),
                        )?,
                        None => super::core::MatchCacheBudgetScope::enter(cache_entries)?,
                    };
                    let mut pairs = super::core::ScopedEcPairs::new(index, threads);
                    let mut acc = EcCounts::new(options.bias);
                    let mut acc_map = EcCountMap::shared(key_pool);
                    while let Ok(batch) = rx.recv() {
                        super::single::pseudoalign_single_end_bifrost_batch_into(
                            index,
                            batch,
                            strand,
                            filter,
                            options,
                            &mut acc,
                            &mut acc_map,
                            &mut pairs,
                        );
                    }
                    Ok((acc, acc_map.into_shared_keys()))
                }),
            );
        }

        drop(work_rx);
        loop {
            match next() {
                Ok(Some(batch)) => {
                    if work_tx.send(batch).is_err() {
                        final_error = Some(Error::InvalidFormat("worker queue closed".into()));
                        break;
                    }
                }
                Ok(None) => break,
                Err(err) => {
                    final_error = Some(err);
                    break;
                }
            }
        }
        drop(work_tx);

        key_pool.clear();
        join_workers(handles, &mut merged, &mut merged_map, &mut final_error);
        timing::add_duration(Stage::Pseudoalign, pseudoalign_start.elapsed());
    });

    if let Some(err) = final_error {
        Err(err)
    } else {
        Ok(merged)
    }
}

type WorkerResult = Result<(EcCounts, super::utils::SharedEcKeys)>;

fn join_workers(
    handles: Vec<thread::ScopedJoinHandle<'_, WorkerResult>>,
    merged: &mut EcCounts,
    merged_map: &mut EcCountMap,
    final_error: &mut Option<Error>,
) {
    for handle in handles {
        match handle.join() {
            Ok(Ok((counts, keys))) => {
                let merge_start = Instant::now();
                super::utils::merge_shared_ec_counts(merged, merged_map, counts, keys);
                timing::add_duration(Stage::EcMerge, merge_start.elapsed());
            }
            Ok(Err(err)) => {
                if final_error.is_none() {
                    *final_error = Some(err);
                }
            }
            Err(_) => {
                if final_error.is_none() {
                    *final_error = Some(Error::InvalidFormat("worker thread panicked".into()));
                }
            }
        }
    }
}

fn mate_mismatch() -> Error {
    Error::InvalidFormat("paired FASTQ length mismatch".into())
}

pub fn pseudoalign_paired_bifrost_with_options_threaded<
    R1: ReadSource + PackedSeqBatchSource + Send,
    R2: ReadSource + PackedSeqBatchSource + Send,
>(
    index: &BifrostIndex,
    reader1: &mut R1,
    reader2: &mut R2,
    strand: Strand,
    options: PseudoalignOptions,
    threads: usize,
) -> Result<EcCounts> {
    // Final-pair diagnostics consume names and qualities and own their budget.
    if threads <= 1 && std::env::var_os("KALLISTORS_PAIRED_FINAL_TRACE").is_some() {
        return super::paired::pseudoalign_paired_bifrost_with_options(
            index, reader1, reader2, strand, options,
        );
    }
    if threads <= 1 {
        // Ordinary mapping parses sequence-only batches without workers.
        return paired_pipeline(
            index,
            &mut || {
                let left = reader1.next_packed_seq_batch(super::BATCH_SIZE / 2)?;
                let right = reader2.next_packed_seq_batch(super::BATCH_SIZE / 2)?;
                match (left, right) {
                    (None, None) => Ok(None),
                    (Some(left), Some(right)) => Ok(Some((left, right))),
                    _ => Err(mate_mismatch()),
                }
            },
            strand,
            options,
            threads,
        );
    }
    // One decoding thread per mate file, zipped batch by batch.
    thread::scope(|scope| {
        let (left_tx, left_rx) = bounded::<PackedSeqBatch>(1);
        let (right_tx, right_rx) = bounded::<PackedSeqBatch>(1);
        let read_start = Instant::now();
        let left = scope.spawn(move || produce_batches(reader1, &left_tx, super::BATCH_SIZE / 2));
        let right = scope.spawn(move || produce_batches(reader2, &right_tx, super::BATCH_SIZE / 2));
        let consumed = paired_pipeline(
            index,
            &mut || match (left_rx.recv().ok(), right_rx.recv().ok()) {
                (None, None) => Ok(None),
                (Some(left), Some(right)) => Ok(Some((left, right))),
                _ => Err(mate_mismatch()),
            },
            strand,
            options,
            threads,
        );
        // Release blocked senders before joining producers after an error.
        drop(left_rx);
        drop(right_rx);
        let mut produced = Ok(());
        for producer in [left, right] {
            let result = producer
                .join()
                .unwrap_or_else(|_| Err(Error::InvalidFormat("producer thread panicked".into())));
            if produced.is_ok() {
                produced = result;
            }
        }
        timing::add_duration(Stage::FastqReadDecompress, read_start.elapsed());
        produced?;
        consumed
    })
}

/// Paired mapping over batches of both mates in input order, from any source.
pub(crate) fn paired_pipeline(
    index: &BifrostIndex,
    next: &mut dyn FnMut() -> Result<Option<(PackedSeqBatch, PackedSeqBatch)>>,
    strand: Strand,
    options: PseudoalignOptions,
    threads: usize,
) -> Result<EcCounts> {
    let cache_entries = options
        .match_cache_size
        .allocation(threads)?
        .entries_per_worker();
    if threads <= 1 {
        let _cache_budget = super::core::MatchCacheBudgetScope::enter(cache_entries)?;
        let mut pairs = super::core::ScopedEcPairs::new(index, threads.max(1));
        let mut counts = EcCounts::new(options.bias);
        let mut ec_map = EcCountMap::default();
        let sampler = super::paired::FragmentLengthSampler::new(options.fragment_length_samples);
        while let Some((left, right)) = next()? {
            if left.len() != right.len() {
                return Err(mate_mismatch());
            }
            super::paired::pseudoalign_paired_bifrost_batch_into(
                index,
                &left,
                &right,
                0..left.len(),
                strand,
                options,
                &mut counts,
                &mut ec_map,
                (options.fragment_length_samples != 0).then_some(&sampler),
                &mut pairs,
            )?;
        }
        return Ok(counts);
    }

    let threads = threads.max(1);
    let mut merged = EcCounts::new(options.bias);
    let mut merged_map: EcCountMap = EcCountMap::default();
    let mut final_error: Option<Error> = None;

    let key_pool = Arc::new(super::utils::SharedEcPool::default());
    // Lag 1: each buffer's budget follows every earlier buffer, so the sample
    // is the first qualifying pairs in input order, as on one thread, for any
    // worker count.
    let sampler = OrderedSampler::new(options.fragment_length_samples, 1);
    let shared_matches = shared_match_table(cache_entries, threads)?;
    let shared_matches = &shared_matches;
    let kmer_filter = shared_kmer_filter(index, threads);
    let kmer_filter = &kmer_filter;
    thread::scope(|scope| {
        let (work_tx, work_rx) = bounded::<(u64, Vec<PairedChunk>)>(0);
        let mut sequence = 0u64;
        let mut handles = Vec::with_capacity(threads);

        let pseudoalign_start = Instant::now();
        for _ in 0..threads {
            let rx = work_rx.clone();
            let sampler = &sampler;
            let key_pool = Arc::clone(&key_pool);
            handles.push(scope.spawn(move || -> WorkerResult {
                let _cache_budget = match &shared_matches {
                    Some(table) => super::core::MatchCacheBudgetScope::enter_shared(
                        Arc::clone(table),
                        kmer_filter.clone(),
                    )?,
                    None => super::core::MatchCacheBudgetScope::enter(cache_entries)?,
                };
                let mut pairs = super::core::ScopedEcPairs::new(index, threads);
                let mut acc = EcCounts::new(options.bias);
                let mut acc_map = EcCountMap::shared(key_pool);
                // A worker that stops early must not leave others waiting.
                let _abort = AbortOnDrop(sampler);
                while let Ok((sequence, chunks)) = rx.recv() {
                    // Record this buffer's observations first and wait for
                    // its budget afterwards, when earlier buffers are done.
                    let recording = options.fragment_length_samples != 0 && !sampler.finished();
                    let local_sampler = if recording {
                        super::paired::FragmentLengthSampler::recording()
                    } else {
                        super::paired::FragmentLengthSampler::new(0)
                    };
                    for chunk in chunks {
                        super::paired::pseudoalign_paired_bifrost_batch_into(
                            index,
                            &chunk.left,
                            &chunk.right,
                            chunk.reads,
                            strand,
                            options,
                            &mut acc,
                            &mut acc_map,
                            (options.fragment_length_samples != 0).then_some(&local_sampler),
                            &mut pairs,
                        )?;
                    }
                    let observed = local_sampler.take_recorded();
                    let budget = if recording {
                        sampler.budget(sequence)
                    } else {
                        0
                    };
                    let kept = budget.min(observed.len());
                    super::paired::add_fragment_lengths(&mut acc, &observed[..kept]);
                    sampler.complete(sequence, kept);
                }
                std::mem::forget(_abort);
                Ok((acc, acc_map.into_shared_keys()))
            }));
        }
        drop(work_rx);

        let mut batch_builder = PairedBatchBuilder::default();
        'input: loop {
            let (left, right) = match next() {
                Ok(Some(pair)) => pair,
                Ok(None) => break,
                Err(err) => {
                    final_error = Some(err);
                    break;
                }
            };
            let limit = if options.fragment_length_samples != 0 && !sampler.finished() {
                PAIRED_SAMPLING_BUFFER_BYTES
            } else {
                0
            };
            match batch_builder.push(left, right, limit) {
                Ok(batches) => {
                    for batch in batches {
                        sequence += 1;
                        if work_tx.send((sequence, batch)).is_err() {
                            final_error = Some(Error::InvalidFormat("worker queue closed".into()));
                            break 'input;
                        }
                    }
                }
                Err(err) => {
                    final_error = Some(err);
                    break;
                }
            }
        }
        if final_error.is_none()
            && let Some(batch) = batch_builder.finish()
            && work_tx.send((sequence + 1, batch)).is_err()
        {
            final_error = Some(Error::InvalidFormat("worker queue closed".into()));
        }
        drop(work_tx);

        key_pool.clear();
        join_workers(handles, &mut merged, &mut merged_map, &mut final_error);
        timing::add_duration(Stage::Pseudoalign, pseudoalign_start.elapsed());
    });

    if let Some(err) = final_error {
        Err(err)
    } else {
        Ok(merged)
    }
}

#[cfg(test)]
mod paired_batch_tests {
    use super::*;
    use crate::io::FastqReader;
    use std::io::Cursor;

    fn packed(lengths: &[usize]) -> PackedSeqBatch {
        let mut fastq = String::new();
        for (i, &length) in lengths.iter().enumerate() {
            fastq.push_str(&format!(
                "@{i}\n{}\n+\n{}\n",
                "A".repeat(length),
                "I".repeat(length)
            ));
        }
        FastqReader::new(Cursor::new(fastq))
            .next_packed_seq_batch(lengths.len())
            .unwrap()
            .unwrap()
    }

    fn lengths(jobs: Vec<Vec<PairedChunk>>) -> Vec<Vec<(usize, usize)>> {
        jobs.into_iter()
            .map(|chunks| {
                chunks
                    .into_iter()
                    .flat_map(|chunk| {
                        chunk
                            .reads
                            .map(move |i| (chunk.left.seq(i).len(), chunk.right.seq(i).len()))
                    })
                    .collect()
            })
            .collect()
    }

    #[test]
    fn paired_sampling_buffers_span_parser_chunks_and_use_both_mate_lengths() {
        let mut builder = PairedBatchBuilder::default();
        assert!(
            builder
                .push(packed(&[10, 20]), packed(&[5, 7]), 100)
                .unwrap()
                .is_empty()
        );
        // First two pairs cost 17 + 29 = 46 bytes. The next costs 54,
        // so exactly 100 fails upstream's strict buffer bound.
        let mut jobs = builder
            .push(packed(&[30, 3]), packed(&[22, 4]), 100)
            .unwrap();
        jobs.push(builder.finish().unwrap());
        assert_eq!(
            lengths(jobs),
            vec![vec![(10, 5), (20, 7)], vec![(30, 22), (3, 4)]]
        );
    }

    #[test]
    fn paired_sampling_buffers_split_one_parser_chunk_without_losing_pairs() {
        let mut builder = PairedBatchBuilder::default();
        let mut jobs = builder
            .push(packed(&[10; 5]), packed(&[10; 5]), 50)
            .unwrap();
        jobs.push(builder.finish().unwrap());
        assert_eq!(
            lengths(jobs),
            vec![vec![(10, 10); 2], vec![(10, 10); 2], vec![(10, 10)]]
        );
        assert!(builder.finish().is_none());
    }

    #[test]
    fn ending_sampling_flushes_pending_reads_before_normal_work() {
        let mut builder = PairedBatchBuilder::default();
        assert!(
            builder
                .push(packed(&[10]), packed(&[12]), 100)
                .unwrap()
                .is_empty()
        );
        let jobs = builder.push(packed(&[20]), packed(&[22]), 0).unwrap();
        assert_eq!(lengths(jobs), vec![vec![(10, 12)], vec![(20, 22)]]);
        assert!(builder.finish().is_none());
    }

    #[test]
    fn oversized_pairs_progress_and_mate_count_errors_are_preserved() {
        let mut builder = PairedBatchBuilder::default();
        let mut jobs = builder
            .push(packed(&[100, 3]), packed(&[100, 4]), 100)
            .unwrap();
        jobs.push(builder.finish().unwrap());
        assert_eq!(lengths(jobs), vec![vec![(100, 100)], vec![(3, 4)]]);
        assert!(builder.push(packed(&[1]), packed(&[1, 2]), 100).is_err());
    }
}

#[cfg(test)]
mod ordered_sampler_tests {
    use super::OrderedSampler;

    /// Budget each buffer gets from the lag rule, computed sequentially.
    fn expected(goal: usize, lag: usize, observed: &[usize]) -> Vec<usize> {
        let mut prefix = vec![0];
        let mut budgets = Vec::new();
        for (index, &available) in observed.iter().enumerate() {
            let budget = if index < lag {
                goal
            } else {
                goal.saturating_sub(prefix[index + 1 - lag])
            };
            budgets.push(budget);
            prefix.push(prefix[index] + budget.min(available));
        }
        budgets
    }

    #[test]
    fn budgets_follow_the_lag_rule_for_any_completion_order() {
        let observed: Vec<usize> = (0..60).map(|i| 200 + (i * 37) % 120).collect();
        for workers in [1usize, 3, 8] {
            let want = expected(10_000, workers, &observed);
            let sampler = OrderedSampler::new(10_000, workers);
            let next = std::sync::atomic::AtomicUsize::new(0);
            let got: Vec<std::sync::atomic::AtomicUsize> =
                (0..observed.len()).map(|_| Default::default()).collect();
            std::thread::scope(|scope| {
                for worker in 0..workers {
                    let (sampler, next, got, observed) = (&sampler, &next, &got, &observed);
                    scope.spawn(move || {
                        loop {
                            let index = next.fetch_add(1, std::sync::atomic::Ordering::SeqCst);
                            let Some(&available) = observed.get(index) else {
                                break;
                            };
                            // Stagger completions so later buffers often finish first.
                            std::thread::sleep(std::time::Duration::from_micros(
                                ((index * 7 + worker * 13) % 5) as u64 * 50,
                            ));
                            let budget = sampler.budget(index as u64 + 1);
                            got[index].store(budget, std::sync::atomic::Ordering::SeqCst);
                            sampler.complete(index as u64 + 1, budget.min(available));
                        }
                    });
                }
            });
            let got: Vec<usize> = got
                .iter()
                .map(|value| value.load(std::sync::atomic::Ordering::SeqCst))
                .collect();
            assert_eq!(got, want, "workers={workers}");
            assert!(sampler.finished());
        }
    }

    #[test]
    fn an_earlier_buffer_finishing_late_keeps_its_budget() {
        let sampler = OrderedSampler::new(10, 1);
        assert_eq!(sampler.budget(1), 10);
        sampler.complete(1, 4);
        // Buffer 3 cannot start until buffer 2 completes, so compute buffer 2 first.
        assert_eq!(sampler.budget(2), 6);
        sampler.complete(2, 6);
        assert_eq!(sampler.budget(3), 0);
        // A late query for an earlier buffer still sees its real budget.
        let late = OrderedSampler::new(10, 2);
        late.complete(1, 10);
        assert_eq!(late.budget(3), 0);
        assert_eq!(late.budget(2), 10);
    }
}
