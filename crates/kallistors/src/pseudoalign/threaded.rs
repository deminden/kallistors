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
    // Bias sampling retains the complete-record route, which owns its budget.
    if threads <= 1 && options.bias {
        return super::single::pseudoalign_single_end_bifrost_with_options(
            index, reader, strand, filter, options,
        );
    }
    let cache_entries = options
        .match_cache_size
        .allocation(threads)?
        .entries_per_worker();
    if threads <= 1 {
        let _cache_budget = super::core::MatchCacheBudgetScope::enter(cache_entries)?;
        // Ordinary synchronous counting consumes only sequences. Reuse the
        // threaded parser and mapper without spawning a producer or worker.
        let mut counts = EcCounts::new(false);
        let mut map = EcCountMap::default();
        while let Some(batch) = reader.next_packed_seq_batch(4_096)? {
            super::single::pseudoalign_single_end_bifrost_batch_into(
                index,
                batch,
                strand,
                filter,
                options,
                &mut counts,
                &mut map,
            );
        }
        return Ok(counts);
    }

    let threads = threads.max(1);
    let mut merged = EcCounts::new(options.bias);
    let mut merged_map: EcCountMap = EcCountMap::default();
    let mut final_error: Option<Error> = None;
    let key_pool = Arc::new(super::utils::SharedEcPool::default());

    thread::scope(|scope| {
        // One prefetched batch overlaps decoding; rendezvous dispatch bounds read
        // storage independently of worker count without dropping or reordering input.
        let (producer_tx, producer_rx) = bounded::<PackedSeqBatch>(1);
        let (work_tx, work_rx) = bounded::<PackedSeqBatch>(0);
        let mut handles = Vec::with_capacity(threads);

        let read_start = Instant::now();
        let producer = scope.spawn(move || produce_batches(reader, &producer_tx, 4_096));

        let pseudoalign_start = Instant::now();
        for _ in 0..threads {
            let rx = work_rx.clone();
            let key_pool = Arc::clone(&key_pool);
            handles.push(
                scope.spawn(move || -> Result<(EcCounts, super::utils::SharedEcKeys)> {
                    let _cache_budget = super::core::MatchCacheBudgetScope::enter(cache_entries)?;
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
                        );
                    }
                    Ok((acc, acc_map.into_shared_keys()))
                }),
            );
        }

        drop(work_rx);
        while let Ok(batch) = producer_rx.recv() {
            if work_tx.send(batch).is_err() {
                final_error = Some(Error::InvalidFormat("worker queue closed".into()));
                break;
            }
        }
        drop(work_tx);
        drop(producer_rx);

        match producer.join() {
            Ok(Ok(())) => {}
            Ok(Err(err)) => {
                if final_error.is_none() {
                    final_error = Some(err);
                }
            }
            Err(_) => {
                if final_error.is_none() {
                    final_error = Some(Error::InvalidFormat("producer thread panicked".into()));
                }
            }
        }
        timing::add_duration(Stage::FastqReadDecompress, read_start.elapsed());

        key_pool.clear();
        for handle in handles {
            match handle.join() {
                Ok(Ok((counts, keys))) => {
                    let merge_start = Instant::now();
                    super::utils::merge_shared_ec_counts(
                        &mut merged,
                        &mut merged_map,
                        counts,
                        keys,
                    );
                    timing::add_duration(Stage::EcMerge, merge_start.elapsed());
                }
                Ok(Err(err)) => {
                    if final_error.is_none() {
                        final_error = Some(err);
                    }
                }
                Err(_) => {
                    if final_error.is_none() {
                        final_error = Some(Error::InvalidFormat("worker thread panicked".into()));
                    }
                }
            }
        }
        timing::add_duration(Stage::Pseudoalign, pseudoalign_start.elapsed());
    });

    if let Some(err) = final_error {
        Err(err)
    } else {
        Ok(merged)
    }
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
    let cache_entries = options
        .match_cache_size
        .allocation(threads)?
        .entries_per_worker();
    if threads <= 1 {
        let _cache_budget = super::core::MatchCacheBudgetScope::enter(cache_entries)?;
        // Ordinary mapping parses sequence-only batches without workers.
        let mut counts = EcCounts::new(options.bias);
        let mut ec_map = EcCountMap::default();
        let sampler = super::paired::FragmentLengthSampler::new(options.fragment_length_samples);
        loop {
            let left = reader1.next_packed_seq_batch(super::BATCH_SIZE / 2)?;
            let right = reader2.next_packed_seq_batch(super::BATCH_SIZE / 2)?;
            match (left, right) {
                (None, None) => return Ok(counts),
                (Some(left), Some(right)) => {
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
                    )?;
                }
                _ => return Err(Error::InvalidFormat("paired FASTQ length mismatch".into())),
            }
        }
    }

    let threads = threads.max(1);
    let mut merged = EcCounts::new(options.bias);
    let mut merged_map: EcCountMap = EcCountMap::default();
    let mut final_error: Option<Error> = None;

    let key_pool = Arc::new(super::utils::SharedEcPool::default());
    let sampler = super::paired::FragmentLengthSampler::new(options.fragment_length_samples);
    thread::scope(|scope| {
        let (left_tx, left_rx) = bounded::<PackedSeqBatch>(1);
        let (right_tx, right_rx) = bounded::<PackedSeqBatch>(1);
        let (work_tx, work_rx) = bounded::<Vec<PairedChunk>>(0);
        let mut handles = Vec::with_capacity(threads);

        let read_start = Instant::now();
        let left_producer =
            scope.spawn(move || produce_batches(reader1, &left_tx, super::BATCH_SIZE / 2));
        let right_producer =
            scope.spawn(move || produce_batches(reader2, &right_tx, super::BATCH_SIZE / 2));

        let pseudoalign_start = Instant::now();
        for _ in 0..threads {
            let rx = work_rx.clone();
            let sampler = &sampler;
            let key_pool = Arc::clone(&key_pool);
            handles.push(
                scope.spawn(move || -> Result<(EcCounts, super::utils::SharedEcKeys)> {
                    let _cache_budget = super::core::MatchCacheBudgetScope::enter(cache_entries)?;
                    let mut acc = EcCounts::new(options.bias);
                    let mut acc_map = EcCountMap::shared(key_pool);
                    while let Ok(chunks) = rx.recv() {
                        let budget = sampler.remaining();
                        let local_sampler = super::paired::FragmentLengthSampler::new(budget);
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
                            )?;
                        }
                        if options.fragment_length_samples != 0 {
                            sampler.submit(budget - local_sampler.remaining());
                        }
                    }
                    Ok((acc, acc_map.into_shared_keys()))
                }),
            );
        }

        let mut batch_builder = PairedBatchBuilder::default();
        loop {
            let left = left_rx.recv().ok();
            let right = right_rx.recv().ok();
            match (left, right) {
                (Some(left_batch), Some(right_batch)) => {
                    let limit = if options.fragment_length_samples != 0 && sampler.remaining() > 0 {
                        PAIRED_SAMPLING_BUFFER_BYTES
                    } else {
                        0
                    };
                    match batch_builder.push(left_batch, right_batch, limit) {
                        Ok(batches) => {
                            for batch in batches {
                                if work_tx.send(batch).is_err() {
                                    final_error =
                                        Some(Error::InvalidFormat("worker queue closed".into()));
                                    break;
                                }
                            }
                            if final_error.is_some() {
                                break;
                            }
                        }
                        Err(err) => {
                            final_error = Some(err);
                            break;
                        }
                    }
                }
                (None, None) => break,
                (None, Some(_)) | (Some(_), None) => {
                    final_error = Some(Error::InvalidFormat("paired FASTQ length mismatch".into()));
                    break;
                }
            }
        }
        if final_error.is_none()
            && let Some(batch) = batch_builder.finish()
            && work_tx.send(batch).is_err()
        {
            final_error = Some(Error::InvalidFormat("worker queue closed".into()));
        }
        drop(work_tx);
        // Release blocked senders before joining producers after an input error.
        drop(left_rx);
        drop(right_rx);

        for producer in [left_producer, right_producer] {
            match producer.join() {
                Ok(Ok(())) => {}
                Ok(Err(err)) => {
                    if final_error.is_none() {
                        final_error = Some(err);
                    }
                }
                Err(_) => {
                    if final_error.is_none() {
                        final_error = Some(Error::InvalidFormat("producer thread panicked".into()));
                    }
                }
            }
        }
        timing::add_duration(Stage::FastqReadDecompress, read_start.elapsed());

        key_pool.clear();
        for handle in handles {
            match handle.join() {
                Ok(Ok((counts, keys))) => {
                    let merge_start = Instant::now();
                    super::utils::merge_shared_ec_counts(
                        &mut merged,
                        &mut merged_map,
                        counts,
                        keys,
                    );
                    timing::add_duration(Stage::EcMerge, merge_start.elapsed());
                }
                Ok(Err(err)) => {
                    if final_error.is_none() {
                        final_error = Some(err);
                    }
                }
                Err(_) => {
                    if final_error.is_none() {
                        final_error = Some(Error::InvalidFormat("worker thread panicked".into()));
                    }
                }
            }
        }
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
