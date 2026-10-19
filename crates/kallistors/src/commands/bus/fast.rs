//! Bounded, ordered BUS counting for sequence-only single-cDNA technologies.

use std::thread::JoinHandle;

use anyhow::{Result, anyhow, bail};
use crossbeam_channel::{Receiver, bounded};
use kallistors::io::{FastqInput, FastqReader, PackedSeqBatch, PackedSeqBatchSource};
use kallistors::pseudoalign::BifrostIndex;

mod query_cache;
use query_cache::QueryCache;
use rayon::prelude::*;

use super::input::ReadGroup;
use super::{
    BusArgs, BusProcessingConfig, BusRecord, BusWriteState, SliceSpec, TechnologySpec, intern_ec,
    record_length, string_to_binary, write_bus_record,
};

const BATCH_READS: usize = 65_536;

pub(super) fn supports(args: &BusArgs, spec: &TechnologySpec) -> bool {
    !args.bam
        && !args.pseudobam
        && !args.genomebam
        && !args.long
        && !args.unmapped
        && !args.aa
        && !args.batch_barcodes
        && !spec.paired
        && !spec.keep_fastq_comments
        && spec.tag_len.is_none()
        && spec.bc.len() == 1
        && spec.umi.len() == 1
        && spec.seq.len() == 1
        && spec.bc[0].file.is_some()
        && spec.umi[0].file.is_some()
        && spec.seq[0].file.is_some()
}

type InputResult = kallistors::Result<Option<PackedSeqBatch>>;

enum Source {
    Inline(FastqReader<FastqInput>),
    Prefetched(Receiver<InputResult>),
}

struct Readers {
    sources: Vec<Source>,
    workers: Vec<JoinHandle<()>>,
    batch_records: usize,
    remaining_records: u64,
}

impl Readers {
    fn new(
        group: &ReadGroup,
        nfiles: usize,
        threaded: bool,
        max_reads: Option<u64>,
    ) -> Result<Self> {
        let multiplier = if group.interleaved { nfiles } else { 1 };
        let record_limit = max_reads
            .unwrap_or(u64::MAX)
            .saturating_mul(multiplier as u64);
        let mut result = Self {
            sources: Vec::new(),
            workers: Vec::new(),
            batch_records: BATCH_READS * multiplier,
            remaining_records: record_limit,
        };
        for path in &group.files {
            if !threaded {
                result
                    .sources
                    .push(Source::Inline(kallistors::io::open_fastq_reader(path)?));
                continue;
            }
            let (sender, receiver) = bounded(1);
            let path = path.clone();
            let batch_records = result.batch_records;
            result.workers.push(std::thread::spawn(move || {
                let mut reader = match kallistors::io::open_fastq_reader(&path) {
                    Ok(reader) => reader,
                    Err(error) => {
                        let _ = sender.send(Err(error));
                        return;
                    }
                };
                let mut remaining = record_limit;
                while remaining > 0 {
                    let count = batch_records.min(usize::try_from(remaining).unwrap_or(usize::MAX));
                    let batch = reader.next_packed_seq_batch(count);
                    if let Ok(Some(batch)) = &batch {
                        remaining = remaining.saturating_sub(batch.len() as u64);
                    }
                    let finished = remaining == 0 || !matches!(&batch, Ok(Some(_)));
                    if sender.send(batch).is_err() || finished {
                        break;
                    }
                }
            }));
            result.sources.push(Source::Prefetched(receiver));
        }
        Ok(result)
    }

    fn next(&mut self, interleaved: bool, nfiles: usize) -> Result<Option<Vec<PackedSeqBatch>>> {
        if self.remaining_records == 0 {
            return Ok(None);
        }
        let batch_records = self
            .batch_records
            .min(usize::try_from(self.remaining_records).unwrap_or(usize::MAX));
        let mut batches = Vec::with_capacity(self.sources.len());
        let mut ended = 0;
        for source in &mut self.sources {
            let next = match source {
                Source::Inline(reader) => reader.next_packed_seq_batch(batch_records),
                Source::Prefetched(receiver) => receiver
                    .recv()
                    .map_err(|_| anyhow!("FASTQ reader stopped before completing its input"))?,
            }
            .map_err(|error| anyhow!("FASTQ read failed: {error}"))?;
            match next {
                Some(batch) => batches.push(batch),
                None => ended += 1,
            }
        }
        if ended == self.sources.len() {
            return Ok(None);
        }
        if ended > 0 || batches.iter().any(|batch| batch.len() != batches[0].len()) {
            bail!("FASTQ files ended at different record counts");
        }
        if interleaved && !batches[0].len().is_multiple_of(nfiles) {
            bail!("interleaved FASTQ ended in the middle of a record group");
        }
        self.remaining_records = self
            .remaining_records
            .saturating_sub(batches[0].len() as u64);
        Ok(Some(batches))
    }
}

impl Drop for Readers {
    fn drop(&mut self) {
        // Disconnect bounded queues before joining, including on errors/--num-reads.
        self.sources.clear();
        for worker in self.workers.drain(..) {
            let _ = worker.join();
        }
    }
}

struct CountedRead<'index> {
    barcode: u64,
    umi: u64,
    flags: u32,
    barcode_len: usize,
    umi_len: usize,
    ec: Option<std::borrow::Cow<'index, [u32]>>,
}

fn sequence(
    batches: &[PackedSeqBatch],
    read: usize,
    slice: SliceSpec,
    interleaved: bool,
    nfiles: usize,
) -> Result<&[u8]> {
    let file = slice.file.expect("validated sequence-only BUS layout");
    if file >= nfiles {
        bail!("technology references missing file {}", file + 1);
    }
    Ok(if interleaved {
        batches[0].seq(read * nfiles + file)
    } else {
        batches[file].seq(read)
    })
}

fn sliced(seq: &[u8], slice: SliceSpec) -> Option<&[u8]> {
    seq.get(slice.start..slice.stop.unwrap_or(seq.len()))
}

fn classify<'index>(
    batches: &[PackedSeqBatch],
    read: usize,
    interleaved: bool,
    spec: &TechnologySpec,
    cache: &QueryCache<'index>,
) -> Result<Option<CountedRead<'index>>> {
    let barcode = sequence(batches, read, spec.bc[0], interleaved, spec.nfiles)?;
    let Some(barcode) = sliced(barcode, spec.bc[0]) else {
        return Ok(None);
    };
    let umi = sequence(batches, read, spec.umi[0], interleaved, spec.nfiles)?;
    let Some(umi) = sliced(umi, spec.umi[0]) else {
        return Ok(None);
    };
    let query = sequence(batches, read, spec.seq[0], interleaved, spec.nfiles)?;
    let query = sliced(query, spec.seq[0]).ok_or_else(|| {
        anyhow!(
            "read length {} is too short for slice {}:{}",
            query.len(),
            spec.seq[0].start,
            spec.seq[0].stop.unwrap_or(query.len())
        )
    })?;
    let mut barcode_flags = 0;
    let mut umi_flags = 0;
    Ok(Some(CountedRead {
        barcode: string_to_binary(barcode, &mut barcode_flags, "barcode")?,
        umi: string_to_binary(umi, &mut umi_flags, "UMI")?,
        flags: barcode_flags | (umi_flags << 8),
        barcode_len: barcode.len(),
        umi_len: umi.len(),
        ec: cache.lookup(query),
    }))
}

fn produce_batches<'index>(
    args: &BusArgs,
    spec: &TechnologySpec,
    index: &'index BifrostIndex,
    groups: &[ReadGroup],
    processing: BusProcessingConfig,
    mut emit: impl FnMut(Vec<Result<Option<CountedRead<'index>>>>) -> Result<bool>,
) -> Result<()> {
    let pool = if args.threads > 1 {
        Some(
            rayon::ThreadPoolBuilder::new()
                .num_threads(args.threads)
                .build()?,
        )
    } else {
        None
    };
    let cache = QueryCache::new(index, processing.options);
    let mut processed = 0u64;
    'groups: for group in groups {
        let remaining = processing
            .max_reads
            .map(|limit| limit.saturating_sub(processed));
        let mut readers = Readers::new(group, spec.nfiles, pool.is_some(), remaining)?;
        while let Some(batches) = readers.next(group.interleaved, spec.nfiles)? {
            let available = batches[0].len() / if group.interleaved { spec.nfiles } else { 1 };
            let remaining = processing
                .max_reads
                .map_or(u64::MAX, |limit| limit.saturating_sub(processed));
            let count = available.min(usize::try_from(remaining).unwrap_or(usize::MAX));
            let map = |read| classify(&batches, read, group.interleaved, spec, &cache);
            let results: Vec<Result<Option<CountedRead>>> = if let Some(pool) = &pool {
                pool.install(|| {
                    (0..count)
                        .into_par_iter()
                        .with_max_len(256)
                        .map(map)
                        .collect()
                })
            } else {
                (0..count).map(map).collect()
            };
            if !emit(results)? {
                break 'groups;
            }
            processed += count as u64;
            if processing.max_reads.is_some_and(|limit| processed >= limit) {
                break 'groups;
            }
        }
    }
    Ok(())
}

pub(super) fn process<W: std::io::Write>(
    args: &BusArgs,
    spec: &TechnologySpec,
    index: &BifrostIndex,
    groups: &[ReadGroup],
    processing: BusProcessingConfig,
    state: &mut BusWriteState<'_, W>,
) -> Result<(u64, u64, u64)> {
    let (mut processed, mut aligned, mut unique) = (0u64, 0u64, 0u64);
    let mut consume = |results: Vec<Result<Option<CountedRead>>>| -> Result<()> {
        for result in results {
            let number = processed;
            processed += 1;
            let Some(read) = result? else { continue };
            record_length(state.bc_len_hist, read.barcode_len);
            record_length(state.umi_len_hist, read.umi_len);
            let Some(ec) = read.ec else { continue };
            aligned += 1;
            unique += u64::from(ec.len() == 1);
            let ec = if let Some(&id) = state.ec_map.get(ec.as_ref()) {
                id
            } else {
                intern_ec(state.ec_map, state.ec_list, ec.into_owned())
            };
            write_bus_record(
                state.writer,
                &BusRecord {
                    barcode: read.barcode,
                    umi: read.umi,
                    ec,
                    count: 1,
                    flags: if args.num {
                        u32::try_from(number).unwrap_or(u32::MAX)
                    } else {
                        read.flags
                    },
                },
            )?;
        }
        Ok(())
    };
    if args.threads > 1 {
        std::thread::scope(|scope| -> Result<()> {
            let (sender, receiver) = bounded(1);
            let mapper = scope.spawn(move || {
                if let Err(error) = produce_batches(args, spec, index, groups, processing, |rows| {
                    Ok(sender.send(Ok(rows)).is_ok())
                }) {
                    let _ = sender.send(Err(error));
                }
            });
            let result = (|| -> Result<()> {
                for batch in &receiver {
                    consume(batch?)?;
                }
                Ok(())
            })();
            // Drop first: a failed output must not leave the producer blocked.
            drop(receiver);
            let joined = mapper.join();
            result?;
            joined.map_err(|_| anyhow!("BUS mapping worker panicked"))?;
            Ok(())
        })?;
    } else {
        produce_batches(args, spec, index, groups, processing, |rows| {
            consume(rows)?;
            Ok(true)
        })?;
    }
    Ok((processed, aligned, unique))
}
