use std::collections::HashMap;
use std::thread;
use std::time::Instant;

use crossbeam_channel::{Receiver, Sender, bounded};

use crate::io::{PackedBatchSource, PackedFastqBatch, ReadSource};
use crate::timing::{self, Stage};
use crate::{Error, Result};

use super::{BifrostIndex, EcCounts, FragmentFilter, PseudoalignOptions, Strand};

fn produce_batches<R: PackedBatchSource>(
    reader: &mut R,
    tx: &Sender<Result<PackedFastqBatch>>,
) -> Result<()> {
    while let Some(batch) = reader.next_packed_batch(super::BATCH_SIZE)? {
        if tx.send(Ok(batch)).is_err() {
            return Ok(());
        }
    }
    Ok(())
}

fn recv_batch(
    rx: &Receiver<Result<PackedFastqBatch>>,
) -> std::result::Result<Option<PackedFastqBatch>, Error> {
    match rx.recv() {
        Ok(Ok(batch)) => Ok(Some(batch)),
        Ok(Err(err)) => Err(err),
        Err(_) => Ok(None),
    }
}

pub fn pseudoalign_single_end_bifrost_with_options_threaded<
    R: ReadSource + PackedBatchSource + Send,
>(
    index: &BifrostIndex,
    reader: &mut R,
    strand: Strand,
    filter: Option<FragmentFilter>,
    options: PseudoalignOptions,
    threads: usize,
) -> Result<EcCounts> {
    if threads <= 1 {
        return super::single::pseudoalign_single_end_bifrost_with_options(
            index, reader, strand, filter, options,
        );
    }

    let threads = threads.max(1);
    let mut merged = EcCounts::new(options.bias);
    let mut merged_map: HashMap<Vec<u32>, usize> = HashMap::new();
    let mut final_error: Option<Error> = None;

    thread::scope(|scope| {
        let (producer_tx, producer_rx) = bounded::<Result<PackedFastqBatch>>(threads * 2);
        let (work_tx, work_rx) = bounded::<PackedFastqBatch>(threads * 2);
        let mut handles = Vec::with_capacity(threads);

        let read_start = Instant::now();
        let producer = scope.spawn(move || produce_batches(reader, &producer_tx));

        let pseudoalign_start = Instant::now();
        for _ in 0..threads {
            let rx = work_rx.clone();
            handles.push(scope.spawn(move || -> Result<EcCounts> {
                let mut acc = EcCounts::new(options.bias);
                let mut acc_map: HashMap<Vec<u32>, usize> = HashMap::new();
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
                Ok(acc)
            }));
        }

        loop {
            match recv_batch(&producer_rx) {
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

        for handle in handles {
            match handle.join() {
                Ok(Ok(counts)) => {
                    let merge_start = Instant::now();
                    super::merge_ec_counts(&mut merged, &mut merged_map, counts);
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
    R1: ReadSource + PackedBatchSource + Send,
    R2: ReadSource + PackedBatchSource + Send,
>(
    index: &BifrostIndex,
    reader1: &mut R1,
    reader2: &mut R2,
    strand: Strand,
    options: PseudoalignOptions,
    threads: usize,
) -> Result<EcCounts> {
    if threads <= 1 {
        return super::paired::pseudoalign_paired_bifrost_with_options(
            index, reader1, reader2, strand, options,
        );
    }

    let threads = threads.max(1);
    let mut merged = EcCounts::new(options.bias);
    let mut merged_map: HashMap<Vec<u32>, usize> = HashMap::new();
    let mut final_error: Option<Error> = None;

    thread::scope(|scope| {
        let (left_tx, left_rx) = bounded::<Result<PackedFastqBatch>>(threads * 2);
        let (right_tx, right_rx) = bounded::<Result<PackedFastqBatch>>(threads * 2);
        let (work_tx, work_rx) = bounded::<(PackedFastqBatch, PackedFastqBatch)>(threads * 2);
        let mut handles = Vec::with_capacity(threads);

        let read_start = Instant::now();
        let left_producer = scope.spawn(move || produce_batches(reader1, &left_tx));
        let right_producer = scope.spawn(move || produce_batches(reader2, &right_tx));

        let pseudoalign_start = Instant::now();
        for _ in 0..threads {
            let rx = work_rx.clone();
            handles.push(scope.spawn(move || -> Result<EcCounts> {
                let mut acc = EcCounts::new(options.bias);
                let mut acc_map: HashMap<Vec<u32>, usize> = HashMap::new();
                while let Ok((left, right)) = rx.recv() {
                    super::paired::pseudoalign_paired_bifrost_batch_into(
                        index,
                        left,
                        right,
                        strand,
                        options,
                        &mut acc,
                        &mut acc_map,
                    )?;
                }
                Ok(acc)
            }));
        }

        loop {
            let left = recv_batch(&left_rx);
            let right = recv_batch(&right_rx);
            match (left, right) {
                (Ok(Some(left_batch)), Ok(Some(right_batch))) => {
                    if work_tx.send((left_batch, right_batch)).is_err() {
                        final_error = Some(Error::InvalidFormat("worker queue closed".into()));
                        break;
                    }
                }
                (Ok(None), Ok(None)) => break,
                (Ok(None), Ok(Some(_))) | (Ok(Some(_)), Ok(None)) => {
                    final_error = Some(Error::InvalidFormat("paired FASTQ length mismatch".into()));
                    break;
                }
                (Err(err), _) | (_, Err(err)) => {
                    final_error = Some(err);
                    break;
                }
            }
        }
        drop(work_tx);

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

        for handle in handles {
            match handle.join() {
                Ok(Ok(counts)) => {
                    let merge_start = Instant::now();
                    super::merge_ec_counts(&mut merged, &mut merged_map, counts);
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
