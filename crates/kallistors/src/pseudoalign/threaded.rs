use std::collections::HashMap;
use std::sync::mpsc;
use std::thread;

use crate::io::{FastqRecord, ReadSource, VecReadSource};
use crate::{Error, Result};

use super::{BifrostIndex, EcCounts, FragmentFilter, PseudoalignOptions, Strand};

pub fn pseudoalign_single_end_bifrost_with_options_threaded<R: ReadSource>(
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
    let mut read_error: Option<Error> = None;

    thread::scope(|scope| {
        let (res_tx, res_rx) = mpsc::channel::<(usize, Result<EcCounts>)>();
        let mut senders = Vec::with_capacity(threads);

        for tid in 0..threads {
            let (tx, rx) = mpsc::channel::<Vec<FastqRecord>>();
            senders.push(tx);
            let res_tx = res_tx.clone();
            scope.spawn(move || {
                let mut acc = EcCounts::new(options.bias);
                let mut acc_map: HashMap<Vec<u32>, usize> = HashMap::new();
                for batch in rx {
                    let mut batch_reader = VecReadSource::new(batch);
                    let res = super::single::pseudoalign_single_end_bifrost_with_options(
                        index,
                        &mut batch_reader,
                        strand,
                        filter,
                        options,
                    );
                    match res {
                        Ok(counts) => super::merge_ec_counts(&mut acc, &mut acc_map, counts),
                        Err(err) => {
                            let _ = res_tx.send((tid, Err(err)));
                            return;
                        }
                    }
                }
                let _ = res_tx.send((tid, Ok(acc)));
            });
        }

        let mut batch = Vec::with_capacity(super::BATCH_SIZE);
        let mut batch_idx = 0usize;
        while let Some(record) = reader.next_record() {
            match record {
                Ok(record) => {
                    batch.push(record);
                    if batch.len() >= super::BATCH_SIZE {
                        let to_send = std::mem::take(&mut batch);
                        if senders[batch_idx % threads].send(to_send).is_err() {
                            read_error = Some(Error::InvalidFormat("worker channel closed".into()));
                            break;
                        }
                        batch = Vec::with_capacity(super::BATCH_SIZE);
                        batch_idx += 1;
                    }
                }
                Err(err) => {
                    read_error = Some(err);
                    break;
                }
            }
        }
        if read_error.is_none() && !batch.is_empty() {
            let to_send = batch;
            if senders[batch_idx % threads].send(to_send).is_err() {
                read_error = Some(Error::InvalidFormat("worker channel closed".into()));
            }
        }
        drop(senders);

        let mut results: Vec<Option<Result<EcCounts>>> = (0..threads).map(|_| None).collect();
        for _ in 0..threads {
            if let Ok((tid, res)) = res_rx.recv() {
                results[tid] = Some(res);
            }
        }

        if read_error.is_none() && results.iter().any(|r| r.is_none()) {
            read_error = Some(Error::InvalidFormat("worker thread failed".into()));
        }

        for res in results.iter_mut().take(threads) {
            if let Some(res) = res.take() {
                match res {
                    Ok(counts) => super::merge_ec_counts(&mut merged, &mut merged_map, counts),
                    Err(err) => {
                        if read_error.is_none() {
                            read_error = Some(err);
                        }
                    }
                }
            }
        }
    });

    if let Some(err) = read_error {
        Err(err)
    } else {
        Ok(merged)
    }
}

pub fn pseudoalign_paired_bifrost_with_options_threaded<R1: ReadSource, R2: ReadSource>(
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
    let mut read_error: Option<Error> = None;

    thread::scope(|scope| {
        let (res_tx, res_rx) = mpsc::channel::<(usize, Result<EcCounts>)>();
        let mut senders = Vec::with_capacity(threads);

        for tid in 0..threads {
            let (tx, rx) = mpsc::channel::<(Vec<FastqRecord>, Vec<FastqRecord>)>();
            senders.push(tx);
            let res_tx = res_tx.clone();
            scope.spawn(move || {
                let mut acc = EcCounts::new(options.bias);
                let mut acc_map: HashMap<Vec<u32>, usize> = HashMap::new();
                for (left, right) in rx {
                    let mut left_reader = VecReadSource::new(left);
                    let mut right_reader = VecReadSource::new(right);
                    let res = super::paired::pseudoalign_paired_bifrost_with_options(
                        index,
                        &mut left_reader,
                        &mut right_reader,
                        strand,
                        options,
                    );
                    match res {
                        Ok(counts) => super::merge_ec_counts(&mut acc, &mut acc_map, counts),
                        Err(err) => {
                            let _ = res_tx.send((tid, Err(err)));
                            return;
                        }
                    }
                }
                let _ = res_tx.send((tid, Ok(acc)));
            });
        }

        let mut left_batch = Vec::with_capacity(super::BATCH_SIZE);
        let mut right_batch = Vec::with_capacity(super::BATCH_SIZE);
        let mut batch_idx = 0usize;
        loop {
            let r1 = reader1.next_record();
            let r2 = reader2.next_record();
            match (r1, r2) {
                (None, None) => break,
                (Some(_), None) | (None, Some(_)) => {
                    read_error = Some(Error::InvalidFormat("paired FASTQ length mismatch".into()));
                    break;
                }
                (Some(a), Some(b)) => match (a, b) {
                    (Ok(a), Ok(b)) => {
                        left_batch.push(a);
                        right_batch.push(b);
                        if left_batch.len() >= super::BATCH_SIZE {
                            let left_send = std::mem::take(&mut left_batch);
                            let right_send = std::mem::take(&mut right_batch);
                            if senders[batch_idx % threads]
                                .send((left_send, right_send))
                                .is_err()
                            {
                                read_error =
                                    Some(Error::InvalidFormat("worker channel closed".into()));
                                break;
                            }
                            left_batch = Vec::with_capacity(super::BATCH_SIZE);
                            right_batch = Vec::with_capacity(super::BATCH_SIZE);
                            batch_idx += 1;
                        }
                    }
                    (Err(err), _) | (_, Err(err)) => {
                        read_error = Some(err);
                        break;
                    }
                },
            }
        }
        if read_error.is_none() && !left_batch.is_empty() {
            let left_send = left_batch;
            let right_send = right_batch;
            if senders[batch_idx % threads]
                .send((left_send, right_send))
                .is_err()
            {
                read_error = Some(Error::InvalidFormat("worker channel closed".into()));
            }
        }
        drop(senders);

        let mut results: Vec<Option<Result<EcCounts>>> = (0..threads).map(|_| None).collect();
        for _ in 0..threads {
            if let Ok((tid, res)) = res_rx.recv() {
                results[tid] = Some(res);
            }
        }

        if read_error.is_none() && results.iter().any(|r| r.is_none()) {
            read_error = Some(Error::InvalidFormat("worker thread failed".into()));
        }

        for res in results.iter_mut().take(threads) {
            if let Some(res) = res.take() {
                match res {
                    Ok(counts) => super::merge_ec_counts(&mut merged, &mut merged_map, counts),
                    Err(err) => {
                        if read_error.is_none() {
                            read_error = Some(err);
                        }
                    }
                }
            }
        }
    });

    if let Some(err) = read_error {
        Err(err)
    } else {
        Ok(merged)
    }
}
