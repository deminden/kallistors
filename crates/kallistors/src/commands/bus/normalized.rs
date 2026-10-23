//! Bounded parallel counting for compound, synthetic and tagged FASTQ layouts.

use super::counting::{Mapper, PreparedRead, ReadTags, write_batch};
use super::*;

const BATCH_READS: usize = 4096;

pub(super) fn supports(args: &BusArgs) -> bool {
    !args.bam && counting::supports(args)
}

fn prepare_read(
    batch: &[kallistors::io::FastqRecord],
    args: &BusArgs,
    spec: &TechnologySpec,
    tag: Option<&TagConfig>,
    batch_index: usize,
) -> Result<Option<PreparedRead>> {
    let Some(barcode) = barcode_sequence(batch, spec, batch_index, args.batch_barcodes)? else {
        return Ok(None);
    };
    let Some(umi) = umi_value(batch, spec, tag)? else {
        return Ok(None);
    };
    let mut flags = 0;
    let barcode_binary = string_to_binary(&barcode, &mut flags, "barcode")?;
    let tags = ReadTags {
        barcode: barcode_binary,
        umi: umi.binary,
        flags: flags | (umi.flag << 8),
        barcode_len: barcode.len(),
        umi_len: umi.len,
    };
    let (sequence, mate) = if spec.paired {
        if spec.seq.len() != 2 {
            bail!("paired BUS technology requires exactly two sequence slices")
        }
        (
            extract_sequence_slice(batch, spec.seq[0], spec, umi.tag_present)?,
            Some(extract_sequence_slice(
                batch,
                spec.seq[1],
                spec,
                umi.tag_present,
            )?),
        )
    } else {
        (concat_sequence_slices(batch, spec, umi.tag_present)?, None)
    };
    Ok(Some(PreparedRead {
        batch_index,
        tags,
        sequence,
        mate,
        disable_strand: umi.disable_strand_specific,
    }))
}

pub(super) fn process<W: Write>(
    args: &BusArgs,
    spec: &TechnologySpec,
    index: &kallistors::pseudoalign::BifrostIndex,
    groups: &[input::ReadGroup],
    processing: BusProcessingConfig,
    tag: Option<&TagConfig>,
    state: &mut BusWriteState<'_, W>,
) -> Result<(u64, u64, u64)> {
    let mapper = Mapper::new(index, processing.options, args.threads)?;
    let mut counts = (0, 0, 0);
    for group in groups {
        let mut readers = group
            .files
            .iter()
            .map(|path| kallistors::io::open_fastq_reader(path))
            .collect::<kallistors::Result<Vec<_>>>()
            .map_err(|e| anyhow!("failed to open FASTQ: {e}"))?;
        loop {
            let mut prepared = Vec::with_capacity(BATCH_READS);
            while prepared.len() < BATCH_READS
                && processing
                    .max_reads
                    .is_none_or(|limit| counts.0 + (prepared.len() as u64) < limit)
            {
                let batch = if group.interleaved {
                    next_interleaved_batch(&mut readers, spec.nfiles)?
                } else {
                    next_parallel_batch(&mut readers)?
                };
                if batch.is_empty() {
                    break;
                }
                prepared.push(prepare_read(&batch, args, spec, tag, group.batch_index));
            }
            if prepared.is_empty() {
                break;
            }
            write_batch(mapper.map(prepared), args.num, &mut counts, state)?;
        }
        if processing.max_reads.is_some_and(|limit| counts.0 >= limit) {
            break;
        }
    }
    Ok(counts)
}
