//! Bounded BAM decoding followed by ordered parallel BUS counting.

use super::bam_input::bam_string_tag;
use super::counting::{Mapper, PreparedRead, ReadTags, write_batch};
use super::*;

const BATCH_READS: usize = 4096;

pub(super) fn process<W: Write>(
    args: &BusArgs,
    index: &kallistors::pseudoalign::BifrostIndex,
    processing: BusProcessingConfig,
    state: &mut BusWriteState<'_, W>,
) -> Result<(u64, u64, u64, u64)> {
    let path = args
        .reads
        .first()
        .ok_or_else(|| anyhow!("--bam expects exactly one BAM file"))?;
    let mut reader = noodles_bam::io::Reader::new(File::open(path)?);
    reader
        .read_header()
        .map_err(|e| anyhow!("failed to read BAM header: {e}"))?;
    let mapper = Mapper::new(index, processing.options, args.threads)?;
    let mut counts = (0, 0, 0);
    let mut input = reader.records();
    loop {
        let mut prepared = Vec::with_capacity(BATCH_READS);
        while prepared.len() < BATCH_READS
            && processing
                .max_reads
                .is_none_or(|limit| counts.0 + (prepared.len() as u64) < limit)
        {
            let Some(record) = input.next() else { break };
            let record = record.map_err(|e| anyhow!("failed to read BAM record: {e}"))?;
            if record.flags().is_secondary() || record.flags().is_supplementary() {
                continue;
            }
            let barcode = bam_string_tag(&record, *b"CR")
                .or_else(|| bam_string_tag(&record, *b"CB").map(normalize_corrected_barcode));
            let umi = bam_string_tag(&record, *b"UR")
                .or_else(|| bam_string_tag(&record, *b"RX"))
                .or_else(|| bam_string_tag(&record, *b"MI").map(normalize_corrected_barcode))
                .or_else(|| bam_string_tag(&record, *b"UB").map(normalize_corrected_barcode));
            prepared.push(match (barcode, umi) {
                (Some(barcode), Some(umi)) => ReadTags::encode(&barcode, &umi).map(|tags| {
                    Some(PreparedRead {
                        batch_index: 0,
                        tags,
                        sequence: record.sequence().iter().collect(),
                        mate: None,
                        disable_strand: false,
                    })
                }),
                _ => Ok(None),
            });
        }
        if prepared.is_empty() {
            break;
        }
        // --num is deliberately ignored for BAM, matching the scalar reader.
        write_batch(mapper.map(prepared), false, &mut counts, state)?;
    }
    Ok((counts.0, counts.1, counts.2, 0))
}
