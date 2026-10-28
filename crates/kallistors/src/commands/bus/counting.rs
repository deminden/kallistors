//! Shared ordered BUS counting for packed FASTQ, normalized layouts and BAM.

use super::fast::QueryCache;
use super::{
    BusArgs, BusRecord, BusWriteState, intern_ec, record_length, record_paired_fragment_length,
    string_to_binary, write_bus_record,
};
use anyhow::Result;
use kallistors::pseudoalign::{BifrostIndex, PseudoalignOptions};
use rayon::prelude::*;
use std::borrow::Cow;

pub(super) fn supports(args: &BusArgs) -> bool {
    !args.pseudobam && !args.genomebam && !args.long && !args.unmapped && !args.aa
}

pub(super) struct ReadTags {
    pub barcode: u64,
    pub umi: u64,
    pub flags: u32,
    pub barcode_len: usize,
    pub umi_len: usize,
}

impl ReadTags {
    pub(super) fn encode(barcode: &[u8], umi: &[u8]) -> Result<Self> {
        let (mut bc_flags, mut umi_flags) = (0, 0);
        Ok(Self {
            barcode: string_to_binary(barcode, &mut bc_flags, "barcode")?,
            umi: string_to_binary(umi, &mut umi_flags, "UMI")?,
            flags: bc_flags | (umi_flags << 8),
            barcode_len: barcode.len(),
            umi_len: umi.len(),
        })
    }
}

pub(super) struct CountedRead<'a, Fragment = Option<(i64, usize)>> {
    pub tags: ReadTags,
    pub ec: Option<Cow<'a, [u32]>>,
    pub fragment_length: Fragment,
}

pub(super) struct PreparedRead {
    pub batch_index: usize,
    pub tags: ReadTags,
    pub sequence: Vec<u8>,
    pub mate: Option<Vec<u8>>,
    pub disable_strand: bool,
}

pub(super) struct Mapper<'a> {
    index: &'a BifrostIndex,
    options: PseudoalignOptions,
    cache: QueryCache<'a>,
    // Tag-free reads can disable strand filtering. Keep a distinct answer cache
    // because equal cDNA bytes need not produce equal ECs under both policies.
    unstranded_cache: QueryCache<'a>,
    pool: Option<rayon::ThreadPool>,
}

impl<'a> Mapper<'a> {
    pub(super) fn new(
        index: &'a BifrostIndex,
        options: PseudoalignOptions,
        threads: usize,
    ) -> Result<Self> {
        let unstranded = PseudoalignOptions {
            strand_specific: None,
            ..options
        };
        Ok(Self {
            index,
            options,
            cache: QueryCache::new(index, options),
            unstranded_cache: QueryCache::new(index, unstranded),
            pool: if threads > 1 {
                Some(
                    rayon::ThreadPoolBuilder::new()
                        .num_threads(threads)
                        .build()?,
                )
            } else {
                None
            },
        })
    }

    pub(super) fn map(
        &self,
        reads: Vec<Result<Option<PreparedRead>>>,
    ) -> Vec<Result<Option<CountedRead<'a>>>> {
        let classify = |read: Result<Option<PreparedRead>>| -> Result<Option<CountedRead<'a>>> {
            let Some(read) = read? else { return Ok(None) };
            let (ec, fragment_length) = if let Some(mate) = &read.mate {
                let options = if read.disable_strand {
                    PseudoalignOptions {
                        strand_specific: None,
                        ..self.options
                    }
                } else {
                    self.options
                };
                match kallistors::pseudoalign::ec_for_pair_bifrost_counting(
                    self.index,
                    &read.sequence,
                    mate,
                    kallistors::pseudoalign::Strand::Unstranded,
                    options,
                ) {
                    Some(mapped) => (Some(Cow::Owned(mapped.ec)), mapped.fragment_length),
                    None => (None, None),
                }
            } else {
                let cache = if read.disable_strand && self.options.strand_specific.is_some() {
                    &self.unstranded_cache
                } else {
                    &self.cache
                };
                (cache.lookup(&read.sequence), None)
            };
            Ok(Some(CountedRead {
                tags: read.tags,
                ec,
                fragment_length: fragment_length.map(|length| (length, read.batch_index)),
            }))
        };
        if let Some(pool) = &self.pool {
            pool.install(|| {
                reads
                    .into_par_iter()
                    .with_max_len(256)
                    .map(classify)
                    .collect()
            })
        } else {
            reads.into_iter().map(classify).collect()
        }
    }
}

pub(super) trait FragmentObservation {
    fn observation(self) -> Option<(i64, usize)>;
}
impl FragmentObservation for () {
    fn observation(self) -> Option<(i64, usize)> {
        None
    }
}
impl FragmentObservation for Option<(i64, usize)> {
    fn observation(self) -> Option<(i64, usize)> {
        self
    }
}

pub(super) fn write_batch<W: std::io::Write, F: FragmentObservation>(
    rows: Vec<Result<Option<CountedRead<'_, F>>>>,
    read_numbers: bool,
    counts: &mut (u64, u64, u64),
    state: &mut BusWriteState<'_, W>,
) -> Result<()> {
    // Indexed parallel collection preserves input order. Intern ECs here so
    // worker completion order cannot change BUS IDs or --num read numbers.
    for row in rows {
        let number = counts.0;
        counts.0 += 1;
        let Some(read) = row? else { continue };
        record_length(state.bc_len_hist, read.tags.barcode_len);
        record_length(state.umi_len_hist, read.tags.umi_len);
        let Some(ec) = read.ec else { continue };
        counts.1 += 1;
        counts.2 += u64::from(ec.len() == 1);
        if let Some((length, batch_index)) = read.fragment_length.observation() {
            record_paired_fragment_length(state.paired_flens, &ec, length);
            if let Some(hist) = state.batch_paired_flens.get_mut(batch_index) {
                record_paired_fragment_length(hist, &ec, length);
            }
        }
        let ec = if let Some(&id) = state.ec_map.get(ec.as_ref()) {
            id
        } else {
            intern_ec(state.ec_map, state.ec_list, ec.into_owned())
        };
        write_bus_record(
            state.writer,
            &BusRecord {
                barcode: read.tags.barcode,
                umi: read.tags.umi,
                ec,
                count: 1,
                flags: if read_numbers {
                    u32::try_from(number).unwrap_or(u32::MAX)
                } else {
                    read.tags.flags
                },
            },
        )?;
    }
    Ok(())
}
