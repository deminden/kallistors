//! Transcriptome/genome BAM output, mate flags and coordinate sorting.

use anyhow::{Result, anyhow};
use noodles_core::Position;
use noodles_sam::Header;
use noodles_sam::alignment::{
    RecordBuf,
    io::Write as _,
    record::{
        Cigar as _, Flags, MappingQuality,
        cigar::{Op, op::Kind},
        data::field::Tag,
    },
    record_buf::{Cigar, QualityScores, Sequence, data::field::Value as SamBufDataValue},
};
use noodles_sam::header::record::value::{
    Map,
    map::{
        Header as SamHeaderRecord, ReferenceSequence,
        header::{sort_order, tag::SORT_ORDER},
    },
};
use std::fs::File;
use std::num::NonZero;
use std::path::{Path, PathBuf};

use super::genome::{GenomeAlignment, GenomeModel};
use super::technology::{
    TechnologySpec, concat_quality_slices, concat_sequence_slices, extract_quality_slice,
    extract_sequence_slice,
};

pub(super) struct PseudoBamWriter {
    writer: noodles_bam::io::Writer<noodles_bgzf::io::Writer<File>>,
    header: Header,
    genome: Option<GenomeModel>,
    pending: Vec<(BamSortKey, RecordBuf)>,
    path: PathBuf,
}

#[derive(Clone, Copy, Eq, Ord, PartialEq, PartialOrd)]
struct BamSortKey {
    reference_sequence_id: usize,
    alignment_start: usize,
    input_order: u64,
}

struct PseudoBamRead<'a> {
    name: Vec<u8>,
    flags: Flags,
    mate_reference_sequence_id: Option<usize>,
    mate_alignment_start: Option<usize>,
    mate_reverse: bool,
    template_length: i32,
    ec_id: i32,
    ec_list: &'a [Vec<u32>],
    placement: Option<kallistors::pseudoalign::PseudoalignPlacement>,
    barcode: &'a [u8],
    umi: &'a [u8],
    sequence: &'a [u8],
    quality: &'a [u8],
    fallback_to_ec_transcript: bool,
}

pub(super) struct PseudoBamBusRead<'a> {
    pub(super) read_id: u64,
    pub(super) ec_id: i32,
    pub(super) ec_list: &'a [Vec<u32>],
    pub(super) batch: &'a [kallistors::io::FastqRecord],
    pub(super) spec: &'a TechnologySpec,
    pub(super) tag_present: bool,
    pub(super) barcode: &'a [u8],
    pub(super) umi: &'a [u8],
    pub(super) placement: Option<kallistors::pseudoalign::PseudoalignPlacement>,
    pub(super) mate_placements: &'a [Option<kallistors::pseudoalign::PseudoalignPlacement>],
    pub(super) mate_has_ec: &'a [bool],
}

impl PseudoBamWriter {
    pub(super) fn create(
        out_dir: &Path,
        index: &kallistors::pseudoalign::BifrostIndex,
        genome: Option<GenomeModel>,
    ) -> Result<Self> {
        let mut builder = Header::builder();
        if genome.is_some() {
            let sam_header = Map::<SamHeaderRecord>::builder()
                .insert(SORT_ORDER, sort_order::COORDINATE)
                .build()
                .map_err(|err| anyhow!("failed to build genomebam header: {err}"))?;
            builder = builder.set_header(sam_header);
        }
        let targets: Vec<(String, usize)> = if let Some(genome) = genome.as_ref() {
            genome.chromosomes.clone()
        } else {
            index
                .transcript_names
                .iter()
                .zip(index.transcript_lengths.iter())
                .map(|(name, length)| (name.clone(), usize::try_from(*length).unwrap_or(1).max(1)))
                .collect()
        };
        for (name, length) in targets {
            let length = NonZero::new(length.max(1))
                .ok_or_else(|| anyhow!("invalid transcript length for {name}"))?;
            builder = builder
                .add_reference_sequence(name.as_str(), Map::<ReferenceSequence>::new(length));
        }
        let header = builder.build();
        let path = out_dir.join("pseudoalignments.bam");
        let file = File::create(&path)?;
        let mut writer = noodles_bam::io::Writer::new(file);
        writer
            .write_header(&header)
            .map_err(|err| anyhow!("failed to write pseudobam header: {err}"))?;
        Ok(Self {
            writer,
            header,
            genome,
            pending: Vec::new(),
            path,
        })
    }

    fn alignment_for(&self, read: &PseudoBamRead<'_>) -> Option<GenomeAlignment> {
        let ec = usize::try_from(read.ec_id)
            .ok()
            .and_then(|idx| read.ec_list.get(idx))?;
        let transcript_id = if let Some(placement) = read.placement {
            placement.transcript_id
        } else if read.fallback_to_ec_transcript {
            *ec.first()?
        } else {
            return None;
        };
        let transcript_alignment = || {
            let reference_sequence_id = usize::try_from(transcript_id).ok()?;
            let start = read
                .placement
                .map(|placement| placement.start)
                .unwrap_or_default();
            let cigar = [Op::new(Kind::Match, read.sequence.len())]
                .into_iter()
                .collect();
            let reverse = read.placement.is_some_and(|placement| placement.reverse);
            Some(GenomeAlignment {
                reference_sequence_id,
                start,
                cigar,
                reverse,
            })
        };
        if let Some(genome) = self.genome.as_ref() {
            read.placement
                .and_then(|placement| genome.project(placement, read.sequence.len()))
        } else {
            transcript_alignment()
        }
    }

    fn write_read(&mut self, read: PseudoBamRead<'_>) -> Result<()> {
        if let Some(alignment) = self.alignment_for(&read) {
            self.write_aligned_read(read, alignment)
        } else {
            self.write_unmapped_read(read)
        }
    }

    fn write_unmapped_read(&mut self, read: PseudoBamRead<'_>) -> Result<()> {
        let barcode = std::str::from_utf8(read.barcode).unwrap_or("");
        let umi = std::str::from_utf8(read.umi).unwrap_or("");
        let mut flags = read.flags | Flags::UNMAPPED;
        if read.flags.is_segmented() && read.mate_reference_sequence_id.is_none() {
            flags |= Flags::MATE_UNMAPPED;
        }
        let mut builder = RecordBuf::builder()
            .set_name(read.name)
            .set_flags(flags)
            .set_mapping_quality(MappingQuality::MIN)
            .set_sequence(Sequence::from(read.sequence.to_vec()))
            .set_quality_scores(QualityScores::from(read.quality.to_vec()))
            .set_data(
                [
                    (Tag::CELL_BARCODE_SEQUENCE, SamBufDataValue::from(barcode)),
                    (Tag::UMI_SEQUENCE, SamBufDataValue::from(umi)),
                ]
                .into_iter()
                .collect(),
            );
        if let (Some(mate_reference_sequence_id), Some(mate_start)) =
            (read.mate_reference_sequence_id, read.mate_alignment_start)
        {
            let mate_alignment_start =
                Position::try_from(mate_start.saturating_add(1)).unwrap_or(Position::MIN);
            builder = builder
                .set_mate_reference_sequence_id(mate_reference_sequence_id)
                .set_mate_alignment_start(mate_alignment_start)
                .set_template_length(read.template_length);
        };
        let record = builder.build();
        if self.genome.is_some() {
            self.pending.push((
                BamSortKey {
                    reference_sequence_id: usize::MAX,
                    alignment_start: usize::MAX,
                    input_order: self.pending.len() as u64,
                },
                record,
            ));
        } else {
            self.writer
                .write_alignment_record(&self.header, &record)
                .map_err(|err| anyhow!("failed to write pseudobam record: {err}"))?;
        }
        Ok(())
    }

    fn write_aligned_read(
        &mut self,
        read: PseudoBamRead<'_>,
        alignment: GenomeAlignment,
    ) -> Result<()> {
        let Some(ec) = usize::try_from(read.ec_id)
            .ok()
            .and_then(|idx| read.ec_list.get(idx))
        else {
            return Ok(());
        };
        let alignment_start =
            Position::try_from(alignment.start.saturating_add(1)).unwrap_or(Position::MIN);
        let flags = read.flags
            | if read.flags.is_segmented() && read.mate_reference_sequence_id.is_none() {
                Flags::MATE_UNMAPPED
            } else {
                Flags::empty()
            }
            | if alignment.reverse {
                Flags::REVERSE_COMPLEMENTED
            } else {
                Flags::empty()
            }
            | if read.mate_reverse {
                Flags::MATE_REVERSE_COMPLEMENTED
            } else {
                Flags::empty()
            };
        let barcode = std::str::from_utf8(read.barcode).unwrap_or("");
        let umi = std::str::from_utf8(read.umi).unwrap_or("");
        let sequence = if alignment.reverse {
            kallistors::util::reverse_complement(read.sequence)
        } else {
            read.sequence.to_vec()
        };
        let quality = if alignment.reverse {
            read.quality.iter().rev().copied().collect::<Vec<_>>()
        } else {
            read.quality.to_vec()
        };
        let mut builder = RecordBuf::builder()
            .set_name(read.name)
            .set_flags(flags)
            .set_reference_sequence_id(alignment.reference_sequence_id)
            .set_alignment_start(alignment_start)
            .set_mapping_quality(MappingQuality::MIN)
            .set_cigar(alignment.cigar)
            .set_sequence(Sequence::from(sequence))
            .set_quality_scores(QualityScores::from(quality))
            .set_data(
                [
                    (Tag::CELL_BARCODE_SEQUENCE, SamBufDataValue::from(barcode)),
                    (Tag::UMI_SEQUENCE, SamBufDataValue::from(umi)),
                    (
                        Tag::ALIGNMENT_HIT_COUNT,
                        SamBufDataValue::from(ec.len() as i32),
                    ),
                ]
                .into_iter()
                .collect(),
            );
        if let (Some(mate_reference_sequence_id), Some(mate_start)) =
            (read.mate_reference_sequence_id, read.mate_alignment_start)
        {
            let mate_alignment_start =
                Position::try_from(mate_start.saturating_add(1)).unwrap_or(Position::MIN);
            builder = builder
                .set_mate_reference_sequence_id(mate_reference_sequence_id)
                .set_mate_alignment_start(mate_alignment_start)
                .set_template_length(read.template_length);
        }
        let record = builder.build();
        if self.genome.is_some() {
            self.pending.push((
                BamSortKey {
                    reference_sequence_id: alignment.reference_sequence_id,
                    alignment_start: alignment.start,
                    input_order: self.pending.len() as u64,
                },
                record,
            ));
        } else {
            self.writer
                .write_alignment_record(&self.header, &record)
                .map_err(|err| anyhow!("failed to write pseudobam record: {err}"))?;
        }
        Ok(())
    }

    pub(super) fn finish(mut self) -> Result<()> {
        self.pending.sort_by_key(|(key, _record)| *key);
        for (_key, record) in self.pending {
            self.writer
                .write_alignment_record(&self.header, &record)
                .map_err(|err| anyhow!("failed to write pseudobam record: {err}"))?;
        }
        self.writer
            .try_finish()
            .map_err(|err| anyhow!("failed to finish pseudobam: {err}"))?;
        if self.genome.is_some() {
            let index = noodles_bam::fs::index(&self.path)
                .map_err(|err| anyhow!("failed to index genomebam: {err}"))?;
            match index {
                noodles_bam::Index::Bai(index) => {
                    noodles_bam::bai::fs::write(self.path.with_extension("bam.bai"), &index)
                        .map_err(|err| anyhow!("failed to write genomebam BAI index: {err}"))?;
                }
                noodles_bam::Index::Csi(index) => {
                    noodles_csi::fs::write(self.path.with_extension("bam.csi"), &index)
                        .map_err(|err| anyhow!("failed to write genomebam CSI index: {err}"))?;
                }
            }
        }
        Ok(())
    }
}

pub(super) fn write_pseudobam_reads(
    writer: &mut PseudoBamWriter,
    read: PseudoBamBusRead<'_>,
) -> Result<()> {
    let read_name = pseudobam_read_name(read.batch, read.spec, read.read_id);
    if read.spec.paired {
        let sequences = read
            .spec
            .seq
            .iter()
            .take(2)
            .map(|seq_spec| {
                extract_sequence_slice(read.batch, *seq_spec, read.spec, read.tag_present)
            })
            .collect::<Result<Vec<_>>>()?;
        let qualities = read
            .spec
            .seq
            .iter()
            .take(2)
            .map(|seq_spec| {
                extract_quality_slice(read.batch, *seq_spec, read.spec, read.tag_present)
            })
            .collect::<Result<Vec<_>>>()?;
        let alignments = sequences
            .iter()
            .enumerate()
            .map(|(mate_idx, sequence)| {
                let pseudo_read = PseudoBamRead {
                    name: read_name.clone(),
                    flags: paired_segment_flags(mate_idx),
                    mate_reference_sequence_id: None,
                    mate_alignment_start: None,
                    mate_reverse: false,
                    template_length: 0,
                    ec_id: read.ec_id,
                    ec_list: read.ec_list,
                    placement: read.mate_placements.get(mate_idx).copied().flatten(),
                    barcode: read.barcode,
                    umi: read.umi,
                    sequence,
                    quality: qualities
                        .get(mate_idx)
                        .map_or(&[][..], std::vec::Vec::as_slice),
                    fallback_to_ec_transcript: read
                        .mate_has_ec
                        .get(mate_idx)
                        .copied()
                        .unwrap_or(false),
                };
                writer.alignment_for(&pseudo_read)
            })
            .collect::<Vec<_>>();
        let template_lengths = if alignments.len() == 2 {
            match (alignments[0].as_ref(), alignments[1].as_ref()) {
                (Some(left), Some(right)) => template_lengths(left, right),
                _ => (0, 0),
            }
        } else {
            (0, 0)
        };
        let both_mates_mapped = alignments.len() == 2 && alignments.iter().all(Option::is_some);
        for (mate_idx, sequence) in sequences.iter().enumerate() {
            let mate_idx_other = usize::from(mate_idx == 0);
            let mate_alignment = alignments.get(mate_idx_other).and_then(Option::as_ref);
            let mut flags = paired_segment_flags(mate_idx);
            if both_mates_mapped {
                flags |= Flags::PROPERLY_SEGMENTED;
            }
            let pseudo_read = PseudoBamRead {
                name: read_name.clone(),
                flags,
                mate_reference_sequence_id: mate_alignment
                    .map(|alignment| alignment.reference_sequence_id),
                mate_alignment_start: mate_alignment.map(|alignment| alignment.start),
                mate_reverse: mate_alignment.is_some_and(|alignment| alignment.reverse),
                template_length: if mate_idx == 0 {
                    template_lengths.0
                } else {
                    template_lengths.1
                },
                ec_id: read.ec_id,
                ec_list: read.ec_list,
                placement: read.mate_placements.get(mate_idx).copied().flatten(),
                barcode: read.barcode,
                umi: read.umi,
                sequence,
                quality: qualities
                    .get(mate_idx)
                    .map_or(&[][..], std::vec::Vec::as_slice),
                fallback_to_ec_transcript: read.mate_has_ec.get(mate_idx).copied().unwrap_or(false),
            };
            writer.write_read(pseudo_read)?;
        }
    } else {
        let sequence = concat_sequence_slices(read.batch, read.spec, read.tag_present)?;
        let quality = concat_quality_slices(read.batch, read.spec, read.tag_present)?;
        writer.write_read(PseudoBamRead {
            name: read_name,
            flags: Flags::empty(),
            mate_reference_sequence_id: None,
            mate_alignment_start: None,
            mate_reverse: false,
            template_length: 0,
            ec_id: read.ec_id,
            ec_list: read.ec_list,
            placement: read.placement,
            barcode: read.barcode,
            umi: read.umi,
            sequence: &sequence,
            quality: &quality,
            fallback_to_ec_transcript: true,
        })?;
    }
    Ok(())
}

fn paired_segment_flags(mate_idx: usize) -> Flags {
    Flags::SEGMENTED
        | if mate_idx == 0 {
            Flags::FIRST_SEGMENT
        } else {
            Flags::LAST_SEGMENT
        }
}

fn template_lengths(left: &GenomeAlignment, right: &GenomeAlignment) -> (i32, i32) {
    // TLEN spans reference coordinates, including intron skips in genome BAM.
    // Read lengths alone would undercount a spliced fragment's outer span.
    if left.reference_sequence_id != right.reference_sequence_id {
        return (0, 0);
    }
    let left_end = left
        .start
        .saturating_add(reference_span(&left.cigar).max(1));
    let right_end = right
        .start
        .saturating_add(reference_span(&right.cigar).max(1));
    let outer_start = left.start.min(right.start);
    let outer_end = left_end.max(right_end);
    let len = i32::try_from(outer_end.saturating_sub(outer_start)).unwrap_or(i32::MAX);
    if left.start <= right.start {
        (len, -len)
    } else {
        (-len, len)
    }
}

fn reference_span(cigar: &Cigar) -> usize {
    cigar
        .iter()
        .filter_map(|result| result.ok())
        .filter(|op| op.kind().consumes_reference())
        .map(|op| op.len())
        .sum()
}

fn pseudobam_read_name(
    batch: &[kallistors::io::FastqRecord],
    spec: &TechnologySpec,
    read_id: u64,
) -> Vec<u8> {
    let record = spec
        .seq
        .first()
        .and_then(|seq_spec| seq_spec.file)
        .and_then(|file| batch.get(file))
        .or_else(|| batch.first());
    let Some(record) = record else {
        return format!("read{read_id}").into_bytes();
    };
    let mut name = record.header.strip_prefix(b"@").unwrap_or(&record.header);
    name = name.trim_ascii_end();
    if let Some(end) = name.iter().position(|base| matches!(base, b' ' | b'\t')) {
        name = &name[..end];
    }
    if name.is_empty() {
        format!("read{read_id}").into_bytes()
    } else {
        name.to_vec()
    }
}
