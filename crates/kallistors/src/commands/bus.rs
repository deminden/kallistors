use anyhow::{Result, anyhow, bail};
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
use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Seek, SeekFrom, Write};
use std::num::NonZero;
use std::path::{Path, PathBuf};
use std::time::Instant;

use kallistors::io::ReadSource;
use noodles_sam::alignment::record::data::field::Value as SamDataValue;

mod output;

use output::{
    BusRecord, RunInfo, fake_barcode, write_batch_long_flens, write_batch_paired_flens,
    write_bus_header, write_bus_record, write_cells, write_long_flens, write_matrix_ec,
    write_novel_read, write_paired_flens, write_run_info, write_sample_barcodes, write_transcripts,
    write_unmapped_ratios,
};

#[derive(Default)]
pub struct BusArgs {
    pub index: Option<PathBuf>,
    pub out_dir: Option<PathBuf>,
    pub technology: Option<String>,
    pub reads: Vec<PathBuf>,
    pub threads: usize,
    pub list: bool,
    pub max_reads: Option<u64>,
    pub batch: Option<PathBuf>,
    pub bam: bool,
    pub pseudobam: bool,
    pub genomebam: bool,
    pub gtf: Option<PathBuf>,
    pub chromosomes: Option<PathBuf>,
    pub batch_barcodes: bool,
    pub interleaved: bool,
    pub tag_sequence: Option<String>,
    pub long: bool,
    pub platform: Option<String>,
    pub error_rate: Option<f64>,
    pub threshold: Option<f64>,
    pub unmapped: bool,
    pub aa: bool,
    pub num: bool,
    pub paired: bool,
    pub unstranded: bool,
    pub fr_stranded: bool,
    pub rf_stranded: bool,
    pub do_union: bool,
    pub no_jump: bool,
    pub verbose: bool,
    pub dfk_onlist: bool,
    pub kallisto_enum: bool,
    pub kallisto_strict: bool,
    pub kallisto_local_fallback: bool,
    pub kallisto_fallback: bool,
    pub discard_special_only: bool,
    pub skip_overcrowded_minimizer: bool,
    pub kallisto_direct_kmer: bool,
    pub kallisto_bifrost_find: bool,
    pub kallisto_sparse_hits: bool,
}

#[derive(Debug, Clone, Copy)]
struct SliceSpec {
    file: Option<usize>,
    start: usize,
    stop: Option<usize>,
}

#[derive(Debug, Clone)]
struct TechnologySpec {
    nfiles: usize,
    seq: Vec<SliceSpec>,
    bc: Vec<SliceSpec>,
    umi: Vec<SliceSpec>,
    paired: bool,
    keep_fastq_comments: bool,
    default_strand: Option<kallistors::pseudoalign::StrandSpecific>,
    tag_len: Option<usize>,
}

impl TechnologySpec {
    fn bulk_like(&self) -> bool {
        !self.keep_fastq_comments && self.umi.first().is_none_or(|spec| spec.file.is_none())
    }

    fn with_default_strand(mut self, strand: kallistors::pseudoalign::StrandSpecific) -> Self {
        self.default_strand = Some(strand);
        self
    }
}

#[derive(Debug, Clone)]
struct TagConfig {
    sequence: Vec<u8>,
    binary: u64,
}

struct ConfiguredTag {
    config: TagConfig,
    inferred_default: bool,
}

struct UmiValue {
    binary: u64,
    flag: u32,
    len: usize,
    seq: Vec<u8>,
    tag_present: bool,
    disable_strand_specific: bool,
}

struct ReadGroup {
    id: Option<String>,
    files: Vec<PathBuf>,
    batch_index: usize,
    interleaved: bool,
}

struct BusWriteState<'a, W: Write> {
    writer: &'a mut W,
    ec_map: &'a mut HashMap<Vec<u32>, i32>,
    ec_list: &'a mut Vec<Vec<u32>>,
    bc_len_hist: &'a mut [u64; 33],
    umi_len_hist: &'a mut [u64; 33],
    long_len_sums: &'a mut [u64],
    long_len_counts: &'a mut [u64],
}

#[derive(Debug, Clone, Copy)]
struct BusProcessingConfig {
    options: kallistors::pseudoalign::PseudoalignOptions,
    threshold: f64,
    max_reads: Option<u64>,
}

struct BusAlignment {
    ec: Option<Vec<u32>>,
    fragment_length: Option<i64>,
    placement: Option<kallistors::pseudoalign::PseudoalignPlacement>,
    mate_placements: Vec<Option<kallistors::pseudoalign::PseudoalignPlacement>>,
    mate_has_ec: Vec<bool>,
    frame_clashes: u64,
}

#[derive(Clone)]
struct GenomeModel {
    chromosomes: Vec<(String, usize)>,
    transcripts: Vec<Option<TranscriptProjection>>,
    genes: Vec<GeneInfo>,
}

#[derive(Clone)]
struct GeneInfo {
    id: String,
    name: String,
}

#[derive(Clone)]
struct TranscriptProjection {
    chromosome_id: usize,
    negative_strand: bool,
    transcript_len: usize,
    exons: Vec<Exon>,
}

#[derive(Clone, Copy)]
struct Exon {
    start: usize,
    end: usize,
}

#[derive(Clone)]
struct GenomeAlignment {
    reference_sequence_id: usize,
    start: usize,
    cigar: Cigar,
    reverse: bool,
}

struct GenomeSegment {
    start: usize,
    end: usize,
}

struct PseudoBamWriter {
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

struct PseudoBamBusRead<'a> {
    read_id: u64,
    ec_id: i32,
    ec_list: &'a [Vec<u32>],
    batch: &'a [kallistors::io::FastqRecord],
    spec: &'a TechnologySpec,
    tag_present: bool,
    barcode: &'a [u8],
    umi: &'a [u8],
    placement: Option<kallistors::pseudoalign::PseudoalignPlacement>,
    mate_placements: &'a [Option<kallistors::pseudoalign::PseudoalignPlacement>],
    mate_has_ec: &'a [bool],
}

impl GenomeModel {
    fn from_gtf(
        gtf: &Path,
        chromosomes: Option<&Path>,
        index: &kallistors::pseudoalign::BifrostIndex,
    ) -> Result<Self> {
        let mut chromosome_ids = HashMap::<String, usize>::new();
        let mut chromosome_lengths = Vec::<usize>::new();
        let mut chromosome_names = Vec::<String>::new();
        if let Some(path) = chromosomes {
            for line in BufReader::new(File::open(path)?).lines() {
                let line = line?;
                let line = line.trim();
                if line.is_empty() || line.starts_with('#') {
                    continue;
                }
                let mut fields = line.split_whitespace();
                let Some(name) = fields.next() else {
                    continue;
                };
                let Some(length) = fields.next() else {
                    bail!("chromosome file line is missing a length: {line}");
                };
                let length = length
                    .parse::<usize>()
                    .map_err(|_| anyhow!("invalid chromosome length in line: {line}"))?;
                if length == 0 {
                    bail!("invalid chromosome length in line: {line}");
                }
                if chromosome_ids.contains_key(name) {
                    bail!("duplicate chromosome name in line: {line}");
                }
                let id = chromosome_names.len();
                chromosome_ids.insert(name.to_string(), id);
                chromosome_names.push(name.to_string());
                chromosome_lengths.push(length);
            }
        }

        let transcript_ids = index
            .transcript_names
            .iter()
            .enumerate()
            .map(|(idx, name)| (name.as_str(), idx))
            .collect::<HashMap<_, _>>();
        let mut raw = std::iter::repeat_with(|| None::<RawTranscriptProjection>)
            .take(index.transcript_names.len())
            .collect::<Vec<_>>();
        let mut gene_ids = HashMap::<String, usize>::new();
        let mut genes = Vec::<GeneInfo>::new();
        for line in BufReader::new(File::open(gtf)?).lines() {
            let line = line?;
            if line.starts_with('#') || line.trim().is_empty() {
                continue;
            }
            let fields = line.split('\t').collect::<Vec<_>>();
            if fields.len() < 9 {
                continue;
            }
            if fields[2] == "gene" {
                if let Some(gene_id) = versioned_gtf_attribute(fields[8], "gene_id", "gene_version")
                    && !gene_ids.contains_key(&gene_id)
                {
                    let gene_name = gtf_attribute(fields[8], "gene_name").unwrap_or_default();
                    gene_ids.insert(gene_id.clone(), genes.len());
                    genes.push(GeneInfo {
                        id: gene_id,
                        name: gene_name,
                    });
                }
                continue;
            }
            if fields[2] != "exon" {
                continue;
            }
            let Some(transcript_name) =
                versioned_gtf_attribute(fields[8], "transcript_id", "transcript_version")
            else {
                continue;
            };
            let transcript_id =
                if let Some(&transcript_id) = transcript_ids.get(transcript_name.as_str()) {
                    transcript_id
                } else if let Some(unversioned) = gtf_attribute(fields[8], "transcript_id") {
                    let Some(&transcript_id) = transcript_ids.get(unversioned.as_str()) else {
                        continue;
                    };
                    transcript_id
                } else {
                    continue;
                };
            if let Some(gene_id) = versioned_gtf_attribute(fields[8], "gene_id", "gene_version")
                && !gene_ids.contains_key(&gene_id)
            {
                let gene_name = gtf_attribute(fields[8], "gene_name").unwrap_or_default();
                gene_ids.insert(gene_id.clone(), genes.len());
                genes.push(GeneInfo {
                    id: gene_id,
                    name: gene_name,
                });
            }
            let start_one_based = fields[3]
                .parse::<usize>()
                .map_err(|_| anyhow!("invalid GTF exon start: {}", fields[3]))?;
            if start_one_based == 0 {
                bail!("invalid GTF exon start: {}", fields[3]);
            }
            let start = start_one_based - 1;
            let end = fields[4]
                .parse::<usize>()
                .map_err(|_| anyhow!("invalid GTF exon end: {}", fields[4]))?;
            if end <= start {
                bail!("invalid GTF exon interval: {}-{}", fields[3], fields[4]);
            }
            let chromosome_id = if let Some(id) = chromosome_ids.get(fields[0]).copied() {
                id
            } else {
                let id = chromosome_names.len();
                chromosome_ids.insert(fields[0].to_string(), id);
                chromosome_names.push(fields[0].to_string());
                chromosome_lengths.push(0);
                id
            };
            chromosome_lengths[chromosome_id] = chromosome_lengths[chromosome_id].max(end);
            let negative_strand = match fields[6] {
                "+" => false,
                "-" => true,
                strand => bail!("invalid GTF exon strand: {strand}"),
            };
            let entry = raw[transcript_id].get_or_insert_with(|| RawTranscriptProjection {
                chromosome_id,
                negative_strand,
                exons: Vec::new(),
            });
            if entry.chromosome_id != chromosome_id {
                bail!("GTF transcript {transcript_name} has exons on multiple chromosomes");
            }
            if entry.negative_strand != negative_strand {
                bail!("GTF transcript {transcript_name} has exons on multiple strands");
            }
            entry.exons.push(Exon { start, end });
        }

        let transcripts = raw
            .into_iter()
            .map(|entry| {
                entry.map(|mut entry| {
                    if entry.negative_strand {
                        entry
                            .exons
                            .sort_by_key(|exon| std::cmp::Reverse(exon.start));
                    } else {
                        entry.exons.sort_by_key(|exon| exon.start);
                    }
                    TranscriptProjection {
                        chromosome_id: entry.chromosome_id,
                        negative_strand: entry.negative_strand,
                        transcript_len: entry
                            .exons
                            .iter()
                            .map(|exon| exon.end.saturating_sub(exon.start))
                            .sum(),
                        exons: entry.exons,
                    }
                })
            })
            .collect();

        Ok(Self {
            chromosomes: chromosome_names
                .into_iter()
                .zip(chromosome_lengths)
                .collect(),
            transcripts,
            genes,
        })
    }

    fn project(
        &self,
        placement: kallistors::pseudoalign::PseudoalignPlacement,
        read_len: usize,
    ) -> Option<GenomeAlignment> {
        let transcript = self
            .transcripts
            .get(placement.transcript_id as usize)?
            .as_ref()?;
        let start = if placement.start.saturating_add(read_len) > transcript.transcript_len {
            transcript.transcript_len.saturating_sub(read_len)
        } else {
            placement.start
        };
        transcript.project(start, read_len, placement.reverse)
    }
}

struct RawTranscriptProjection {
    chromosome_id: usize,
    negative_strand: bool,
    exons: Vec<Exon>,
}

impl TranscriptProjection {
    fn project(
        &self,
        transcript_start: usize,
        read_len: usize,
        read_reverse: bool,
    ) -> Option<GenomeAlignment> {
        let mut offset = transcript_start;
        let mut remaining = read_len;
        let mut segments = Vec::new();

        for exon in &self.exons {
            let exon_len = exon.end.saturating_sub(exon.start);
            if offset >= exon_len {
                offset -= exon_len;
                continue;
            }
            let take = remaining.min(exon_len - offset);
            let segment_start = if self.negative_strand {
                exon.end.saturating_sub(offset + take)
            } else {
                exon.start + offset
            };
            let segment_end = segment_start + take;
            segments.push(GenomeSegment {
                start: segment_start,
                end: segment_end,
            });
            remaining -= take;
            offset = 0;
            if remaining == 0 {
                break;
            }
        }
        if remaining > 0 {
            return None;
        }
        segments.sort_by_key(|segment| segment.start);
        let mut ops = Vec::new();
        let mut previous_segment_end = None;
        for segment in &segments {
            if let Some(previous_end) = previous_segment_end {
                let skip = segment.start.saturating_sub(previous_end);
                if skip > 0 {
                    ops.push(Op::new(Kind::Skip, skip));
                }
            }
            ops.push(Op::new(
                Kind::Match,
                segment.end.saturating_sub(segment.start),
            ));
            previous_segment_end = Some(segment.end);
        }
        Some(GenomeAlignment {
            reference_sequence_id: self.chromosome_id,
            start: segments.first()?.start,
            cigar: ops.into_iter().collect(),
            reverse: read_reverse ^ self.negative_strand,
        })
    }
}

fn gtf_attribute(attributes: &str, key: &str) -> Option<String> {
    for attribute in attributes.split(';') {
        let attribute = attribute.trim();
        let mut parts = attribute.splitn(2, char::is_whitespace);
        let Some(name) = parts.next() else {
            continue;
        };
        if name != key {
            continue;
        }
        let value = parts.next()?.trim().trim_matches('"');
        if !value.is_empty() {
            return Some(value.to_string());
        }
    }
    None
}

fn versioned_gtf_attribute(attributes: &str, key: &str, version_key: &str) -> Option<String> {
    let mut value = gtf_attribute(attributes, key)?;
    if !value.contains('.')
        && let Some(version) = gtf_attribute(attributes, version_key)
    {
        value.push('.');
        value.push_str(&version);
    }
    Some(value)
}

impl PseudoBamWriter {
    fn create(
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

    fn finish(mut self) -> Result<()> {
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
            noodles_bam::bai::fs::write(self.path.with_extension("bam.bai"), &index)
                .map_err(|err| anyhow!("failed to write genomebam index: {err}"))?;
        }
        Ok(())
    }
}

struct ThresholdChoice {
    value: f64,
    report: Option<ThresholdReport>,
}

enum ThresholdReport {
    Computed(f64),
    InvalidSupplied,
    InvalidComputed(f64),
}

fn long_read_threshold(
    threshold: Option<f64>,
    error_rate: Option<f64>,
    k: usize,
) -> ThresholdChoice {
    if let Some(threshold) = threshold {
        if 0.0 < threshold && threshold < 1.0 {
            return ThresholdChoice {
                value: threshold,
                report: None,
            };
        }
        return ThresholdChoice {
            value: 0.8,
            report: Some(ThresholdReport::InvalidSupplied),
        };
    }
    if let Some(error_rate) = error_rate {
        let computed = (1.0 / error_rate - 2.0 * k as f64) * error_rate;
        if 0.0 < computed && computed < 1.0 {
            return ThresholdChoice {
                value: computed,
                report: Some(ThresholdReport::Computed(computed)),
            };
        }
        return ThresholdChoice {
            value: 0.8,
            report: Some(ThresholdReport::InvalidComputed(computed)),
        };
    }
    ThresholdChoice {
        value: 0.8,
        report: None,
    }
}

fn report_threshold_choice(report: &ThresholdReport) {
    match report {
        ThresholdReport::Computed(threshold) => {
            eprintln!("Using computed threshold {threshold}");
        }
        ThresholdReport::InvalidSupplied => {
            eprintln!(
                "Threshold not in (0,1). Setting default threshold for unmapped kmers to 0.8"
            );
        }
        ThresholdReport::InvalidComputed(threshold) => {
            eprintln!(
                "Supplied and computed threshold are invalid, using default value of 0.8 (computed: {threshold})"
            );
        }
    }
}

fn bus_pseudoalign_options(
    args: &BusArgs,
    strand_specific: Option<kallistors::pseudoalign::StrandSpecific>,
) -> kallistors::pseudoalign::PseudoalignOptions {
    kallistors::pseudoalign::PseudoalignOptions {
        kallisto_enum: args.kallisto_enum,
        kallisto_strict: args.kallisto_strict,
        kallisto_local_fallback: args.kallisto_local_fallback,
        kallisto_fallback: args.kallisto_fallback,
        discard_special_only: args.discard_special_only,
        skip_overcrowded_minimizer: args.skip_overcrowded_minimizer,
        kallisto_direct_kmer: args.kallisto_direct_kmer,
        kallisto_bifrost_find: args.kallisto_bifrost_find,
        kallisto_sparse_hits: args.kallisto_sparse_hits,
        strand_specific,
        do_union: args.do_union,
        no_jump: args.no_jump,
        dfk_onlist: args.dfk_onlist || args.aa,
        investigation: super::investigation_options_from_env(),
        ..kallistors::pseudoalign::PseudoalignOptions::default()
    }
}

fn bus_strand_specific(
    args: &BusArgs,
    spec: &TechnologySpec,
    index: &kallistors::pseudoalign::BifrostIndex,
) -> Option<kallistors::pseudoalign::StrandSpecific> {
    if args.unstranded {
        None
    } else if args.fr_stranded {
        Some(kallistors::pseudoalign::StrandSpecific::FR)
    } else if args.rf_stranded {
        Some(kallistors::pseudoalign::StrandSpecific::RF)
    } else if index.has_strand_annotations() {
        spec.default_strand
    } else {
        None
    }
}

fn user_specified_strand(args: &BusArgs) -> bool {
    args.unstranded || args.fr_stranded || args.rf_stranded
}

fn report_default_strand(strand: kallistors::pseudoalign::StrandSpecific) {
    let option = match strand {
        kallistors::pseudoalign::StrandSpecific::FR => "--fr-stranded",
        kallistors::pseudoalign::StrandSpecific::RF => "--rf-stranded",
    };
    eprintln!(
        "[bus] Note: Strand option was not specified; setting it to {option} for specified technology"
    );
}

fn report_missing_strand_annotations() {
    eprintln!(
        "[bus] Note: Strand option was not specified; index has no strand annotations, processing as --unstranded"
    );
}

pub fn run(args: BusArgs) -> Result<()> {
    if args.list {
        print_technology_list();
        return Ok(());
    }
    let technology = args.technology.as_deref().unwrap_or("Bulk");
    let index_path = args
        .index
        .as_deref()
        .ok_or_else(|| anyhow!("kallisto index file missing"))?;
    if !index_path.exists() {
        bail!("kallisto index file not found {}", index_path.display());
    }
    let out_dir = args
        .out_dir
        .as_deref()
        .ok_or_else(|| anyhow!("need to specify output directory"))?;
    if out_dir.exists() && !out_dir.is_dir() {
        bail!("file {} exists and is not a directory", out_dir.display());
    }
    if args.batch.is_none() && args.reads.is_empty() {
        if args.technology.is_none() {
            bail!(
                "the technology must be specified via -x, use \"bulk\" for regular RNA-seq reads"
            );
        }
        bail!("missing read files");
    }
    if args.threads == 0 {
        bail!("invalid number of threads 0");
    }
    let write_pseudobam = args.pseudobam || args.genomebam;
    if args.pseudobam && args.genomebam {
        bail!("--pseudobam and --genomebam are mutually exclusive");
    }
    if args.genomebam && args.gtf.is_none() {
        bail!("--genomebam requires --gtf");
    }
    if args.gtf.is_some() && !args.genomebam {
        bail!("--gtf is only valid with --genomebam");
    }
    if args.chromosomes.is_some() && !args.genomebam {
        bail!("--chromosomes is only valid with --genomebam");
    }
    if let Some(gtf) = args.gtf.as_deref()
        && !gtf.exists()
    {
        bail!("GTF file not found: {}", gtf.display());
    }
    if let Some(chromosomes) = args.chromosomes.as_deref()
        && !chromosomes.exists()
    {
        bail!("Chromosome file not found: {}", chromosomes.display());
    }
    if args.batch.is_some() && !args.reads.is_empty() {
        bail!("cannot specify batch mode and supply read files");
    }
    if let Some(batch) = args.batch.as_deref()
        && !batch.exists()
    {
        bail!("file not found {}", batch.display());
    }
    if args.batch.is_some() {
        eprintln!("[bus] will try running read files supplied in batch file");
        if args.technology.is_none() && args.paired {
            eprintln!(
                "[bus] --paired ignored; single/paired-end is inferred from number of files supplied"
            );
        }
    }
    if args.bam && args.batch.is_some() {
        bail!("--bam cannot be combined with batch mode");
    }
    if args.bam && args.interleaved {
        bail!("--bam cannot be combined with interleaved FASTQ input");
    }
    if args.bam && args.reads.len() != 1 {
        bail!("--bam expects exactly one BAM file");
    }
    if args.bam && args.num {
        eprintln!("Warning: --bam option was used, so --num option will be ignored");
    }
    if args.bam && args.batch_barcodes {
        bail!("--batch-barcodes requires batch mode and cannot be used with --bam");
    }
    if args.bam && args.paired && !args.long && !args.aa {
        bail!("Paired reads are not compatible with the specified technology");
    }
    if args.batch_barcodes && args.batch.is_none() {
        bail!("--batch-barcodes requires batch mode");
    }
    if args.interleaved && args.batch.is_some() {
        bail!("interleaved input cannot be specified with a batch file");
    }
    if args.interleaved && args.reads.len() != 1 {
        bail!("interleaved input expects exactly one FASTQ file");
    }
    if args.fr_stranded && args.rf_stranded {
        bail!("--fr-stranded and --rf-stranded are mutually exclusive");
    }
    if args.unstranded && (args.fr_stranded || args.rf_stranded) {
        bail!("--unstranded cannot be combined with --fr-stranded or --rf-stranded");
    }
    if args.long
        && let Some(error_rate) = args.error_rate
        && error_rate <= 0.0
    {
        bail!("--error-rate must be greater than 0");
    }
    if args.long
        && let Some(platform) = args.platform.as_deref()
        && !platform.eq_ignore_ascii_case("PACBIO")
        && !platform.eq_ignore_ascii_case("ONT")
    {
        bail!("--platform must be PACBIO or ONT");
    }
    if (args.do_union || args.no_jump) && (args.long || args.aa) {
        bail!("--union and --no-jump are not compatible with --long or --aa");
    }
    let mut spec = technology_spec(technology)?;
    let technology_base = technology_base_upper(technology);
    if technology_base == "BULK" {
        if write_pseudobam {
            bail!("Pseudobam not supported yet in this mode");
        }
        if args.bam {
            bail!("--bam not supported in this mode");
        }
        if args.tag_sequence.is_some() {
            bail!("--tag not supported in this mode");
        }
    }
    let tag_config = configure_tag(technology, args.tag_sequence.as_deref(), &mut spec)?;
    if let Some(tag) = tag_config.as_ref()
        && tag.inferred_default
    {
        eprintln!(
            "[bus] Using {} as UMI tag sequence",
            String::from_utf8_lossy(&tag.config.sequence)
        );
    }
    if args.bam && tag_config.is_some() {
        bail!("--tag is not supported with --bam input");
    }
    if write_pseudobam && args.bam {
        bail!("--pseudobam/--genomebam is only supported for FASTQ BUS input");
    }
    if args.technology.is_none()
        && !args.aa
        && !args.long
        && !spec.paired
        && should_infer_bulk_paired(&args)?
    {
        mark_technology_paired(&mut spec, "Bulk input")?;
    }
    if args.paired && !args.aa && !args.long && !spec.paired {
        mark_technology_paired(&mut spec, "--paired")?;
    } else if args.paired && args.long && !args.aa && technology_base == "SMARTSEQ2" {
        add_unpaired_mate_slice(&mut spec, "--paired")?;
    }
    if args.aa && args.paired && !spec.paired {
        eprintln!("[bus] --paired ignored; --aa only supports single-end reads");
    }
    if args.aa && spec.paired {
        bail!("--aa BUS mode currently supports single cDNA read technologies");
    }
    if write_pseudobam && !pseudobam_compatible(&spec) {
        bail!(
            "BAM output is currently only supported for technologies with a single cDNA read file or paired-end reads"
        );
    }
    let read_groups = read_groups(&args, spec.nfiles)?;
    if read_groups.is_empty() {
        bail!("missing read files");
    }
    validate_read_files(&args, &read_groups)?;
    if args.verbose {
        report_verbose_read_groups(&read_groups);
    }
    let synthetic_bc = spec.bc.first().is_some_and(|v| v.file.is_none());
    let fixed_bc_len = total_len(&spec.bc);
    if let Some(fixed_bc_len) = fixed_bc_len
        && fixed_bc_len > 32
    {
        bail!("barcode length {fixed_bc_len} exceeds BUS limit of 32 bases");
    }
    let fixed_umi_len = total_len(&spec.umi);
    if let Some(fixed_umi_len) = fixed_umi_len
        && fixed_umi_len > 32
    {
        bail!("UMI length {fixed_umi_len} exceeds BUS limit of 32 bases");
    }
    let batch_barcode_prefix_len = if args.batch_barcodes {
        let Some(fixed_bc_len) = fixed_bc_len else {
            bail!("--batch-barcodes requires a bounded barcode length");
        };
        if fixed_bc_len >= 32 {
            bail!("--batch-barcodes requires barcode length shorter than 32 bases");
        }
        Some(32 - fixed_bc_len)
    } else {
        None
    };
    if let Some(prefix_len) = batch_barcode_prefix_len
        && prefix_len == 0
    {
        bail!("--batch-barcodes requires barcode length shorter than 32 bases");
    }
    let bc_len = if args.batch_barcodes && fixed_bc_len.is_some_and(|v| v > 0) {
        32
    } else if synthetic_bc {
        16
    } else {
        fixed_bc_len.unwrap_or(0)
    };
    if !args.bam
        && args.batch.is_none()
        && !args.interleaved
        && !args.reads.len().is_multiple_of(spec.nfiles)
    {
        bail!(
            "technology {} expects files in groups of {}, got {}",
            technology,
            spec.nfiles,
            args.reads.len()
        );
    }

    let start = Instant::now();
    let index = if args.kallisto_direct_kmer {
        kallistors::pseudoalign::build_bifrost_index_with_kmer_pos(index_path, false)
    } else if args.kallisto_fallback {
        kallistors::pseudoalign::build_bifrost_index_with_kmer(index_path)
    } else {
        kallistors::pseudoalign::build_bifrost_index_with_positions_threaded(
            index_path,
            false,
            args.threads.max(1),
        )
    }
    .map_err(|err| anyhow!("pseudoalign failed: {err}"))?;
    let threshold = if args.long {
        let choice = long_read_threshold(args.threshold, args.error_rate, index.k);
        if let Some(report) = choice.report.as_ref() {
            report_threshold_choice(report);
        }
        choice.value
    } else {
        0.8
    };
    let strand_specific = bus_strand_specific(&args, &spec, &index);
    if strand_specific.is_some() && !strand_specific_compatible(&spec) {
        bail!(
            "Strand-specific read processing is only supported for technologies with a single cDNA read file or paired-end reads"
        );
    }
    if !user_specified_strand(&args)
        && let Some(strand) = strand_specific
        && spec.default_strand == Some(strand)
    {
        report_default_strand(strand);
    } else if !user_specified_strand(&args)
        && strand_specific.is_none()
        && spec.default_strand.is_some()
        && !index.has_strand_annotations()
    {
        report_missing_strand_annotations();
    }
    let options = bus_pseudoalign_options(&args, strand_specific);
    let processing = BusProcessingConfig {
        options,
        threshold,
        max_reads: args.max_reads.filter(|&value| value != 0),
    };

    fs::create_dir_all(out_dir)?;
    let mut bus_writer = BufWriter::new(File::create(out_dir.join("output.bus"))?);
    let genome_model = if args.genomebam {
        Some(GenomeModel::from_gtf(
            args.gtf.as_deref().expect("validated"),
            args.chromosomes.as_deref(),
            &index,
        )?)
    } else {
        None
    };
    if let Some(model) = genome_model.as_ref() {
        write_gene_list(&out_dir.join("matrix.genelist.txt"), model)?;
    }
    let mut pseudo_bam = if write_pseudobam {
        Some(PseudoBamWriter::create(out_dir, &index, genome_model)?)
    } else {
        None
    };
    let umi_len = if args.bam {
        0
    } else if spec.bulk_like() {
        1
    } else {
        fixed_umi_len.unwrap_or(0)
    };
    let header_bc_len = if args.bam { 0 } else { bc_len };
    write_bus_header(&mut bus_writer, header_bc_len, umi_len)?;

    let mut ec_map: HashMap<Vec<u32>, i32> = HashMap::new();
    let mut ec_list: Vec<Vec<u32>> = Vec::new();
    let mut records = 0u64;
    let mut aligned = 0u64;
    let mut unique = 0u64;
    let mut bc_len_hist = [0u64; 33];
    let mut umi_len_hist = [0u64; 33];
    let mut unmapped_ratios = Vec::new();
    let mut long_len_sums = vec![0u64; index.transcript_names.len()];
    let mut long_len_counts = vec![0u64; index.transcript_names.len()];
    let mut paired_flens = vec![0u32; kallistors::pseudoalign::MAX_FRAG_LEN as usize];
    let batch_count = read_groups
        .iter()
        .map(|group| group.batch_index)
        .max()
        .map_or(0, |idx| idx + 1);
    let mut batch_long_len_sums = vec![vec![0u64; index.transcript_names.len()]; batch_count];
    let mut batch_long_len_counts = vec![vec![0u64; index.transcript_names.len()]; batch_count];
    let mut batch_paired_flens =
        vec![vec![0u32; kallistors::pseudoalign::MAX_FRAG_LEN as usize]; batch_count];
    let mut frame_clashes = 0u64;
    let mut novel_writer = if args.long && args.unmapped {
        Some(BufWriter::new(File::create(out_dir.join("novel.fastq"))?))
    } else {
        None
    };

    if args.bam {
        let mut state = BusWriteState {
            writer: &mut bus_writer,
            ec_map: &mut ec_map,
            ec_list: &mut ec_list,
            bc_len_hist: &mut bc_len_hist,
            umi_len_hist: &mut umi_len_hist,
            long_len_sums: &mut long_len_sums,
            long_len_counts: &mut long_len_counts,
        };
        let bam_records = process_bam_records(
            &args,
            &index,
            processing,
            &mut state,
            &mut unmapped_ratios,
            &mut novel_writer,
        )?;
        records = bam_records.0;
        aligned = bam_records.1;
        unique = bam_records.2;
        frame_clashes = bam_records.3;
    } else {
        'groups: for group in &read_groups {
            let mut readers = group
                .files
                .iter()
                .map(|path| kallistors::io::open_fastq_reader(path))
                .collect::<kallistors::Result<Vec<_>>>()
                .map_err(|err| anyhow!("failed to open FASTQ: {err}"))?;

            loop {
                let batch = if group.interleaved {
                    next_interleaved_batch(&mut readers, spec.nfiles)?
                } else {
                    next_parallel_batch(&mut readers)?
                };
                if batch.is_empty() {
                    break;
                }
                if let Some(max_reads) = processing.max_reads
                    && records >= max_reads
                {
                    break 'groups;
                }
                records += 1;

                let Some(barcode_seq) =
                    barcode_sequence(&batch, &spec, group.batch_index, args.batch_barcodes)?
                else {
                    record_unmapped_ratio_for_skipped_read(
                        args.unmapped,
                        &mut unmapped_ratios,
                        &index,
                        &batch,
                        &spec,
                        options,
                    )?;
                    continue;
                };
                let Some(umi_value) =
                    umi_value(&batch, &spec, tag_config.as_ref().map(|tag| &tag.config))?
                else {
                    record_unmapped_ratio_for_skipped_read(
                        args.unmapped,
                        &mut unmapped_ratios,
                        &index,
                        &batch,
                        &spec,
                        options,
                    )?;
                    continue;
                };
                let umi_seq_len = umi_value.len;
                record_length(&mut bc_len_hist, barcode_seq.len());
                record_length(&mut umi_len_hist, umi_seq_len);
                let mut barcode_flag = 0;
                let barcode = string_to_binary(&barcode_seq, &mut barcode_flag, "barcode")?;

                let query = if args.long || args.unmapped {
                    Some(bus_query_sequence(&batch, &spec, umi_value.tag_present)?)
                } else {
                    None
                };
                let query_quality = if args.long || args.unmapped {
                    Some(bus_query_quality(&batch, &spec, umi_value.tag_present)?)
                } else {
                    None
                };
                let too_many_empty_kmers = if let Some(query) = query.as_deref() {
                    let ratio = kallistors::pseudoalign::unmapped_kmer_ratio_bifrost(
                        &index, query, options,
                    );
                    if args.unmapped {
                        unmapped_ratios.push(ratio);
                    }
                    args.long && ratio > threshold
                } else {
                    false
                };
                let mut read_options = options;
                if umi_value.disable_strand_specific {
                    read_options.strand_specific = None;
                }
                let alignment = pseudoalign_bus_read(
                    &index,
                    &batch,
                    &spec,
                    read_options,
                    umi_value.tag_present,
                    args.aa,
                )?;
                let disjoint_intersect = alignment.ec.is_none();
                let novel_long_read = args.long && (too_many_empty_kmers || disjoint_intersect);
                if args.long
                    && args.unmapped
                    && disjoint_intersect
                    && let (Some(writer), Some(query), Some(quality)) = (
                        novel_writer.as_mut(),
                        query.as_deref(),
                        query_quality.as_deref(),
                    )
                {
                    write_novel_read(writer, "unmapped", query, quality)?;
                }
                if novel_long_read
                    && let (Some(writer), Some(query), Some(quality)) = (
                        novel_writer.as_mut(),
                        query.as_deref(),
                        query_quality.as_deref(),
                    )
                {
                    let label = if disjoint_intersect {
                        "novel_disjointIntersect"
                    } else {
                        "novel_tooManyEmptyKmers"
                    };
                    write_novel_read(writer, label, query, quality)?;
                }
                if args.long
                    && !novel_long_read
                    && let Some(ec) = alignment.ec.as_deref()
                {
                    record_long_lengths(
                        &mut long_len_sums,
                        &mut long_len_counts,
                        ec,
                        query.as_ref().map_or(0, Vec::len),
                    );
                    if let (Some(batch_sums), Some(batch_counts)) = (
                        batch_long_len_sums.get_mut(group.batch_index),
                        batch_long_len_counts.get_mut(group.batch_index),
                    ) {
                        record_long_lengths(
                            batch_sums,
                            batch_counts,
                            ec,
                            query.as_ref().map_or(0, Vec::len),
                        );
                    }
                }
                if !args.long
                    && spec.paired
                    && let (Some(ec), Some(fragment_length)) =
                        (alignment.ec.as_deref(), alignment.fragment_length)
                {
                    record_paired_fragment_length(&mut paired_flens, ec, fragment_length);
                    if let Some(batch_flens) = batch_paired_flens.get_mut(group.batch_index) {
                        record_paired_fragment_length(batch_flens, ec, fragment_length);
                    }
                }
                let pseudo_bam_placement = alignment.placement;
                let pseudo_bam_mate_placements = alignment.mate_placements.clone();
                frame_clashes = frame_clashes.saturating_add(alignment.frame_clashes);
                let ec = if novel_long_read {
                    -1
                } else {
                    alignment
                        .ec
                        .map(|ec| intern_ec(&mut ec_map, &mut ec_list, ec))
                        .unwrap_or(-1)
                };
                if ec >= 0 {
                    aligned += 1;
                    if is_unique_ec(&ec_list, ec) {
                        unique += 1;
                    }
                    let flags = if args.num {
                        u32::try_from(records - 1).unwrap_or(u32::MAX)
                    } else {
                        barcode_flag | (umi_value.flag << 8)
                    };
                    write_bus_record(
                        &mut bus_writer,
                        &BusRecord {
                            barcode,
                            umi: umi_value.binary,
                            ec,
                            count: 1,
                            flags,
                        },
                    )?;
                    if let Some(writer) = pseudo_bam.as_mut() {
                        write_pseudobam_reads(
                            writer,
                            PseudoBamBusRead {
                                read_id: records - 1,
                                ec_id: ec,
                                ec_list: &ec_list,
                                batch: &batch,
                                spec: &spec,
                                tag_present: umi_value.tag_present,
                                barcode: &barcode_seq,
                                umi: &umi_value.seq,
                                placement: pseudo_bam_placement,
                                mate_placements: &pseudo_bam_mate_placements,
                                mate_has_ec: &alignment.mate_has_ec,
                            },
                        )?;
                    }
                }
            }
        }
    }
    bus_writer.flush()?;
    let mut bus_file = bus_writer.into_inner()?;
    if header_bc_len == 0 {
        let observed_bc_len = modal_length(&bc_len_hist).or(fixed_bc_len);
        if let Some(observed_bc_len) = observed_bc_len {
            bus_file.seek(SeekFrom::Start(8))?;
            bus_file.write_all(&observed_bc_len.to_le_bytes())?;
        }
    }
    if umi_len == 0 {
        let observed_umi_len = modal_length(&umi_len_hist).or(fixed_umi_len);
        if let Some(observed_umi_len) = observed_umi_len {
            bus_file.seek(SeekFrom::Start(12))?;
            bus_file.write_all(&observed_umi_len.to_le_bytes())?;
        }
    }

    write_matrix_ec(&out_dir.join("matrix.ec"), &ec_list)?;
    if args.unmapped {
        write_unmapped_ratios(
            &out_dir.join("unmapped_ratio.txt"),
            &unmapped_ratios,
            args.batch.is_some(),
        )?;
    }
    if args.long {
        if args.batch.is_some() {
            write_batch_long_flens(
                &out_dir.join("flens.txt"),
                &index.transcript_lengths,
                &batch_long_len_sums,
                &batch_long_len_counts,
                index.k,
            )?;
        } else {
            write_long_flens(
                &out_dir.join("flens.txt"),
                &index.transcript_lengths,
                &long_len_sums,
                &long_len_counts,
                index.k,
            )?;
        }
    } else if spec.paired {
        if args.batch.is_some() {
            write_batch_paired_flens(&out_dir.join("flens.txt"), &batch_paired_flens)?;
        } else {
            write_paired_flens(&out_dir.join("flens.txt"), &paired_flens)?;
        }
    }
    if args.long || spec.paired || no_umi(&spec) {
        fs::copy(index_path, out_dir.join("index.saved"))?;
    }
    if let Some(mut writer) = novel_writer {
        writer.flush()?;
    }
    if let Some(writer) = pseudo_bam {
        writer.finish()?;
    }
    if args.batch.is_some() || synthetic_bc {
        write_cells(
            &out_dir.join("matrix.cells"),
            read_groups
                .iter()
                .enumerate()
                .map(|(idx, group)| (idx, group.id.as_deref())),
        )?;
        if synthetic_bc || args.batch_barcodes {
            let sample_barcode_len = batch_barcode_prefix_len.unwrap_or(16);
            write_sample_barcodes(
                &out_dir.join("matrix.sample.barcodes"),
                read_groups.iter().map(|group| group.batch_index as u64),
                sample_barcode_len,
            )?;
        }
    }
    write_transcripts(
        &out_dir.join("transcripts.txt"),
        &index.transcript_names,
        index.onlist.as_deref(),
    )?;
    write_run_info(
        &out_dir.join("run_info.json"),
        RunInfo {
            processed: records,
            aligned,
            unique,
            targets: onlist_target_count(&index.transcript_names, index.onlist.as_deref()),
            k: index.k,
            frame_clashes: args.aa.then_some(frame_clashes),
        },
    )?;

    println!(
        "processed: {records}, aligned: {aligned}, ecs: {}, time: {:.2?}",
        ec_list.len(),
        start.elapsed()
    );
    if let Some(max_reads) = processing.max_reads
        && records < max_reads
    {
        eprintln!(
            "Note: Number of reads processed is less than --numReads: {max_reads}, returning 1"
        );
        bail!("Number of reads processed is less than --numReads: {max_reads}, returning 1");
    }
    if aligned == 0 {
        bail!("zero reads pseudoaligned");
    }
    Ok(())
}

fn process_bam_records(
    args: &BusArgs,
    index: &kallistors::pseudoalign::BifrostIndex,
    processing: BusProcessingConfig,
    state: &mut BusWriteState<'_, impl Write>,
    unmapped_ratios: &mut Vec<f64>,
    novel_writer: &mut Option<BufWriter<File>>,
) -> Result<(u64, u64, u64, u64)> {
    let bam_path = args
        .reads
        .first()
        .ok_or_else(|| anyhow!("--bam expects exactly one BAM file"))?;
    let file = File::open(bam_path)?;
    let mut reader = noodles_bam::io::Reader::new(file);
    reader
        .read_header()
        .map_err(|err| anyhow!("failed to read BAM header: {err}"))?;

    let mut records = 0u64;
    let mut aligned = 0u64;
    let mut unique = 0u64;
    let mut frame_clashes = 0u64;

    for result in reader.records() {
        let record = result.map_err(|err| anyhow!("failed to read BAM record: {err}"))?;
        if record.flags().is_secondary() || record.flags().is_supplementary() {
            continue;
        }
        if let Some(max_reads) = processing.max_reads
            && records >= max_reads
        {
            break;
        }
        records += 1;

        let sequence = record.sequence().iter().collect::<Vec<_>>();
        let quality = record
            .quality_scores()
            .iter()
            .map(|score| score.saturating_add(b'!'))
            .collect::<Vec<_>>();
        let quality = if quality.len() == sequence.len() {
            quality
        } else {
            vec![b'!'; sequence.len()]
        };
        let Some(barcode_seq) = bam_string_tag(&record, [b'C', b'R'])
            .or_else(|| bam_string_tag(&record, [b'C', b'B']).map(normalize_corrected_barcode))
        else {
            record_bam_outputs_for_skipped_read(
                args,
                unmapped_ratios,
                novel_writer,
                index,
                &sequence,
                &quality,
                processing.options,
            )?;
            continue;
        };
        let Some(umi_seq) = bam_string_tag(&record, [b'U', b'R'])
            .or_else(|| bam_string_tag(&record, [b'R', b'X']))
            .or_else(|| bam_string_tag(&record, [b'M', b'I']).map(normalize_corrected_barcode))
            .or_else(|| bam_string_tag(&record, [b'U', b'B']).map(normalize_corrected_barcode))
        else {
            record_bam_outputs_for_skipped_read(
                args,
                unmapped_ratios,
                novel_writer,
                index,
                &sequence,
                &quality,
                processing.options,
            )?;
            continue;
        };
        let too_many_empty_kmers = if args.long || args.unmapped {
            let ratio = kallistors::pseudoalign::unmapped_kmer_ratio_bifrost(
                index,
                &sequence,
                processing.options,
            );
            if args.unmapped {
                unmapped_ratios.push(ratio);
            }
            args.long && ratio > processing.threshold
        } else {
            false
        };

        record_length(state.bc_len_hist, barcode_seq.len());
        record_length(state.umi_len_hist, umi_seq.len());

        let mut barcode_flag = 0;
        let barcode = string_to_binary(&barcode_seq, &mut barcode_flag, "barcode")?;
        let mut umi_flag = 0;
        let umi = string_to_binary(&umi_seq, &mut umi_flag, "UMI")?;

        let alignment = if args.aa {
            aa_alignment_for_query(index, &sequence, processing.options)
        } else {
            BusAlignment {
                ec: kallistors::pseudoalign::ec_for_sequence_bifrost(
                    index,
                    &sequence,
                    kallistors::pseudoalign::Strand::Unstranded,
                    processing.options,
                ),
                fragment_length: None,
                placement: None,
                mate_placements: Vec::new(),
                mate_has_ec: Vec::new(),
                frame_clashes: 0,
            }
        };
        let disjoint_intersect = alignment.ec.is_none();
        let novel_long_read = args.long && (too_many_empty_kmers || disjoint_intersect);
        if args.long
            && args.unmapped
            && disjoint_intersect
            && let Some(writer) = novel_writer.as_mut()
        {
            write_novel_read(writer, "unmapped", &sequence, &quality)?;
        }
        if novel_long_read && let Some(writer) = novel_writer.as_mut() {
            let label = if disjoint_intersect {
                "novel_disjointIntersect"
            } else {
                "novel_tooManyEmptyKmers"
            };
            write_novel_read(writer, label, &sequence, &quality)?;
        }
        frame_clashes = frame_clashes.saturating_add(alignment.frame_clashes);
        if args.long
            && !novel_long_read
            && let Some(ec) = alignment.ec.as_deref()
        {
            record_long_lengths(
                state.long_len_sums,
                state.long_len_counts,
                ec,
                sequence.len(),
            );
        }
        let ec = if novel_long_read {
            -1
        } else {
            alignment
                .ec
                .map(|ec| intern_ec(state.ec_map, state.ec_list, ec))
                .unwrap_or(-1)
        };

        if ec >= 0 {
            aligned += 1;
            if is_unique_ec(state.ec_list, ec) {
                unique += 1;
            }
            let flags = barcode_flag | (umi_flag << 8);
            write_bus_record(
                state.writer,
                &BusRecord {
                    barcode,
                    umi,
                    ec,
                    count: 1,
                    flags,
                },
            )?;
        }
    }

    Ok((records, aligned, unique, frame_clashes))
}

fn is_unique_ec(ec_list: &[Vec<u32>], ec: i32) -> bool {
    usize::try_from(ec)
        .ok()
        .and_then(|idx| ec_list.get(idx))
        .is_some_and(|ec| ec.len() == 1)
}

fn no_umi(spec: &TechnologySpec) -> bool {
    spec.umi.first().is_some_and(|spec| spec.file.is_none())
}

fn strand_specific_compatible(spec: &TechnologySpec) -> bool {
    spec.seq.len() == 1 || (spec.seq.len() == 2 && spec.paired)
}

fn pseudobam_compatible(spec: &TechnologySpec) -> bool {
    strand_specific_compatible(spec)
}

fn record_length(hist: &mut [u64; 33], len: usize) {
    hist[len.min(32)] = hist[len.min(32)].saturating_add(1);
}

fn modal_length(hist: &[u64; 33]) -> Option<u32> {
    hist.iter()
        .enumerate()
        .max_by_key(|&(len, count)| (*count, std::cmp::Reverse(len)))
        .and_then(|(len, &count)| (count > 0).then_some(len as u32))
}

fn write_pseudobam_reads(writer: &mut PseudoBamWriter, read: PseudoBamBusRead<'_>) -> Result<()> {
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

fn intern_ec(
    ec_map: &mut HashMap<Vec<u32>, i32>,
    ec_list: &mut Vec<Vec<u32>>,
    ec: Vec<u32>,
) -> i32 {
    if let Some(id) = ec_map.get(&ec) {
        *id
    } else {
        let id = i32::try_from(ec_list.len()).unwrap_or(i32::MAX);
        ec_map.insert(ec.clone(), id);
        ec_list.push(ec);
        id
    }
}

fn record_long_lengths(sums: &mut [u64], counts: &mut [u64], ec: &[u32], read_len: usize) {
    if ec.len() != 1 {
        return;
    }
    let Ok(transcript_id) = usize::try_from(ec[0]) else {
        return;
    };
    let Some(sum) = sums.get_mut(transcript_id) else {
        return;
    };
    *sum = sum.saturating_add(read_len as u64);
    if let Some(count) = counts.get_mut(transcript_id) {
        *count = count.saturating_add(1);
    }
}

fn record_paired_fragment_length(flens: &mut [u32], ec: &[u32], fragment_length: i64) {
    if ec.len() != 1 {
        return;
    }
    let Ok(idx) = usize::try_from(fragment_length) else {
        return;
    };
    let Some(bin) = flens.get_mut(idx) else {
        return;
    };
    *bin = bin.saturating_add(1);
}

fn bam_string_tag(record: &noodles_bam::Record, tag: [u8; 2]) -> Option<Vec<u8>> {
    match record.data().get(&tag)?.ok()? {
        SamDataValue::String(value) | SamDataValue::Hex(value) => {
            let bytes: &[u8] = value.as_ref();
            Some(bytes.to_vec())
        }
        _ => None,
    }
}

fn normalize_corrected_barcode(mut barcode: Vec<u8>) -> Vec<u8> {
    if let Some(pos) = barcode.iter().position(|base| *base == b'-') {
        barcode.truncate(pos);
    }
    barcode
}

fn print_technology_list() {
    println!("List of supported single-cell technologies");
    println!();
    println!("short name       description");
    println!("----------       -----------");
    println!("10xv1            10x version 1 chemistry");
    println!("10xv2            10x version 2 chemistry");
    println!("10xv3            10x version 3 chemistry");
    println!("10xv4            10x version 4 chemistry");
    println!("Bulk             Bulk RNA-seq");
    println!("ParseV3          Parse Evercode V3");
    println!("SmartSeq2        Smart-seq2 (multiplexed)");
    println!("BDWTA            BD Rhapsody WTA");
    println!("CELSeq           CEL-Seq");
    println!("CELSeq2          CEL-Seq version 2");
    println!("DropSeq          DropSeq");
    println!("inDropsv1        inDrops version 1 chemistry");
    println!("inDropsv2        inDrops version 2 chemistry");
    println!("inDropsv3        inDrops version 3 chemistry");
    println!("MATQSEQ          MATQ-SEQ");
    println!("PETRISEQ         PETRI-SEQ");
    println!("SCRBSeq          SCRB-Seq");
    println!("SmartSeq3        Smart-seq3");
    println!("SPLiT-seq        SPLiT-seq");
    println!("STORM-seq        STORM-seq");
    println!("SureCell         SureCell for ddSEQ");
    println!("VASA-seq         VASA-seq");
    println!("Visium           10x Visium Spatial Transcriptomics");
}

fn read_groups(args: &BusArgs, nfiles: usize) -> Result<Vec<ReadGroup>> {
    if nfiles == 0 {
        bail!("technology must require at least one input file");
    }
    if let Some(batch) = args.batch.as_deref() {
        parse_batch_file(batch, nfiles)
    } else if args.interleaved {
        Ok(vec![ReadGroup {
            id: None,
            files: args.reads.clone(),
            batch_index: 0,
            interleaved: true,
        }])
    } else {
        Ok(args
            .reads
            .chunks(nfiles)
            .enumerate()
            .map(|(idx, files)| ReadGroup {
                id: None,
                files: files.to_vec(),
                batch_index: idx,
                interleaved: false,
            })
            .collect())
    }
}

fn validate_read_files(args: &BusArgs, groups: &[ReadGroup]) -> Result<()> {
    let total_files = groups.iter().map(|group| group.files.len()).sum::<usize>();
    let allow_stdin = !args.bam && args.batch.is_none() && total_files == 1;
    for group in groups {
        for path in &group.files {
            if path == Path::new("-") {
                if allow_stdin {
                    continue;
                }
                bail!("file not found {}", path.display());
            }
            if !path.exists() {
                bail!("file not found {}", path.display());
            }
        }
    }
    Ok(())
}

fn report_verbose_read_groups(groups: &[ReadGroup]) {
    for (sample_idx, group) in groups.iter().enumerate() {
        let label = group
            .id
            .as_deref()
            .map_or_else(|| (sample_idx + 1).to_string(), str::to_string);
        eprintln!("[bus] will process sample {label}:");
        for path in &group.files {
            eprintln!("[bus]   {}", path.display());
        }
    }
}

fn next_parallel_batch<R: ReadSource>(
    readers: &mut [R],
) -> Result<Vec<kallistors::io::FastqRecord>> {
    let mut batch = Vec::with_capacity(readers.len());
    for reader in readers {
        match reader.next_record() {
            Some(record) => batch.push(record.map_err(|err| anyhow!("FASTQ read failed: {err}"))?),
            None => {
                if batch.is_empty() {
                    break;
                }
                bail!("FASTQ files ended at different record counts");
            }
        }
    }
    Ok(batch)
}

fn next_interleaved_batch<R: ReadSource>(
    readers: &mut [R],
    nfiles: usize,
) -> Result<Vec<kallistors::io::FastqRecord>> {
    let Some(reader) = readers.first_mut() else {
        return Ok(Vec::new());
    };
    let mut batch = Vec::with_capacity(nfiles);
    for _ in 0..nfiles {
        match reader.next_record() {
            Some(record) => batch.push(record.map_err(|err| anyhow!("FASTQ read failed: {err}"))?),
            None => {
                if batch.is_empty() {
                    break;
                }
                bail!("interleaved FASTQ ended in the middle of a record group");
            }
        }
    }
    Ok(batch)
}

fn parse_batch_file(path: &Path, nfiles: usize) -> Result<Vec<ReadGroup>> {
    let text = fs::read_to_string(path)?;
    let mut groups = Vec::new();
    let mut batch_ids = HashMap::<String, usize>::new();
    let batch_dir = path.parent().unwrap_or_else(|| Path::new("."));
    for (line_idx, line) in text.lines().enumerate() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let parts = line.split_whitespace().collect::<Vec<_>>();
        if parts.len() != nfiles + 1 {
            bail!(
                "batch file line {} has {} files, expected {}",
                line_idx + 1,
                parts.len().saturating_sub(1),
                nfiles
            );
        }
        let id = parts[0].to_string();
        let next_batch_index = batch_ids.len();
        let batch_index = *batch_ids.entry(id.clone()).or_insert(next_batch_index);
        groups.push(ReadGroup {
            id: Some(id),
            files: parts[1..]
                .iter()
                .map(|part| resolve_batch_read_path(batch_dir, part))
                .collect(),
            batch_index,
            interleaved: false,
        });
    }
    if groups.is_empty() {
        bail!("batch file contains no read groups");
    }
    Ok(groups)
}

fn should_infer_bulk_paired(args: &BusArgs) -> Result<bool> {
    if let Some(batch) = args.batch.as_deref() {
        return Ok(infer_batch_nfiles(batch)? == 2);
    }
    if args.interleaved {
        return Ok(true);
    }
    Ok(args.reads.len() > 1 && args.reads.len().is_multiple_of(2))
}

fn infer_batch_nfiles(path: &Path) -> Result<usize> {
    let text = fs::read_to_string(path)?;
    let mut inferred = None;
    for (line_idx, line) in text.lines().enumerate() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let parts = line.split_whitespace().collect::<Vec<_>>();
        let files = parts.len().saturating_sub(1);
        if files != 1 && files != 2 {
            bail!(
                "batch file line {} has {} files, expected 1 or 2",
                line_idx + 1,
                files
            );
        }
        match inferred {
            Some(expected) if expected != files => {
                bail!(
                    "batch file line {} has {} files, expected {}",
                    line_idx + 1,
                    files,
                    expected
                );
            }
            Some(_) => {}
            None => inferred = Some(files),
        }
    }
    inferred.ok_or_else(|| anyhow!("batch file contains no read groups"))
}

fn resolve_batch_read_path(batch_dir: &Path, value: &str) -> PathBuf {
    let path = PathBuf::from(value);
    if value == "-" || path.is_absolute() || path.exists() {
        path
    } else {
        batch_dir.join(path)
    }
}

fn barcode_sequence(
    batch: &[kallistors::io::FastqRecord],
    spec: &TechnologySpec,
    batch_index: usize,
    batch_barcodes: bool,
) -> Result<Option<Vec<u8>>> {
    let extracted = if spec.bc.first().is_some_and(|v| v.file.is_none()) {
        fake_barcode(batch_index as u64, 16)
    } else {
        let Some(barcode) = concat_slices_optional(batch, &spec.bc, false)? else {
            return Ok(None);
        };
        barcode
    };
    if batch_barcodes && spec.bc.first().is_none_or(|v| v.file.is_some()) {
        if extracted.len() >= 32 {
            bail!("--batch-barcodes requires barcode length shorter than 32 bases");
        }
        let prefix_len = 32usize.saturating_sub(extracted.len()).min(32);
        let mut prefixed = fake_barcode(batch_index as u64, prefix_len);
        prefixed.extend_from_slice(&extracted);
        Ok(Some(prefixed))
    } else {
        Ok(Some(extracted))
    }
}

fn configure_tag(
    technology: &str,
    tag_sequence: Option<&str>,
    spec: &mut TechnologySpec,
) -> Result<Option<ConfiguredTag>> {
    let inferred_default = tag_sequence.is_none() && technology.eq_ignore_ascii_case("SMARTSEQ3");
    let tag = tag_sequence
        .map(str::to_string)
        .or_else(|| inferred_default.then(|| "ATTGCGCAATG".to_string()));
    let Some(tag) = tag else {
        return Ok(None);
    };
    validate_tag_sequence(&tag)?;
    let Some(first_umi) = spec.umi.first_mut() else {
        bail!("tag sequence requires a UMI slice");
    };
    let Some(stop) = first_umi.stop else {
        bail!("tag sequence requires a bounded UMI slice");
    };
    first_umi.start = first_umi.start.saturating_add(tag.len());
    if first_umi.start >= stop {
        bail!("tag sequence must be shorter than UMI sequence");
    }
    spec.tag_len = Some(tag.len());
    let mut flag = 0;
    Ok(Some(ConfiguredTag {
        config: TagConfig {
            binary: string_to_binary(tag.as_bytes(), &mut flag, "tag sequence")?,
            sequence: tag.into_bytes(),
        },
        inferred_default,
    }))
}

fn validate_tag_sequence(tag: &str) -> Result<()> {
    for (idx, base) in tag.bytes().enumerate() {
        if !matches!(base, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't') {
            bail!(
                "tag sequence contains invalid base {} at position {}",
                char::from(base),
                idx + 1
            );
        }
    }
    Ok(())
}

fn umi_value(
    batch: &[kallistors::io::FastqRecord],
    spec: &TechnologySpec,
    tag: Option<&TagConfig>,
) -> Result<Option<UmiValue>> {
    if spec.bulk_like() {
        let mut flag = 0;
        return Ok(Some(UmiValue {
            binary: string_to_binary(b"A", &mut flag, "UMI")?,
            flag,
            len: 1,
            seq: b"A".to_vec(),
            tag_present: false,
            disable_strand_specific: false,
        }));
    }
    if spec.keep_fastq_comments {
        let Some(umi_seq) = rx_umi_from_header(
            batch
                .first()
                .ok_or_else(|| anyhow!("RX UMI extraction requires at least one FASTQ record"))?,
        ) else {
            return Ok(None);
        };
        let mut flag = 0;
        return Ok(Some(UmiValue {
            binary: string_to_binary(&umi_seq, &mut flag, "UMI")?,
            flag,
            len: umi_seq.len(),
            seq: umi_seq,
            tag_present: false,
            disable_strand_specific: false,
        }));
    }

    let Some(tag) = tag else {
        let Some(umi_seq) = concat_slices_optional(batch, &spec.umi, false)? else {
            return Ok(None);
        };
        let mut flag = 0;
        return Ok(Some(UmiValue {
            binary: string_to_binary(&umi_seq, &mut flag, "UMI")?,
            flag,
            len: umi_seq.len(),
            seq: umi_seq,
            tag_present: false,
            disable_strand_specific: false,
        }));
    };

    let Some(first_umi) = spec.umi.first() else {
        bail!("tag sequence requires a UMI slice");
    };
    let Some(mut tagged_umi) = extract_tagged_umi_optional(batch, *first_umi, tag.sequence.len())?
    else {
        return Ok(None);
    };
    for spec in spec.umi.iter().skip(1) {
        let Some(slice) = extract_slice_optional(batch, *spec)? else {
            return Ok(None);
        };
        tagged_umi.extend_from_slice(&slice);
    }
    let mut flag = 0;
    let tagged_binary = string_to_binary(&tagged_umi, &mut flag, "tagged UMI")?;
    let tag_len = tag.sequence.len();
    let hamming_threshold = usize::from(tag_len > 5);
    let tag_present = tagged_umi.len() >= tag_len
        && hamming(
            tag.binary,
            tagged_binary >> (2 * (tagged_umi.len() - tag_len)),
            tag_len,
        ) <= hamming_threshold;
    if tag_present {
        let umi_len = tagged_umi.len() - tag_len;
        let binary = if umi_len == 0 {
            0
        } else {
            tagged_binary & !(u64::MAX << (2 * umi_len))
        };
        Ok(Some(UmiValue {
            binary,
            flag,
            len: umi_len,
            seq: tagged_umi[tag_len..].to_vec(),
            tag_present: true,
            disable_strand_specific: false,
        }))
    } else {
        Ok(Some(UmiValue {
            binary: u64::MAX,
            flag: 0,
            len: 0,
            seq: Vec::new(),
            tag_present: false,
            disable_strand_specific: true,
        }))
    }
}

fn extract_tagged_umi_optional(
    records: &[kallistors::io::FastqRecord],
    spec: SliceSpec,
    tag_len: usize,
) -> Result<Option<Vec<u8>>> {
    let mut tagged_spec = spec;
    tagged_spec.start = tagged_spec.start.saturating_sub(tag_len);
    extract_slice_optional(records, tagged_spec)
}

fn rx_umi_from_header(record: &kallistors::io::FastqRecord) -> Option<Vec<u8>> {
    let header = record.header.strip_prefix(b"@").unwrap_or(&record.header);
    let pos = header.windows(5).position(|window| window == b"RX:Z:")?;
    let start = pos + 5;
    let end = header[start..]
        .iter()
        .position(|base| matches!(base, b' ' | b'\t'))
        .map_or(header.len(), |offset| start + offset);
    (end > start).then(|| header[start..end].to_vec())
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

fn record_unmapped_ratio_for_skipped_read(
    enabled: bool,
    ratios: &mut Vec<f64>,
    index: &kallistors::pseudoalign::BifrostIndex,
    batch: &[kallistors::io::FastqRecord],
    spec: &TechnologySpec,
    options: kallistors::pseudoalign::PseudoalignOptions,
) -> Result<()> {
    if enabled {
        let query = bus_query_sequence(batch, spec, false)?;
        ratios.push(kallistors::pseudoalign::unmapped_kmer_ratio_bifrost(
            index, &query, options,
        ));
    }
    Ok(())
}

fn record_bam_outputs_for_skipped_read(
    args: &BusArgs,
    ratios: &mut Vec<f64>,
    novel_writer: &mut Option<BufWriter<File>>,
    index: &kallistors::pseudoalign::BifrostIndex,
    sequence: &[u8],
    quality: &[u8],
    options: kallistors::pseudoalign::PseudoalignOptions,
) -> Result<()> {
    if args.unmapped {
        ratios.push(kallistors::pseudoalign::unmapped_kmer_ratio_bifrost(
            index, sequence, options,
        ));
    }
    if args.long
        && args.unmapped
        && let Some(writer) = novel_writer.as_mut()
    {
        write_novel_read(writer, "skipped_missingTags", sequence, quality)?;
    }
    Ok(())
}

fn pseudoalign_bus_read(
    index: &kallistors::pseudoalign::BifrostIndex,
    batch: &[kallistors::io::FastqRecord],
    spec: &TechnologySpec,
    options: kallistors::pseudoalign::PseudoalignOptions,
    tag_present: bool,
    aa: bool,
) -> Result<BusAlignment> {
    if spec.paired {
        if spec.seq.len() != 2 {
            bail!("paired BUS technology requires exactly two sequence slices");
        }
        let left = extract_sequence_slice(batch, spec.seq[0], spec, tag_present)?;
        let right = extract_sequence_slice(batch, spec.seq[1], spec, tag_present)?;
        let trace = kallistors::pseudoalign::trace_read_pair_bifrost(
            index,
            &left,
            &right,
            kallistors::pseudoalign::Strand::Unstranded,
            options,
        );
        Ok(BusAlignment {
            ec: trace.merged_ec,
            fragment_length: trace.estimated_fragment_length,
            placement: trace.placement,
            mate_placements: vec![trace.left.placement, trace.right.placement],
            mate_has_ec: vec![
                trace.left.ec_after_strand_filter.is_some(),
                trace.right.ec_after_strand_filter.is_some(),
            ],
            frame_clashes: 0,
        })
    } else {
        let query = concat_sequence_slices(batch, spec, tag_present)?;
        if aa {
            return Ok(aa_alignment_for_query(index, &query, options));
        }
        let trace = kallistors::pseudoalign::trace_read_bifrost(
            index,
            &query,
            kallistors::pseudoalign::Strand::Unstranded,
            None,
            options,
        );
        Ok(BusAlignment {
            ec: trace.ec_after_strand_filter,
            fragment_length: None,
            placement: trace.placement,
            mate_placements: Vec::new(),
            mate_has_ec: Vec::new(),
            frame_clashes: 0,
        })
    }
}

fn aa_alignment_for_query(
    index: &kallistors::pseudoalign::BifrostIndex,
    query: &[u8],
    options: kallistors::pseudoalign::PseudoalignOptions,
) -> BusAlignment {
    let mut best_ec: Option<Vec<u32>> = None;
    let mut best_cardinality = usize::MAX;
    let mut frame_clashes = 0u64;
    let mut frame = Vec::new();
    for seq in [query.to_vec(), kallistors::util::reverse_complement(query)] {
        for offset in 0..3 {
            if offset >= seq.len() {
                continue;
            }
            frame.clear();
            frame.extend_from_slice(&kallistors::util::nucleotide_to_comma_free(&seq[offset..]));
            if frame.len() < index.k {
                continue;
            }
            if let Some(ec) = kallistors::pseudoalign::ec_for_sequence_bifrost(
                index,
                &frame,
                kallistors::pseudoalign::Strand::Unstranded,
                options,
            ) {
                match ec.len().cmp(&best_cardinality) {
                    std::cmp::Ordering::Less => {
                        best_cardinality = ec.len();
                        best_ec = Some(ec);
                    }
                    std::cmp::Ordering::Equal => {
                        frame_clashes = frame_clashes.saturating_add(1);
                    }
                    std::cmp::Ordering::Greater => {}
                }
            }
        }
    }
    BusAlignment {
        ec: best_ec,
        fragment_length: None,
        placement: None,
        mate_placements: Vec::new(),
        mate_has_ec: Vec::new(),
        frame_clashes,
    }
}

fn bus_query_sequence(
    batch: &[kallistors::io::FastqRecord],
    spec: &TechnologySpec,
    tag_present: bool,
) -> Result<Vec<u8>> {
    if spec.paired {
        if spec.seq.len() != 2 {
            bail!("paired BUS technology requires exactly two sequence slices");
        }
        let mut query = extract_sequence_slice(batch, spec.seq[0], spec, tag_present)?;
        query.push(b'N');
        query.extend_from_slice(&extract_sequence_slice(
            batch,
            spec.seq[1],
            spec,
            tag_present,
        )?);
        Ok(query)
    } else {
        concat_sequence_slices(batch, spec, tag_present)
    }
}

fn bus_query_quality(
    batch: &[kallistors::io::FastqRecord],
    spec: &TechnologySpec,
    tag_present: bool,
) -> Result<Vec<u8>> {
    if spec.paired {
        if spec.seq.len() != 2 {
            bail!("paired BUS technology requires exactly two sequence slices");
        }
        let mut quality = extract_quality_slice(batch, spec.seq[0], spec, tag_present)?;
        quality.push(b'!');
        quality.extend_from_slice(&extract_quality_slice(
            batch,
            spec.seq[1],
            spec,
            tag_present,
        )?);
        Ok(quality)
    } else {
        concat_quality_slices(batch, spec, tag_present)
    }
}

fn technology_spec(name: &str) -> Result<TechnologySpec> {
    let (base, suffix) = name
        .split_once('%')
        .map_or((name, None), |(base, suffix)| (base, Some(suffix)));
    let upper = base.to_ascii_uppercase();
    let mut spec = match upper.as_str() {
        "BULK" => tech(
            1,
            vec![slice(0, 0, None)],
            vec![sentinel_slice()],
            vec![sentinel_slice()],
            false,
        ),
        "10XV1" => tech(
            3,
            vec![slice(2, 0, None)],
            vec![slice(1, 0, Some(10))],
            vec![slice(0, 0, Some(14))],
            false,
        )
        .with_default_strand(kallistors::pseudoalign::StrandSpecific::FR),
        "10XV2" => tech(
            2,
            vec![slice(1, 0, None)],
            vec![slice(0, 16, Some(26))],
            vec![slice(0, 0, Some(16))],
            false,
        )
        .with_default_strand(kallistors::pseudoalign::StrandSpecific::FR),
        "10XV3" | "10XV4" | "VISIUM" => tech(
            2,
            vec![slice(1, 0, None)],
            vec![slice(0, 16, Some(28))],
            vec![slice(0, 0, Some(16))],
            false,
        )
        .with_default_strand(kallistors::pseudoalign::StrandSpecific::FR),
        "SMARTSEQ3" => tech(
            4,
            vec![slice(2, 22, None), slice(3, 0, None)],
            vec![slice(2, 0, Some(19))],
            vec![slice(0, 0, None), slice(1, 0, None)],
            true,
        )
        .with_default_strand(kallistors::pseudoalign::StrandSpecific::FR),
        "SMARTSEQ2" => tech(
            3,
            vec![slice(2, 0, None)],
            vec![sentinel_slice()],
            vec![slice(0, 0, None), slice(1, 0, None)],
            false,
        ),
        "DROPSEQ" => tech(
            2,
            vec![slice(1, 0, None)],
            vec![slice(0, 12, Some(20))],
            vec![slice(0, 0, Some(12))],
            false,
        ),
        "CELSEQ" => tech(
            2,
            vec![slice(1, 0, None)],
            vec![slice(0, 8, Some(12))],
            vec![slice(0, 0, Some(8))],
            false,
        )
        .with_default_strand(kallistors::pseudoalign::StrandSpecific::FR),
        "CELSEQ2" => tech(
            2,
            vec![slice(1, 0, None)],
            vec![slice(0, 0, Some(6))],
            vec![slice(0, 6, Some(12))],
            false,
        )
        .with_default_strand(kallistors::pseudoalign::StrandSpecific::FR),
        "SCRBSEQ" => tech(
            2,
            vec![slice(1, 0, None)],
            vec![slice(0, 6, Some(16))],
            vec![slice(0, 0, Some(6))],
            false,
        ),
        "SURECELL" => tech(
            2,
            vec![slice(1, 0, None)],
            vec![slice(0, 51, Some(59))],
            vec![
                slice(0, 0, Some(6)),
                slice(0, 21, Some(27)),
                slice(0, 42, Some(48)),
            ],
            false,
        )
        .with_default_strand(kallistors::pseudoalign::StrandSpecific::FR),
        "INDROPSV1" => tech(
            2,
            vec![slice(1, 0, None)],
            vec![slice(0, 42, Some(48))],
            vec![slice(0, 0, Some(11)), slice(0, 30, Some(38))],
            false,
        ),
        "INDROPSV2" => tech(
            2,
            vec![slice(0, 0, None)],
            vec![slice(1, 42, Some(48))],
            vec![slice(1, 0, Some(11)), slice(1, 30, Some(38))],
            false,
        ),
        "INDROPSV3" => tech(
            3,
            vec![slice(2, 0, None)],
            vec![slice(1, 8, Some(14))],
            vec![slice(0, 0, Some(8)), slice(1, 0, Some(8))],
            false,
        ),
        "PARSEV3" => tech(
            2,
            vec![slice(0, 0, None)],
            vec![slice(1, 0, Some(10))],
            vec![
                slice(1, 10, Some(18)),
                slice(1, 30, Some(38)),
                slice(1, 50, Some(58)),
            ],
            false,
        )
        .with_default_strand(kallistors::pseudoalign::StrandSpecific::FR),
        "PETRISEQ" => tech(
            2,
            vec![slice(1, 0, Some(17))],
            vec![slice(0, 0, Some(7))],
            vec![
                slice(0, 7, Some(14)),
                slice(0, 29, Some(36)),
                slice(0, 50, Some(58)),
            ],
            false,
        ),
        "SPLIT-SEQ" => tech(
            2,
            vec![slice(0, 0, None)],
            vec![slice(1, 0, Some(10))],
            vec![
                slice(1, 10, Some(18)),
                slice(1, 48, Some(56)),
                slice(1, 78, Some(86)),
            ],
            false,
        )
        .with_default_strand(kallistors::pseudoalign::StrandSpecific::FR),
        "STORM-SEQ" => tech(
            2,
            vec![slice(0, 0, None), slice(1, 14, None)],
            vec![slice(1, 0, Some(8))],
            vec![sentinel_slice()],
            true,
        )
        .with_default_strand(kallistors::pseudoalign::StrandSpecific::RF),
        "MATQSEQ" => tech(
            1,
            vec![slice(0, 8, None)],
            vec![sentinel_slice()],
            vec![slice(0, 0, Some(8))],
            false,
        ),
        "BDWTA" => tech(
            2,
            vec![slice(1, 0, None)],
            vec![slice(0, 52, Some(60))],
            vec![
                slice(0, 0, Some(9)),
                slice(0, 21, Some(30)),
                slice(0, 43, Some(52)),
            ],
            false,
        )
        .with_default_strand(kallistors::pseudoalign::StrandSpecific::FR),
        "VASA-SEQ" => tech(
            1,
            vec![slice(0, 14, None)],
            vec![slice(0, 0, Some(6))],
            vec![slice(0, 6, Some(14))],
            false,
        )
        .with_default_strand(kallistors::pseudoalign::StrandSpecific::FR),
        other => parse_custom_technology(other)?,
    };
    apply_technology_suffix(&mut spec, suffix)?;
    Ok(spec)
}

fn apply_technology_suffix(spec: &mut TechnologySpec, suffix: Option<&str>) -> Result<()> {
    let Some(suffix) = suffix else {
        return Ok(());
    };
    let mut fields = suffix.split('%');
    if let Some(strand) = fields.next() {
        let strand = strand.to_ascii_uppercase();
        if strand.is_empty() || strand == "NONE" {
        } else if strand.starts_with("FORWARD") {
            spec.default_strand = Some(kallistors::pseudoalign::StrandSpecific::FR);
        } else if strand.starts_with("REVERSE") {
            spec.default_strand = Some(kallistors::pseudoalign::StrandSpecific::RF);
        } else {
            bail!("invalid technology strand suffix: {strand}");
        }
    }
    if let Some(parity) = fields.next() {
        let parity = parity.to_ascii_uppercase();
        if parity.is_empty() || parity == "NONE" {
        } else if parity.starts_with("PAIRED") {
            if !spec.paired {
                mark_technology_paired(spec, "%PAIRED")?;
            }
        } else {
            bail!("invalid technology pairing suffix: {parity}");
        }
    }
    if let Some(extra) = fields.next() {
        bail!("unexpected technology suffix field: {extra}");
    }
    Ok(())
}

fn mark_technology_paired(spec: &mut TechnologySpec, source: &str) -> Result<()> {
    match spec.seq.len() {
        1 => {
            let file = spec.nfiles;
            spec.seq.push(slice(file, 0, None));
            spec.nfiles += 1;
            spec.paired = true;
        }
        2 => {
            spec.paired = true;
        }
        _ => bail!("{source} requires a technology with one or two sequence slices"),
    }
    Ok(())
}

fn add_unpaired_mate_slice(spec: &mut TechnologySpec, source: &str) -> Result<()> {
    match spec.seq.len() {
        1 => {
            let file = spec.nfiles;
            spec.seq.push(slice(file, 0, None));
            spec.nfiles += 1;
        }
        2 => {}
        _ => bail!("{source} requires a technology with one or two sequence slices"),
    }
    Ok(())
}

fn technology_base_upper(name: &str) -> String {
    name.split_once('%')
        .map_or(name, |(base, _)| base)
        .to_ascii_uppercase()
}

fn tech(
    nfiles: usize,
    seq: Vec<SliceSpec>,
    umi: Vec<SliceSpec>,
    bc: Vec<SliceSpec>,
    paired: bool,
) -> TechnologySpec {
    TechnologySpec {
        nfiles,
        seq,
        bc,
        umi,
        paired,
        keep_fastq_comments: false,
        default_strand: None,
        tag_len: None,
    }
}

fn slice(file: usize, start: usize, stop: Option<usize>) -> SliceSpec {
    SliceSpec {
        file: Some(file),
        start,
        stop,
    }
}

fn sentinel_slice() -> SliceSpec {
    SliceSpec {
        file: None,
        start: 0,
        stop: None,
    }
}

fn parse_custom_technology(name: &str) -> Result<TechnologySpec> {
    let colon_count = name.as_bytes().iter().filter(|&&byte| byte == b':').count();
    if colon_count != 2 {
        let detail = match colon_count {
            0 => "none found".to_string(),
            1 => "only one found".to_string(),
            3 => "three found".to_string(),
            count => format!("{count} found"),
        };
        bail!("custom technology must contain two colons (:), {detail}: \"{name}\"");
    }
    let parts = name.split(':').collect::<Vec<_>>();
    let bc = parse_slice_list(parts[0])?;
    if bc.is_empty() {
        bail!("custom technology barcode list is empty");
    }
    let mut keep_fastq_comments = false;
    let umi = if parts[1].eq_ignore_ascii_case("RX") {
        keep_fastq_comments = true;
        vec![sentinel_slice()]
    } else {
        parse_slice_list(parts[1])?
    };
    if umi.is_empty() {
        bail!("custom technology UMI list is empty");
    }
    let seq = parse_slice_list(parts[2])?;
    if seq.is_empty() {
        bail!("custom technology sequence list is empty");
    }
    let nfiles = bc
        .iter()
        .chain(umi.iter())
        .chain(seq.iter())
        .filter_map(|spec| spec.file)
        .max()
        .map_or(1, |file| file + 1);
    Ok(TechnologySpec {
        nfiles,
        paired: false,
        keep_fastq_comments,
        seq,
        bc,
        umi,
        default_strand: None,
        tag_len: None,
    })
}

fn parse_slice_list(value: &str) -> Result<Vec<SliceSpec>> {
    if value.is_empty() {
        return Ok(Vec::new());
    }
    let nums = value
        .split(',')
        .map(|part| {
            part.parse::<i32>()
                .map_err(|_| anyhow!("invalid technology integer {part:?}"))
        })
        .collect::<Result<Vec<_>>>()?;
    if nums.len() % 3 != 0 {
        bail!("technology slice lists must contain triples of file,start,stop");
    }
    let mut specs = Vec::with_capacity(nums.len() / 3);
    for triple in nums.chunks_exact(3) {
        let file = triple[0];
        let start = triple[1];
        let stop = triple[2];
        if file < -1 {
            bail!("invalid technology file number {file}");
        }
        if file == -1 {
            specs.push(sentinel_slice());
            continue;
        }
        if start < 0 {
            bail!("invalid technology start {start}");
        }
        if stop != 0 && stop <= start {
            bail!("invalid technology stop {stop} before start {start}");
        }
        specs.push(SliceSpec {
            file: Some(usize::try_from(file).map_err(|_| anyhow!("invalid file number {file}"))?),
            start: usize::try_from(start).map_err(|_| anyhow!("invalid start {start}"))?,
            stop: (stop != 0)
                .then(|| usize::try_from(stop).map_err(|_| anyhow!("invalid stop {stop}")))
                .transpose()?,
        });
    }
    Ok(specs)
}

fn total_len(specs: &[SliceSpec]) -> Option<u32> {
    let mut total = 0u32;
    for spec in specs {
        spec.file?;
        let stop = spec.stop?;
        total = total.saturating_add(u32::try_from(stop.saturating_sub(spec.start)).ok()?);
    }
    Some(total)
}

#[cfg(test)]
fn concat_slices(
    records: &[kallistors::io::FastqRecord],
    specs: &[SliceSpec],
    separate_with_n: bool,
) -> Result<Vec<u8>> {
    let mut out = Vec::new();
    for (idx, spec) in specs.iter().enumerate() {
        if idx > 0 && separate_with_n {
            out.push(b'N');
        }
        out.extend_from_slice(&extract_slice(records, *spec)?);
    }
    Ok(out)
}

fn concat_slices_optional(
    records: &[kallistors::io::FastqRecord],
    specs: &[SliceSpec],
    separate_with_n: bool,
) -> Result<Option<Vec<u8>>> {
    let mut out = Vec::new();
    for (idx, spec) in specs.iter().enumerate() {
        if idx > 0 && separate_with_n {
            out.push(b'N');
        }
        let Some(slice) = extract_slice_optional(records, *spec)? else {
            return Ok(None);
        };
        out.extend_from_slice(&slice);
    }
    Ok(Some(out))
}

fn concat_sequence_slices(
    records: &[kallistors::io::FastqRecord],
    spec: &TechnologySpec,
    tag_present: bool,
) -> Result<Vec<u8>> {
    let mut out = Vec::new();
    for (idx, seq_spec) in spec.seq.iter().enumerate() {
        if idx > 0 {
            out.push(b'N');
        }
        out.extend_from_slice(&extract_sequence_slice(
            records,
            *seq_spec,
            spec,
            tag_present,
        )?);
    }
    Ok(out)
}

fn concat_quality_slices(
    records: &[kallistors::io::FastqRecord],
    spec: &TechnologySpec,
    tag_present: bool,
) -> Result<Vec<u8>> {
    let mut out = Vec::new();
    for (idx, seq_spec) in spec.seq.iter().enumerate() {
        if idx > 0 {
            out.push(b'!');
        }
        out.extend_from_slice(&extract_quality_slice(
            records,
            *seq_spec,
            spec,
            tag_present,
        )?);
    }
    Ok(out)
}

fn extract_sequence_slice(
    records: &[kallistors::io::FastqRecord],
    seq_spec: SliceSpec,
    spec: &TechnologySpec,
    tag_present: bool,
) -> Result<Vec<u8>> {
    let Some(umi_spec) = spec.umi.first() else {
        return extract_slice(records, seq_spec);
    };
    if !tag_present
        && spec.tag_len.is_some()
        && umi_spec.file == seq_spec.file
        && umi_spec.file.is_some()
    {
        let mut adjusted = seq_spec;
        adjusted.start = umi_spec.start.saturating_sub(spec.tag_len.unwrap_or(0));
        adjusted.stop = seq_spec.stop;
        extract_slice(records, adjusted)
    } else {
        extract_slice(records, seq_spec)
    }
}

fn extract_quality_slice(
    records: &[kallistors::io::FastqRecord],
    seq_spec: SliceSpec,
    spec: &TechnologySpec,
    tag_present: bool,
) -> Result<Vec<u8>> {
    let Some(umi_spec) = spec.umi.first() else {
        return extract_slice_from(
            records,
            seq_spec,
            |record| record.qual.as_slice(),
            "quality",
        );
    };
    if !tag_present
        && spec.tag_len.is_some()
        && umi_spec.file == seq_spec.file
        && umi_spec.file.is_some()
    {
        let mut adjusted = seq_spec;
        adjusted.start = umi_spec.start.saturating_sub(spec.tag_len.unwrap_or(0));
        adjusted.stop = seq_spec.stop;
        extract_slice_from(
            records,
            adjusted,
            |record| record.qual.as_slice(),
            "quality",
        )
    } else {
        extract_slice_from(
            records,
            seq_spec,
            |record| record.qual.as_slice(),
            "quality",
        )
    }
}

fn extract_slice(records: &[kallistors::io::FastqRecord], spec: SliceSpec) -> Result<Vec<u8>> {
    extract_slice_from(records, spec, |record| record.seq.as_slice(), "read")
}

fn extract_slice_optional(
    records: &[kallistors::io::FastqRecord],
    spec: SliceSpec,
) -> Result<Option<Vec<u8>>> {
    extract_slice_from_optional(records, spec, |record| record.seq.as_slice())
}

fn extract_slice_from(
    records: &[kallistors::io::FastqRecord],
    spec: SliceSpec,
    get: impl Fn(&kallistors::io::FastqRecord) -> &[u8],
    label: &str,
) -> Result<Vec<u8>> {
    let Some(file) = spec.file else {
        return Ok(Vec::new());
    };
    let seq = get(records
        .get(file)
        .ok_or_else(|| anyhow!("technology references missing file {}", file + 1))?);
    let end = spec.stop.unwrap_or(seq.len());
    if spec.start > end || end > seq.len() {
        bail!(
            "{label} length {} is too short for slice {}:{}",
            seq.len(),
            spec.start,
            end
        );
    }
    Ok(seq[spec.start..end].to_vec())
}

fn extract_slice_from_optional(
    records: &[kallistors::io::FastqRecord],
    spec: SliceSpec,
    get: impl Fn(&kallistors::io::FastqRecord) -> &[u8],
) -> Result<Option<Vec<u8>>> {
    let Some(file) = spec.file else {
        return Ok(Some(Vec::new()));
    };
    let seq = get(records
        .get(file)
        .ok_or_else(|| anyhow!("technology references missing file {}", file + 1))?);
    let end = spec.stop.unwrap_or(seq.len());
    if spec.start > end || end > seq.len() {
        return Ok(None);
    }
    Ok(Some(seq[spec.start..end].to_vec()))
}

fn onlist_target_count(names: &[String], onlist: Option<&[bool]>) -> usize {
    match onlist {
        Some(onlist) => names
            .iter()
            .enumerate()
            .filter(|(idx, _)| onlist.get(*idx).copied().unwrap_or(false))
            .count(),
        None => names.len(),
    }
}

fn write_gene_list(path: &Path, model: &GenomeModel) -> Result<()> {
    let mut writer = BufWriter::new(File::create(path)?);
    for (idx, gene) in model.genes.iter().enumerate() {
        writeln!(writer, "{idx}\t{}\t{}", gene.id, gene.name)?;
    }
    Ok(())
}

fn string_to_binary(seq: &[u8], flag: &mut u32, label: &str) -> Result<u64> {
    if seq.len() > 32 {
        bail!("{label} length {} exceeds BUS limit of 32 bases", seq.len());
    }
    *flag = 0;
    let mut result = 0u64;
    let mut num_n = 0u32;
    let mut pos_n = 0u32;
    for (idx, base) in seq.iter().enumerate() {
        let x = (base & 4) >> 1;
        if (base & 3) == 2 {
            if num_n == 0 {
                pos_n = u32::try_from(idx).unwrap_or(0);
            }
            num_n += 1;
        }
        result <<= 2;
        result |= u64::from(x + ((x ^ (base & 2)) >> 1));
    }
    if num_n > 0 {
        *flag = (num_n.min(3) & 3) | ((pos_n & 31) << 2);
    }
    Ok(result)
}

fn hamming(left: u64, right: u64, len: usize) -> usize {
    let diff = left ^ right;
    let mut distance = 0usize;
    for idx in 0..len {
        let shift = 2 * (len - idx - 1);
        if ((diff >> shift) & 0x03) != 0 {
            distance += 1;
        }
    }
    distance
}

#[cfg(test)]
mod tests {
    use super::{
        BusArgs, bus_pseudoalign_options, bus_query_quality, concat_sequence_slices, concat_slices,
        mark_technology_paired, onlist_target_count, read_groups, string_to_binary,
        technology_spec, total_len, write_batch_long_flens, write_batch_paired_flens,
        write_matrix_ec, write_transcripts,
    };
    use kallistors::io::FastqRecord;
    use kallistors::pseudoalign::StrandSpecific;

    #[test]
    fn bus_binary_encoding_matches_kallisto_order() {
        let mut flag = 99;
        assert_eq!(
            string_to_binary(b"ACGT", &mut flag, "test").unwrap(),
            0b00_01_10_11
        );
        assert_eq!(flag, 0);
    }

    #[test]
    fn tenx_v3_lengths_match_preset() {
        let spec = technology_spec("10xv3").unwrap();
        assert_eq!(spec.nfiles, 2);
        assert_eq!(total_len(&spec.bc), Some(16));
        assert_eq!(total_len(&spec.umi), Some(12));
        assert_eq!(spec.default_strand, Some(StrandSpecific::FR));
    }

    #[test]
    fn bulk_is_single_end_by_default_and_pairable() {
        let mut spec = technology_spec("Bulk").unwrap();
        assert_eq!(spec.nfiles, 1);
        assert_eq!(spec.seq.len(), 1);
        assert_eq!(total_len(&spec.bc), None);
        assert_eq!(total_len(&spec.umi), None);
        assert!(!spec.paired);

        mark_technology_paired(&mut spec, "--paired").unwrap();

        assert_eq!(spec.nfiles, 2);
        assert_eq!(spec.seq.len(), 2);
        assert!(spec.paired);
    }

    #[test]
    fn read_groups_rejects_zero_input_file_technologies() {
        let args = BusArgs {
            reads: vec!["reads.fastq".into()],
            ..BusArgs::default()
        };
        let err = match read_groups(&args, 0) {
            Ok(_) => panic!("zero-input-file technology should fail"),
            Err(err) => err.to_string(),
        };
        assert!(err.contains("at least one input file"), "{err}");
    }

    #[test]
    fn bus_forwards_dfk_onlist_to_pseudoalign_options() {
        let args = BusArgs {
            dfk_onlist: true,
            ..BusArgs::default()
        };
        let options = bus_pseudoalign_options(&args, None);

        assert!(options.dfk_onlist);
    }

    #[test]
    fn bus_aa_enables_dfk_onlist_by_default() {
        let args = BusArgs {
            aa: true,
            ..BusArgs::default()
        };
        let options = bus_pseudoalign_options(&args, None);

        assert!(options.dfk_onlist);
    }

    #[test]
    fn parse_v3_lengths_and_strand_match_preset() {
        let spec = technology_spec("PARSEV3").unwrap();
        assert_eq!(spec.nfiles, 2);
        assert_eq!(total_len(&spec.bc), Some(24));
        assert_eq!(total_len(&spec.umi), Some(10));
        assert_eq!(spec.default_strand, Some(StrandSpecific::FR));
    }

    #[test]
    fn compound_preset_slices_match_kallisto_layouts() {
        let records = vec![
            FastqRecord {
                header: b"@r1".to_vec(),
                seq: b"0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz".to_vec(),
                plus: b"+".to_vec(),
                qual: vec![b'I'; 62],
            },
            FastqRecord {
                header: b"@r2".to_vec(),
                seq: b"abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789".to_vec(),
                plus: b"+".to_vec(),
                qual: vec![b'I'; 62],
            },
            FastqRecord {
                header: b"@r3".to_vec(),
                seq: b"QRST".to_vec(),
                plus: b"+".to_vec(),
                qual: vec![b'I'; 4],
            },
        ];

        let parse = technology_spec("PARSEV3").unwrap();
        assert_eq!(
            concat_slices(&records, &parse.bc, false).unwrap(),
            b"klmnopqrEFGHIJKLYZ012345"
        );
        assert_eq!(
            concat_slices(&records, &parse.umi, false).unwrap(),
            b"abcdefghij"
        );
        assert_eq!(
            concat_sequence_slices(&records, &parse, false).unwrap(),
            records[0].seq
        );

        let petri = technology_spec("PETRISEQ").unwrap();
        assert_eq!(
            concat_slices(&records, &petri.bc, false).unwrap(),
            b"789ABCDTUVWXYZopqrstuv"
        );
        assert_eq!(
            concat_slices(&records, &petri.umi, false).unwrap(),
            b"0123456"
        );
        assert_eq!(
            concat_sequence_slices(&records, &petri, false).unwrap(),
            b"abcdefghijklmnopq"
        );

        let surecell = technology_spec("SURECELL").unwrap();
        assert_eq!(
            concat_slices(&records, &surecell.bc, false).unwrap(),
            b"012345LMNOPQghijkl"
        );
        assert_eq!(
            concat_slices(&records, &surecell.umi, false).unwrap(),
            b"pqrstuvw"
        );

        let bdwta = technology_spec("BDWTA").unwrap();
        assert_eq!(
            concat_slices(&records, &bdwta.bc, false).unwrap(),
            b"012345678LMNOPQRSThijklmnop"
        );
        assert_eq!(
            concat_slices(&records, &bdwta.umi, false).unwrap(),
            b"qrstuvwx"
        );

        let vasa = technology_spec("VASA-SEQ").unwrap();
        assert_eq!(
            concat_slices(&records, &vasa.bc, false).unwrap(),
            b"6789ABCD"
        );
        assert_eq!(
            concat_slices(&records, &vasa.umi, false).unwrap(),
            b"012345"
        );
        assert_eq!(
            concat_sequence_slices(&records, &vasa, false).unwrap(),
            b"EFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"
        );

        let dropseq = technology_spec("DROPSEQ").unwrap();
        assert_eq!(
            concat_slices(&records, &dropseq.bc, false).unwrap(),
            b"0123456789AB"
        );
        assert_eq!(
            concat_slices(&records, &dropseq.umi, false).unwrap(),
            b"CDEFGHIJ"
        );

        let indrops_v1 = technology_spec("INDROPSV1").unwrap();
        assert_eq!(
            concat_slices(&records, &indrops_v1.bc, false).unwrap(),
            b"0123456789AUVWXYZab"
        );
        assert_eq!(
            concat_slices(&records, &indrops_v1.umi, false).unwrap(),
            b"ghijkl"
        );

        let indrops_v2 = technology_spec("INDROPSV2").unwrap();
        assert_eq!(
            concat_slices(&records, &indrops_v2.bc, false).unwrap(),
            b"abcdefghijkEFGHIJKL"
        );
        assert_eq!(
            concat_slices(&records, &indrops_v2.umi, false).unwrap(),
            b"QRSTUV"
        );

        let indrops_v3 = technology_spec("INDROPSV3").unwrap();
        assert_eq!(
            concat_slices(&records, &indrops_v3.bc, false).unwrap(),
            b"01234567abcdefgh"
        );
        assert_eq!(
            concat_slices(&records, &indrops_v3.umi, false).unwrap(),
            b"ijklmn"
        );
        assert_eq!(
            concat_sequence_slices(&records, &indrops_v3, false).unwrap(),
            b"QRST"
        );

        let celseq = technology_spec("CELSEQ").unwrap();
        assert_eq!(
            concat_slices(&records, &celseq.bc, false).unwrap(),
            b"01234567"
        );
        assert_eq!(
            concat_slices(&records, &celseq.umi, false).unwrap(),
            b"89AB"
        );

        let celseq2 = technology_spec("CELSEQ2").unwrap();
        assert_eq!(
            concat_slices(&records, &celseq2.bc, false).unwrap(),
            b"6789AB"
        );
        assert_eq!(
            concat_slices(&records, &celseq2.umi, false).unwrap(),
            b"012345"
        );

        let scrbseq = technology_spec("SCRBSEQ").unwrap();
        assert_eq!(
            concat_slices(&records, &scrbseq.bc, false).unwrap(),
            b"012345"
        );
        assert_eq!(
            concat_slices(&records, &scrbseq.umi, false).unwrap(),
            b"6789ABCDEF"
        );
    }

    #[test]
    fn smartseq3_uses_full_length_barcodes() {
        let spec = technology_spec("SMARTSEQ3").unwrap();
        assert_eq!(spec.nfiles, 4);
        assert_eq!(total_len(&spec.bc), None);
        assert_eq!(total_len(&spec.umi), Some(19));
        assert!(spec.paired);
    }

    #[test]
    fn custom_technology_parser_accepts_kallisto_triples() {
        let spec = technology_spec("0,0,16:0,16,28:1,0,0").unwrap();
        assert_eq!(spec.nfiles, 2);
        assert_eq!(total_len(&spec.bc), Some(16));
        assert_eq!(total_len(&spec.umi), Some(12));
        assert!(!spec.paired);
    }

    #[test]
    fn custom_technology_parser_keeps_at_least_one_input_file() {
        let spec = technology_spec("-1,-1,-1:-1,-1,-1:-1,-1,-1").unwrap();
        assert_eq!(spec.nfiles, 1);
        assert_eq!(total_len(&spec.bc), None);
        assert_eq!(total_len(&spec.umi), None);
    }

    #[test]
    fn custom_rx_technology_reads_umi_from_fastq_comments() {
        let spec = technology_spec("0,0,16:RX:1,0,0").unwrap();
        assert!(spec.keep_fastq_comments);
        assert!(!spec.bulk_like());
        assert_eq!(total_len(&spec.umi), None);
    }

    #[test]
    fn custom_technology_parser_reports_colon_count_and_empty_lists() {
        let no_colons = technology_spec("0,0,16").unwrap_err().to_string();
        assert!(no_colons.contains("none found"), "{no_colons}");

        let one_colon = technology_spec("0,0,16:0,16,28").unwrap_err().to_string();
        assert!(one_colon.contains("only one found"), "{one_colon}");

        let three_colons = technology_spec("0,0,16:0,16,28:1,0,0:extra")
            .unwrap_err()
            .to_string();
        assert!(three_colons.contains("three found"), "{three_colons}");

        let empty_barcode = technology_spec(":0,16,28:1,0,0").unwrap_err().to_string();
        assert!(
            empty_barcode.contains("custom technology barcode list is empty"),
            "{empty_barcode}"
        );
    }

    #[test]
    fn matrix_ec_uses_comma_separated_transcript_lists() {
        let dir = tempfile::tempdir().expect("tempdir");
        let path = dir.path().join("matrix.ec");

        write_matrix_ec(&path, &[vec![0, 2], vec![1]]).expect("write matrix.ec");

        assert_eq!(std::fs::read_to_string(path).unwrap(), "0\t0,2\n1\t1\n");
    }

    #[test]
    fn transcripts_sidecar_uses_only_onlist_targets() {
        let dir = tempfile::tempdir().expect("tempdir");
        let path = dir.path().join("transcripts.txt");
        let names = vec![
            "tx0".to_string(),
            "tx1".to_string(),
            "d_list.0".to_string(),
            "tx2".to_string(),
        ];
        let onlist = vec![true, true, false, true];

        write_transcripts(&path, &names, Some(&onlist)).expect("write transcripts");

        assert_eq!(std::fs::read_to_string(path).unwrap(), "tx0\ntx1\ntx2\n");
        assert_eq!(onlist_target_count(&names, Some(&onlist)), 3);
        assert_eq!(onlist_target_count(&names, None), 4);
    }

    #[test]
    fn all_listed_technology_names_are_accepted() {
        for name in [
            "10xv1",
            "10xv2",
            "10xv3",
            "10xv4",
            "Bulk",
            "ParseV3",
            "SmartSeq2",
            "BDWTA",
            "CELSeq",
            "CELSeq2",
            "DropSeq",
            "inDropsv1",
            "inDropsv2",
            "inDropsv3",
            "MATQSEQ",
            "PETRISEQ",
            "SCRBSeq",
            "SmartSeq3",
            "SPLiT-seq",
            "STORM-seq",
            "SureCell",
            "VASA-seq",
            "Visium",
        ] {
            let spec = technology_spec(name).unwrap_or_else(|err| {
                panic!("listed technology {name} should parse: {err}");
            });
            assert!(spec.nfiles > 0, "listed technology {name} should use reads");
            assert!(
                !spec.seq.is_empty(),
                "listed technology {name} should define cDNA slices"
            );
            assert!(
                !spec.bc.is_empty(),
                "listed technology {name} should define barcode slices"
            );
            assert!(
                !spec.umi.is_empty(),
                "listed technology {name} should define UMI slices"
            );
        }
    }

    #[test]
    fn custom_two_sequence_slices_are_paired_only_when_requested() {
        let unpaired = technology_spec("-1,-1,-1:0,0,4:0,4,8,1,0,4").unwrap();
        assert!(!unpaired.paired);
        assert_eq!(unpaired.nfiles, 2);
        assert_eq!(unpaired.seq.len(), 2);

        let paired = technology_spec("-1,-1,-1:0,0,4:0,4,8,1,0,4%NONE%PAIRED").unwrap();
        assert!(paired.paired);
        assert_eq!(paired.nfiles, 2);
        assert_eq!(paired.seq.len(), 2);
    }

    #[test]
    fn technology_suffix_applies_strand_and_paired_mode() {
        let forward = technology_spec("0,0,16:0,16,28:1,0,0%FORWARD").unwrap();
        assert_eq!(forward.default_strand, Some(StrandSpecific::FR));
        assert!(!forward.paired);

        let none_paired = technology_spec("0,0,16:0,16,28:1,0,0%NONE%PAIRED").unwrap();
        assert_eq!(none_paired.default_strand, None);
        assert!(none_paired.paired);

        let reverse_paired = technology_spec("0,0,16:0,16,28:1,0,0%REVERSE%PAIRED").unwrap();
        assert_eq!(reverse_paired.default_strand, Some(StrandSpecific::RF));
        assert!(reverse_paired.paired);
        assert_eq!(reverse_paired.nfiles, 3);
        assert_eq!(reverse_paired.seq.len(), 2);
        assert_eq!(reverse_paired.seq[1].file, Some(2));

        let invalid_strand = technology_spec("0,0,16:0,16,28:1,0,0%SIDEWAYS")
            .unwrap_err()
            .to_string();
        assert!(
            invalid_strand.contains("invalid technology strand suffix: SIDEWAYS"),
            "{invalid_strand}"
        );
        let invalid_pairing = technology_spec("0,0,16:0,16,28:1,0,0%FORWARD%MATED")
            .unwrap_err()
            .to_string();
        assert!(
            invalid_pairing.contains("invalid technology pairing suffix: MATED"),
            "{invalid_pairing}"
        );
        let extra = technology_spec("0,0,16:0,16,28:1,0,0%FORWARD%PAIRED%EXTRA")
            .unwrap_err()
            .to_string();
        assert!(
            extra.contains("unexpected technology suffix field: EXTRA"),
            "{extra}"
        );
    }

    #[test]
    fn multiple_sequence_slices_are_separated_by_n() {
        let records = vec![
            FastqRecord {
                header: b"@r1".to_vec(),
                seq: b"AAAACCCC".to_vec(),
                plus: b"+".to_vec(),
                qual: b"IIIIIIII".to_vec(),
            },
            FastqRecord {
                header: b"@r2".to_vec(),
                seq: b"GGGGTTTT".to_vec(),
                plus: b"+".to_vec(),
                qual: b"IIIIIIII".to_vec(),
            },
        ];
        let spec = technology_spec("-1,-1,-1:0,0,4:0,4,8,1,0,4").unwrap();
        assert_eq!(
            concat_sequence_slices(&records, &spec, false).unwrap(),
            b"CCCCNGGGG"
        );
        assert_eq!(
            bus_query_quality(&records, &spec, false).unwrap(),
            b"IIII!IIII"
        );
    }

    #[test]
    fn batch_paired_flens_writer_preserves_each_batch_histogram() {
        let dir = tempfile::tempdir().expect("tempdir");
        let path = dir.path().join("flens.txt");
        let mut first = vec![0u32; 8];
        let mut second = vec![0u32; 8];
        first[3] = 2;
        second[5] = 7;

        write_batch_paired_flens(&path, &[first, second]).expect("write flens");

        assert_eq!(
            std::fs::read_to_string(path).unwrap(),
            "0 0 0 2 0 0 0 0\n0 0 0 0 0 7 0 0\n"
        );
    }

    #[test]
    fn batch_long_flens_writer_preserves_each_batch_estimate() {
        let dir = tempfile::tempdir().expect("tempdir");
        let path = dir.path().join("flens.txt");
        let transcript_lengths = vec![100, 120];
        let sums = vec![vec![80, 0], vec![0, 90]];
        let counts = vec![vec![2, 0], vec![0, 3]];

        write_batch_long_flens(&path, &transcript_lengths, &sums, &counts, 31)
            .expect("write long flens");

        assert_eq!(std::fs::read_to_string(path).unwrap(), "9 89\n69 1\n");
    }
}
