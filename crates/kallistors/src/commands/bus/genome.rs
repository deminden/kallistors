//! GTF parsing and transcript-to-genome projection.

use anyhow::{Result, anyhow, bail};
use noodles_sam::alignment::{
    record::cigar::{Op, op::Kind},
    record_buf::Cigar,
};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

#[derive(Clone)]
pub(super) struct GenomeModel {
    pub(super) chromosomes: Vec<(String, usize)>,
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
pub(super) struct GenomeAlignment {
    pub(super) reference_sequence_id: usize,
    pub(super) start: usize,
    pub(super) cigar: Cigar,
    pub(super) reverse: bool,
}

struct GenomeSegment {
    start: usize,
    end: usize,
}

impl GenomeModel {
    pub(super) fn from_gtf(
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
            // GTF intervals are one-based and closed. Internally retain
            // zero-based half-open intervals, so the numeric end is unchanged.
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

    pub(super) fn project(
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
        // Exons were walked in transcript order, which runs backward on minus
        // strands. BAM CIGAR operations always follow increasing genome positions.
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

pub(super) fn write_gene_list(path: &Path, model: &GenomeModel) -> Result<()> {
    let mut writer = BufWriter::new(File::create(path)?);
    for (idx, gene) in model.genes.iter().enumerate() {
        writeln!(writer, "{idx}\t{}\t{}", gene.id, gene.name)?;
    }
    Ok(())
}
