use anyhow::{Result, anyhow, bail};
use kallistors::pseudoalign::StrandSpecific;

#[derive(Debug, Clone, Copy)]
pub(crate) struct SliceSpec {
    pub(crate) file: Option<usize>,
    pub(crate) start: usize,
    pub(crate) stop: Option<usize>,
}

#[derive(Debug, Clone)]
pub(crate) struct TechnologySpec {
    pub(crate) nfiles: usize,
    pub(crate) seq: Vec<SliceSpec>,
    pub(crate) bc: Vec<SliceSpec>,
    pub(crate) umi: Vec<SliceSpec>,
    pub(crate) paired: bool,
    pub(crate) keep_fastq_comments: bool,
    pub(crate) default_strand: Option<StrandSpecific>,
    pub(crate) tag_len: Option<usize>,
}

impl TechnologySpec {
    pub(crate) fn bulk_like(&self) -> bool {
        !self.keep_fastq_comments && self.umi.first().is_none_or(|spec| spec.file.is_none())
    }

    fn with_default_strand(mut self, strand: StrandSpecific) -> Self {
        self.default_strand = Some(strand);
        self
    }
}

pub(crate) fn technology_spec(name: &str) -> Result<TechnologySpec> {
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
        .with_default_strand(StrandSpecific::FR),
        "10XV2" => tech(
            2,
            vec![slice(1, 0, None)],
            vec![slice(0, 16, Some(26))],
            vec![slice(0, 0, Some(16))],
            false,
        )
        .with_default_strand(StrandSpecific::FR),
        "10XV3" | "10XV4" | "VISIUM" => tech(
            2,
            vec![slice(1, 0, None)],
            vec![slice(0, 16, Some(28))],
            vec![slice(0, 0, Some(16))],
            false,
        )
        .with_default_strand(StrandSpecific::FR),
        "SMARTSEQ3" => tech(
            4,
            vec![slice(2, 22, None), slice(3, 0, None)],
            vec![slice(2, 0, Some(19))],
            vec![slice(0, 0, None), slice(1, 0, None)],
            true,
        )
        .with_default_strand(StrandSpecific::FR),
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
        .with_default_strand(StrandSpecific::FR),
        "CELSEQ2" => tech(
            2,
            vec![slice(1, 0, None)],
            vec![slice(0, 0, Some(6))],
            vec![slice(0, 6, Some(12))],
            false,
        )
        .with_default_strand(StrandSpecific::FR),
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
        .with_default_strand(StrandSpecific::FR),
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
        .with_default_strand(StrandSpecific::FR),
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
        .with_default_strand(StrandSpecific::FR),
        "STORM-SEQ" => tech(
            2,
            vec![slice(0, 0, None), slice(1, 14, None)],
            vec![slice(1, 0, Some(8))],
            vec![sentinel_slice()],
            true,
        )
        .with_default_strand(StrandSpecific::RF),
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
        .with_default_strand(StrandSpecific::FR),
        "VASA-SEQ" => tech(
            1,
            vec![slice(0, 14, None)],
            vec![slice(0, 0, Some(6))],
            vec![slice(0, 6, Some(14))],
            false,
        )
        .with_default_strand(StrandSpecific::FR),
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
            spec.default_strand = Some(StrandSpecific::FR);
        } else if strand.starts_with("REVERSE") {
            spec.default_strand = Some(StrandSpecific::RF);
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

pub(crate) fn mark_technology_paired(spec: &mut TechnologySpec, source: &str) -> Result<()> {
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

pub(crate) fn add_unpaired_mate_slice(spec: &mut TechnologySpec, source: &str) -> Result<()> {
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

pub(crate) fn technology_base_upper(name: &str) -> String {
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
            if !((start == -1 && stop == -1) || (start == 0 && stop == 0)) {
                bail!("sentinel technology slices must be -1,-1,-1 or -1,0,0");
            }
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

pub(crate) fn total_len(specs: &[SliceSpec]) -> Option<u32> {
    let mut total = 0u32;
    for spec in specs {
        spec.file?;
        let stop = spec.stop?;
        total = total.saturating_add(u32::try_from(stop.saturating_sub(spec.start)).ok()?);
    }
    Some(total)
}

#[cfg(test)]
pub(crate) fn concat_slices(
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

pub(crate) fn concat_slices_optional(
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

pub(crate) fn concat_sequence_slices(
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

pub(crate) fn concat_quality_slices(
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

pub(crate) fn extract_sequence_slice(
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

pub(crate) fn extract_quality_slice(
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

pub(crate) fn extract_slice(
    records: &[kallistors::io::FastqRecord],
    spec: SliceSpec,
) -> Result<Vec<u8>> {
    extract_slice_from(records, spec, |record| record.seq.as_slice(), "read")
}

pub(crate) fn extract_slice_optional(
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
