use anyhow::{Result, anyhow, bail};
use kallistors::io::ReadSource;
use std::collections::HashMap;
use std::fs;
use std::path::{Path, PathBuf};

#[derive(Debug, Clone)]
pub(crate) struct ReadGroup {
    pub(crate) id: Option<String>,
    pub(crate) files: Vec<PathBuf>,
    pub(crate) batch_index: usize,
    pub(crate) interleaved: bool,
}

pub(crate) fn read_groups(args: &super::BusArgs, nfiles: usize) -> Result<Vec<ReadGroup>> {
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

pub(crate) fn validate_read_files(args: &super::BusArgs, groups: &[ReadGroup]) -> Result<()> {
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

pub(crate) fn report_verbose_read_groups(groups: &[ReadGroup]) {
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

pub(crate) fn next_parallel_batch<R: ReadSource>(
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

pub(crate) fn next_interleaved_batch<R: ReadSource>(
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

pub(crate) fn should_infer_bulk_paired(args: &super::BusArgs) -> Result<bool> {
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
