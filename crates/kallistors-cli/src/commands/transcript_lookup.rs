use anyhow::{Result, anyhow};
use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

fn parse_ids_from_list(path: &Path) -> Result<HashSet<u32>> {
    let file = File::open(path)?;
    let mut ids = HashSet::new();
    for line in BufReader::new(file).lines() {
        let line = line?;
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        let value: u32 = trimmed
            .parse()
            .map_err(|_| anyhow!("invalid transcript id: {trimmed}"))?;
        ids.insert(value);
    }
    Ok(ids)
}

fn parse_ids_from_hits(path: &Path) -> Result<HashSet<u32>> {
    let file = File::open(path)?;
    let mut reader = BufReader::new(file);
    let mut header = String::new();
    if reader.read_line(&mut header)? == 0 {
        return Ok(HashSet::new());
    }
    let header = header.trim_end();
    let cols: Vec<&str> = header.split('\t').collect();
    let ec_idx = cols
        .iter()
        .position(|&c| c == "ec")
        .ok_or_else(|| anyhow!("hits file missing 'ec' column"))?;
    let mut ids = HashSet::new();
    let mut line = String::new();
    while reader.read_line(&mut line)? > 0 {
        let row = line.trim_end();
        if row.is_empty() {
            line.clear();
            continue;
        }
        let parts: Vec<&str> = row.split('\t').collect();
        if let Some(ec_field) = parts.get(ec_idx)
            && *ec_field != "-"
            && !ec_field.is_empty()
        {
            for val in ec_field.split(',') {
                if val.is_empty() {
                    continue;
                }
                let id: u32 = val
                    .parse()
                    .map_err(|_| anyhow!("invalid transcript id in hits file: {val}"))?;
                ids.insert(id);
            }
        }
        line.clear();
    }
    Ok(ids)
}

pub fn run(
    index: &Path,
    ids: Option<&Path>,
    hits: Option<&Path>,
    out: Option<&Path>,
) -> Result<()> {
    if ids.is_none() && hits.is_none() {
        return Err(anyhow!("--ids or --hits is required"));
    }
    let mut all_ids = HashSet::new();
    if let Some(path) = ids {
        all_ids.extend(parse_ids_from_list(path)?);
    }
    if let Some(path) = hits {
        all_ids.extend(parse_ids_from_hits(path)?);
    }

    let mut id_list: Vec<u32> = all_ids.into_iter().collect();
    id_list.sort_unstable();

    let index = kallistors::index::Index::load(index)
        .map_err(|err| anyhow!("transcript-lookup failed: {err}"))?;

    let writer: Box<dyn Write> = match out {
        Some(path) => Box::new(BufWriter::new(File::create(path)?)),
        None => Box::new(BufWriter::new(std::io::stdout())),
    };
    let mut writer = writer;
    writeln!(writer, "id\tname\tlength")?;
    for id in id_list {
        if let Some(info) = index.transcripts.get(id as usize) {
            writeln!(writer, "{}\t{}\t{}", id, info.name, info.length)?;
        } else {
            writeln!(writer, "{}\t-\t0", id)?;
        }
    }
    Ok(())
}
