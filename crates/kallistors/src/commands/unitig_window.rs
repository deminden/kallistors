use anyhow::{Result, anyhow};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

fn parse_minimizer_info(path: &Path) -> Result<HashMap<(String, usize), (usize, String)>> {
    let file = File::open(path)?;
    let mut reader = BufReader::new(file);
    let mut header = String::new();
    if reader.read_line(&mut header)? == 0 {
        return Ok(HashMap::new());
    }
    let header = header.trim_end();
    let cols: Vec<&str> = header.split('\t').collect();
    let read_idx = cols
        .iter()
        .position(|&c| c == "read")
        .ok_or_else(|| anyhow!("minimizer file missing 'read' column"))?;
    let pos_idx = cols
        .iter()
        .position(|&c| c == "read_pos")
        .ok_or_else(|| anyhow!("minimizer file missing 'read_pos' column"))?;
    let min_pos_idx = cols
        .iter()
        .position(|&c| c == "min_pos")
        .ok_or_else(|| anyhow!("minimizer file missing 'min_pos' column"))?;
    let kmer_idx = cols
        .iter()
        .position(|&c| c == "kmer")
        .ok_or_else(|| anyhow!("minimizer file missing 'kmer' column"))?;

    let mut map = HashMap::new();
    let mut line = String::new();
    while reader.read_line(&mut line)? > 0 {
        let row = line.trim_end();
        if row.is_empty() {
            line.clear();
            continue;
        }
        let parts: Vec<&str> = row.split('\t').collect();
        let read = parts.get(read_idx).unwrap_or(&"").trim();
        let pos = parts.get(pos_idx).unwrap_or(&"").trim();
        let min_pos = parts.get(min_pos_idx).unwrap_or(&"").trim();
        let kmer = parts.get(kmer_idx).unwrap_or(&"").trim();
        if read.is_empty() || pos.is_empty() || min_pos.is_empty() {
            line.clear();
            continue;
        }
        let pos: usize = match pos.parse() {
            Ok(v) => v,
            Err(_) => {
                line.clear();
                continue;
            }
        };
        let min_pos: usize = match min_pos.parse() {
            Ok(v) => v,
            Err(_) => {
                line.clear();
                continue;
            }
        };
        map.entry((read.to_string(), pos))
            .or_insert((min_pos, kmer.to_string()));
        line.clear();
    }
    Ok(map)
}

pub fn run(
    index_path: &Path,
    positions_path: &Path,
    minimizer_path: &Path,
    out: &Path,
) -> Result<()> {
    let index = kallistors::pseudoalign::build_bifrost_index(index_path)?;
    let diff = index.k.saturating_sub(index.g) as isize;
    let read_info = parse_minimizer_info(minimizer_path)?;

    let file = File::open(positions_path)?;
    let mut reader = BufReader::new(file);
    let mut header = String::new();
    if reader.read_line(&mut header)? == 0 {
        return Ok(());
    }
    let header = header.trim_end();
    let cols: Vec<&str> = header.split('\t').collect();
    let read_idx = cols
        .iter()
        .position(|&c| c == "read")
        .ok_or_else(|| anyhow!("positions file missing 'read' column"))?;
    let pos_idx = cols
        .iter()
        .position(|&c| c == "read_pos")
        .ok_or_else(|| anyhow!("positions file missing 'read_pos' column"))?;
    let unitig_idx = cols
        .iter()
        .position(|&c| c == "unitig_id")
        .ok_or_else(|| anyhow!("positions file missing 'unitig_id' column"))?;
    let rel_idx = cols
        .iter()
        .position(|&c| c == "rel_pos")
        .ok_or_else(|| anyhow!("positions file missing 'rel_pos' column"))?;
    let kind_idx = cols
        .iter()
        .position(|&c| c == "kind")
        .ok_or_else(|| anyhow!("positions file missing 'kind' column"))?;
    let note_idx = cols
        .iter()
        .position(|&c| c == "note")
        .ok_or_else(|| anyhow!("positions file missing 'note' column"))?;

    let mut writer = BufWriter::new(File::create(out)?);
    writeln!(
        writer,
        "read\tread_pos\tunitig_id\tkind\trel_pos\tmin_pos\tstart\tend\tunitig_len\tin_bounds\tunitig_kmer\tread_kmer\tnote\tscan_start\tscan_end\tscan_matches\tscan_match_count"
    )?;

    let mut line = String::new();
    while reader.read_line(&mut line)? > 0 {
        let row = line.trim_end();
        if row.is_empty() {
            line.clear();
            continue;
        }
        let parts: Vec<&str> = row.split('\t').collect();
        let read = parts.get(read_idx).unwrap_or(&"").trim().to_string();
        let pos: usize = match parts.get(pos_idx).unwrap_or(&"").trim().parse() {
            Ok(v) => v,
            Err(_) => {
                line.clear();
                continue;
            }
        };
        let unitig_id: usize = match parts.get(unitig_idx).unwrap_or(&"").trim().parse() {
            Ok(v) => v,
            Err(_) => {
                line.clear();
                continue;
            }
        };
        let rel_pos: isize = match parts.get(rel_idx).unwrap_or(&"").trim().parse() {
            Ok(v) => v,
            Err(_) => {
                line.clear();
                continue;
            }
        };
        let kind = parts.get(kind_idx).unwrap_or(&"").trim();
        let note = parts.get(note_idx).unwrap_or(&"").trim();

        let Some((min_pos, read_kmer)) = read_info.get(&(read.clone(), pos)).cloned() else {
            line.clear();
            continue;
        };
        let start = rel_pos - min_pos as isize;
        let end = start + index.k as isize;

        let seq: &[u8] = if kind == "k" {
            index
                .km_unitigs
                .get(unitig_id)
                .map(|v| v.as_slice())
                .unwrap_or(&[])
        } else {
            index
                .unitigs
                .get(unitig_id)
                .map(|v| v.as_slice())
                .unwrap_or(&[])
        };
        let unitig_len = seq.len();
        let in_bounds = start >= 0 && end >= 0 && (end as usize) <= unitig_len;
        let unitig_kmer = if in_bounds {
            let s = start as usize;
            String::from_utf8_lossy(&seq[s..s + index.k]).into_owned()
        } else {
            "-".to_string()
        };
        let scan_start = rel_pos.saturating_sub(diff).max(0) as usize;
        let scan_end = (rel_pos + diff).min(unitig_len as isize) as usize;
        let mut scan_matches = Vec::new();
        if unitig_len >= index.k && scan_start < unitig_len {
            let max_start = scan_end.saturating_sub(index.k);
            let scan_end_pos = max_start.min(unitig_len.saturating_sub(index.k));
            for s in scan_start..=scan_end_pos {
                if &seq[s..s + index.k] == read_kmer.as_bytes() {
                    scan_matches.push(s);
                }
            }
        }
        let scan_match_text = if scan_matches.is_empty() {
            "-".to_string()
        } else {
            scan_matches
                .iter()
                .map(|v| v.to_string())
                .collect::<Vec<_>>()
                .join(",")
        };
        let scan_match_count = scan_matches.len();

        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            read,
            pos,
            unitig_id,
            kind,
            rel_pos,
            min_pos,
            start,
            end,
            unitig_len,
            if in_bounds { "1" } else { "0" },
            unitig_kmer,
            read_kmer,
            note,
            scan_start,
            scan_end,
            scan_match_text,
            scan_match_count
        )?;
        line.clear();
    }
    Ok(())
}
