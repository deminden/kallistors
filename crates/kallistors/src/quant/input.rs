//! EC-count and transcript-sequence input parsing.

use super::*;

pub fn load_ec_counts(path: &Path) -> Result<EcCountsInput> {
    let file = File::open(path).map_err(|_| Error::MissingFile(path.to_path_buf()))?;
    let reader = BufReader::new(file);
    let mut classes: Vec<Vec<u32>> = Vec::new();
    let mut counts: Vec<u32> = Vec::new();

    for (line_no, line) in reader.lines().enumerate() {
        let line = line?;
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }
        let mut parts = line.splitn(3, '\t');
        let ec_id: usize = parts
            .next()
            .ok_or_else(|| Error::InvalidFormat(format!("missing EC id at line {}", line_no + 1)))?
            .parse()
            .map_err(|_| Error::InvalidFormat(format!("invalid EC id at line {}", line_no + 1)))?;
        let count: u32 = parts
            .next()
            .ok_or_else(|| {
                Error::InvalidFormat(format!("missing EC count at line {}", line_no + 1))
            })?
            .parse()
            .map_err(|_| {
                Error::InvalidFormat(format!("invalid EC count at line {}", line_no + 1))
            })?;
        let list_str = parts
            .next()
            .ok_or_else(|| {
                Error::InvalidFormat(format!("missing EC list at line {}", line_no + 1))
            })?
            .trim();
        let mut ids = Vec::new();
        if !list_str.is_empty() {
            for token in list_str.split(',') {
                if token.is_empty() {
                    continue;
                }
                let val: u32 = token.parse().map_err(|_| {
                    Error::InvalidFormat(format!("invalid transcript id at line {}", line_no + 1))
                })?;
                ids.push(val);
            }
        }
        if ec_id != classes.len() {
            return Err(Error::InvalidFormat(format!(
                "EC id out of order at line {}",
                line_no + 1
            )));
        }
        classes.push(ids);
        counts.push(count);
    }

    Ok(EcCountsInput {
        ec_list: EcList { classes },
        counts,
    })
}

pub fn load_transcript_sequences(fasta: &Path, names: &[String]) -> Result<Vec<Vec<u8>>> {
    let file = File::open(fasta).map_err(|_| Error::MissingFile(fasta.to_path_buf()))?;
    let reader = BufReader::new(file);
    let mut seqs: HashMap<String, Vec<u8>> = HashMap::new();
    let mut current_name: Option<String> = None;
    let mut current_seq: Vec<u8> = Vec::new();
    for line in reader.lines() {
        let line = line?;
        if let Some(rest) = line.strip_prefix('>') {
            if let Some(name) = current_name.take() {
                seqs.insert(name, std::mem::take(&mut current_seq));
            }
            current_name = Some(rest.trim().to_string());
        } else {
            current_seq.extend_from_slice(line.trim().as_bytes());
        }
    }
    if let Some(name) = current_name.take() {
        seqs.insert(name, current_seq);
    }

    let mut out = Vec::with_capacity(names.len());
    for name in names {
        let seq = seqs
            .get(name)
            .cloned()
            .ok_or_else(|| Error::InvalidFormat(format!("missing transcript {}", name)))?;
        out.push(seq);
    }
    Ok(out)
}
