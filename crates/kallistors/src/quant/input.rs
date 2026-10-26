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
    let mut reader = BufReader::new(file);
    enum SequenceEntry {
        Missing,
        Owned(Vec<u8>),
        Emitted(usize),
    }
    // Header keys borrow stable requested names. Each unique request owns one
    // slot, while the request order records duplicates without cloning strings.
    let mut wanted: HashMap<&str, usize> = HashMap::with_capacity(names.len());
    let mut entries = Vec::with_capacity(names.len());
    let mut requests = Vec::with_capacity(names.len());
    for name in names {
        let slot = match wanted.entry(name.as_str()) {
            std::collections::hash_map::Entry::Occupied(entry) => *entry.get(),
            std::collections::hash_map::Entry::Vacant(entry) => {
                let slot = entries.len();
                entries.push(SequenceEntry::Missing);
                entry.insert(slot);
                slot
            }
        };
        requests.push(slot);
    }
    let mut current_slot = None;
    let mut before_first_header = true;
    let mut current_seq = Vec::new();
    let mut line = String::new();
    loop {
        if line.capacity() > 64 * 1024 {
            line = String::new();
        } else {
            line.clear();
        }
        if reader.read_line(&mut line)? == 0 {
            break;
        }
        if let Some(rest) = line.strip_prefix('>') {
            if let Some(slot) = current_slot.take() {
                entries[slot] = SequenceEntry::Owned(std::mem::take(&mut current_seq));
            }
            current_slot = wanted.get(rest.trim()).copied();
            before_first_header = false;
            if current_slot.is_none() {
                current_seq.clear();
            }
        } else if before_first_header || current_slot.is_some() {
            // Retain the parser's historical pre-header sequence prefix. Even
            // skipped records are read as UTF-8 above, including empty requests.
            current_seq.extend_from_slice(line.trim().as_bytes());
        }
    }
    if let Some(slot) = current_slot {
        entries[slot] = SequenceEntry::Owned(current_seq);
    }
    let mut out: Vec<Vec<u8>> = Vec::with_capacity(names.len());
    for (&slot, name) in requests.iter().zip(names) {
        let entry = &mut entries[slot];
        let sequence = match entry {
            SequenceEntry::Missing => {
                return Err(Error::InvalidFormat(format!("missing transcript {}", name)));
            }
            SequenceEntry::Owned(bytes) => {
                let bytes = std::mem::take(bytes);
                *entry = SequenceEntry::Emitted(out.len());
                bytes
            }
            SequenceEntry::Emitted(index) => out[*index].clone(),
        };
        out.push(sequence);
    }
    Ok(out)
}

#[cfg(test)]
mod sequence_input_tests {
    use super::*;

    #[test]
    fn transcript_buffers_preserve_requested_duplicates_last_header_and_input_errors() {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("transcripts.fa");
        // Last duplicate FASTA header wins; requested order/duplicates,
        // case and trimming across CRLF/unterminated lines are observable.
        std::fs::write(&path, b">a\r\nACGT\r\n>b\r\n ttaa \r\n>a \r\n GG \r\nCC").unwrap();
        let names = vec!["a".to_owned(), "b".to_owned(), "a".to_owned()];
        let mut sequences = load_transcript_sequences(&path, &names).unwrap();
        assert_eq!(
            sequences,
            [b"GGCC".to_vec(), b"ttaa".to_vec(), b"GGCC".to_vec()]
        );
        sequences[0][0] = b'T';
        assert_eq!(sequences[2].as_slice(), b"GGCC");
        assert!(
            matches!(load_transcript_sequences(&path, &["b".into(), "missing".into()]),
            Err(Error::InvalidFormat(message)) if message == "missing transcript missing")
        );

        // Unrequested records, empty requested records, pre-header prefixes and
        // complete header whitespace retain the original parser contract.
        std::fs::write(
            &path,
            b"prefix\n> a \nAC\n>ignored\nGGGG\n>b\n>a\nT\n>unused\nCC\n",
        )
        .unwrap();
        assert_eq!(
            load_transcript_sequences(&path, &["b".into(), "a".into(), "a".into()]).unwrap(),
            [Vec::<u8>::new(), b"T".to_vec(), b"T".to_vec()]
        );
        std::fs::write(&path, b" pre \n> a \nFIX\n>b\n").unwrap();
        assert_eq!(
            load_transcript_sequences(&path, &["a".into(), "b".into()]).unwrap(),
            [b"preFIX".to_vec(), Vec::<u8>::new()]
        );
        std::fs::write(&path, b">ignored\nAC\xffGT\n>a\nT\n").unwrap();
        assert!(matches!(load_transcript_sequences(&path,&["a".into()]),
            Err(Error::Io(error)) if error.kind()==std::io::ErrorKind::InvalidData));

        // A long record followed by a small record must retain both contents.
        let long = "A".repeat(128 * 1024);
        std::fs::write(&path, format!(">a\n{long}\n>b\nG\n")).unwrap();
        let loaded = load_transcript_sequences(&path, &["a".into(), "b".into()]).unwrap();
        assert_eq!(loaded[0].as_slice(), long.as_bytes());
        assert_eq!(loaded[1].as_slice(), b"G");
        // Even an empty requested table still reads/validates the source.
        std::fs::write(&path, b">a\nAC\xffGT\n").unwrap();
        assert!(matches!(load_transcript_sequences(&path, &[]),
            Err(Error::Io(error)) if error.kind() == std::io::ErrorKind::InvalidData));
    }
}
