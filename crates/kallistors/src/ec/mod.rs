//! Equivalence class handling.

use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use crate::{Error, Result};

/// An EC list in kallisto text format.
#[derive(Debug, Clone)]
pub struct EcList {
    pub classes: Vec<Vec<u32>>,
}

impl EcList {
    /// Load a kallisto EC text file (e.g. `output.ec.txt`).
    pub fn load_text(path: &Path) -> Result<Self> {
        let file = File::open(path).map_err(|_| Error::MissingFile(path.to_path_buf()))?;
        let reader = BufReader::new(file);
        let mut classes: Vec<Vec<u32>> = Vec::new();

        for (line_no, line) in reader.lines().enumerate() {
            let line = line?;
            if line.trim().is_empty() {
                continue;
            }
            let mut parts = line.splitn(2, '\t');
            let ec_id_str = parts.next().ok_or_else(|| {
                Error::InvalidFormat(format!("missing EC id at line {}", line_no + 1))
            })?;
            let ec_id: usize = ec_id_str.parse().map_err(|_| {
                Error::InvalidFormat(format!("invalid EC id at line {}", line_no + 1))
            })?;

            let list_str = parts.next().ok_or_else(|| {
                Error::InvalidFormat(format!("missing EC list at line {}", line_no + 1))
            })?;
            let mut ids: Vec<u32> = Vec::new();
            if !list_str.is_empty() {
                for token in list_str.split(',') {
                    if token.is_empty() {
                        continue;
                    }
                    let val: u32 = token.parse().map_err(|_| {
                        Error::InvalidFormat(format!(
                            "invalid transcript id at line {}",
                            line_no + 1
                        ))
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
        }

        Ok(Self { classes })
    }

    /// Return basic stats: (num_ecs, max_ec_size).
    pub fn stats(&self) -> (usize, usize) {
        let max_size = self.classes.iter().map(|c| c.len()).max().unwrap_or(0);
        (self.classes.len(), max_size)
    }

    /// Write a kallisto EC text file.
    pub fn write_text(&self, path: &Path) -> Result<()> {
        let file = File::create(path)?;
        let mut writer = BufWriter::new(file);
        for (ec_id, ids) in self.classes.iter().enumerate() {
            write!(writer, "{}", ec_id)?;
            write!(writer, "\t")?;
            for (i, val) in ids.iter().enumerate() {
                if i > 0 {
                    write!(writer, ",")?;
                }
                write!(writer, "{}", val)?;
            }
            writeln!(writer)?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn load_ec_list_text() {
        let mut tmp = NamedTempFile::new().unwrap();
        writeln!(tmp, "0\t1,2,3").unwrap();
        writeln!(tmp, "1\t4").unwrap();
        writeln!(tmp, "2\t").unwrap();

        let ec = EcList::load_text(tmp.path()).unwrap();
        assert_eq!(ec.classes.len(), 3);
        assert_eq!(ec.classes[0], vec![1, 2, 3]);
        assert_eq!(ec.classes[1], vec![4]);
        assert_eq!(ec.classes[2], Vec::<u32>::new());
    }

    #[test]
    fn roundtrip_ec_list_text() {
        let mut tmp = NamedTempFile::new().unwrap();
        writeln!(tmp, "0\t2,3").unwrap();
        writeln!(tmp, "1\t").unwrap();
        writeln!(tmp, "2\t5,6,7").unwrap();

        let ec = EcList::load_text(tmp.path()).unwrap();
        let out = NamedTempFile::new().unwrap();
        ec.write_text(out.path()).unwrap();

        let reloaded = EcList::load_text(out.path()).unwrap();
        assert_eq!(ec.classes, reloaded.classes);
    }
}
