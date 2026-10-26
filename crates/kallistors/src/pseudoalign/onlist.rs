use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom, Write};

use crate::index::bifrost::deserialize_roaring_to_vec;
use crate::{Error, Result};

use super::io::{read_i32_le, read_u64_le};

pub(super) struct OnlistData {
    pub(super) onlist: Option<Vec<bool>>,
    pub(super) lengths: Vec<u32>,
    pub(super) names: Vec<String>,
    pub(super) shade_to_color: Vec<Option<u32>>,
    pub(super) shade_sequences: Vec<bool>,
    pub(super) use_shade: bool,
}

impl OnlistData {
    fn empty() -> Self {
        Self {
            onlist: None,
            lengths: Vec::new(),
            names: Vec::new(),
            shade_to_color: Vec::new(),
            shade_sequences: Vec::new(),
            use_shade: false,
        }
    }
}
/// Names retained eagerly for compatibility loaders, or read from the same
/// open index file after mapping and EM have released their working storage.
#[doc(hidden)]
pub struct TranscriptNames {
    eager: Vec<String>,
    deferred: Option<(File, u64, usize)>,
}
impl TranscriptNames {
    pub(super) fn eager(names: Vec<String>) -> Self {
        Self {
            eager: names,
            deferred: None,
        }
    }
    pub(super) fn deferred(file: File, start: u64, count: usize) -> Self {
        Self {
            eager: Vec::new(),
            deferred: Some((file, start, count)),
        }
    }
    /// Borrow already loaded names; deferred owners keep an empty eager view.
    #[doc(hidden)]
    pub fn eager_slice(&self) -> &[String] {
        &self.eager
    }
    /// Stream selected transcript names without retaining a second name table.
    pub fn write_lines(self, writer: &mut impl Write, onlist: Option<&[bool]>) -> Result<usize> {
        let selected = |id: usize| onlist.is_none_or(|mask| mask.get(id).copied().unwrap_or(false));
        let mut written = 0;
        if let Some((file, start, count)) = self.deferred {
            let mut reader = BufReader::with_capacity(1 << 20, file);
            reader.seek(SeekFrom::Start(start))?;
            let mut bytes = Vec::new();
            for id in 0..count {
                let len = read_u64_le(&mut reader)? as usize;
                bytes.resize(len, 0);
                reader.read_exact(&mut bytes)?;
                if selected(id) {
                    writer.write_all(String::from_utf8_lossy(&bytes).as_bytes())?;
                    writer.write_all(b"\n")?;
                    written += 1;
                }
            }
        } else {
            for (id, name) in self.eager.iter().enumerate() {
                if selected(id) {
                    writer.write_all(name.as_bytes())?;
                    writer.write_all(b"\n")?;
                    written += 1;
                }
            }
        }
        Ok(written)
    }
    pub fn load(self) -> Result<Vec<String>> {
        let Some((file, start, count)) = self.deferred else {
            return Ok(self.eager);
        };
        let mut reader = BufReader::with_capacity(1 << 20, file);
        reader.seek(SeekFrom::Start(start))?;
        let mut names = Vec::with_capacity(count);
        for _ in 0..count {
            names.push(read_name(&mut reader)?);
        }
        Ok(names)
    }
}
fn read_name<R: Read>(reader: &mut R) -> Result<String> {
    let len = read_u64_le(reader)? as usize;
    let mut buf = vec![0; len];
    reader.read_exact(&mut buf)?;
    Ok(match String::from_utf8(buf) {
        Ok(name) => name,
        Err(err) => String::from_utf8_lossy(err.as_bytes()).into_owned(),
    })
}
fn read_lengths<R: Read>(reader: &mut R) -> Result<Option<Vec<u32>>> {
    let num_trans = match read_i32_le(reader) {
        Ok(n) => n,
        Err(_) => return Ok(None),
    };
    if num_trans < 0 {
        return Err(Error::InvalidFormat("negative transcript count".into()));
    }
    let mut lengths = Vec::with_capacity(num_trans as usize);
    for _ in 0..num_trans {
        let len = read_i32_le(reader)?;
        if len < 0 {
            return Err(Error::InvalidFormat("negative transcript length".into()));
        }
        lengths.push(len as u32);
    }
    Ok(Some(lengths))
}
fn read_bitmap<R: Read>(reader: &mut R, count: usize) -> Result<Option<Vec<bool>>> {
    let size = read_u64_le(reader)?;
    if size == 0 {
        return Ok(None);
    }
    let mut buf = vec![0; size as usize];
    reader.read_exact(&mut buf)?;
    let vals = unsafe { deserialize_roaring_to_vec(&buf) }
        .ok_or_else(|| Error::InvalidFormat("invalid onlist bitmap".into()))?;
    let mut onlist = vec![false; count];
    for v in vals {
        if let Some(slot) = onlist.get_mut(v as usize) {
            *slot = true;
        }
    }
    Ok(Some(onlist))
}
// Scan every complete name for shade annotations before deferring. Shade indices
// use the eager parser, including its encounter-order color resolution.
pub(super) fn read_onlist_scanned<R: BufRead + Seek>(
    reader: &mut R,
) -> Result<(OnlistData, Option<(u64, usize)>)> {
    let start = reader.stream_position()?;
    let Some(lengths) = read_lengths(reader)? else {
        return Ok((OnlistData::empty(), None));
    };
    let count = lengths.len();
    let names_start = reader.stream_position()?;
    let finder = memchr::memmem::Finder::new(b"_shade_");
    let mut scratch = Vec::new();
    for _ in 0..count {
        let len = read_u64_le(reader)? as usize;
        let available = reader.fill_buf()?;
        let shade = if available.len() >= len {
            let shade = finder.find(&available[..len]).is_some();
            reader.consume(len);
            shade
        } else {
            scratch.resize(len, 0);
            reader.read_exact(&mut scratch)?;
            finder.find(&scratch).is_some()
        };
        if shade {
            reader.seek(SeekFrom::Start(start))?;
            return read_onlist_with_lengths(reader).map(|data| (data, None));
        }
    }
    let onlist = read_bitmap(reader, count)?;
    Ok((
        OnlistData {
            onlist,
            lengths,
            names: Vec::new(),
            shade_to_color: Vec::new(),
            shade_sequences: Vec::new(),
            use_shade: false,
        },
        Some((names_start, count)),
    ))
}

fn split_shade_name(name: &str) -> Option<(&str, &str)> {
    let needle = "_shade_";
    name.find(needle)
        .map(|idx| (&name[..idx], &name[idx + needle.len()..]))
}

pub(super) fn read_onlist<R: Read>(reader: &mut R) -> Result<Option<Vec<bool>>> {
    read_onlist_with_lengths(reader).map(|data| data.onlist)
}

pub(super) fn read_onlist_with_lengths<R: Read>(reader: &mut R) -> Result<OnlistData> {
    let Some(lengths) = read_lengths(reader)? else {
        return Ok(OnlistData::empty());
    };
    let num_trans = lengths.len();
    let mut names: Vec<String> = Vec::with_capacity(num_trans);
    let mut shade_to_color = Vec::new();
    let mut shade_sequences = Vec::new();
    let mut name_to_index: Option<HashMap<String, usize>> = None;
    let mut use_shade = false;
    for idx in 0..num_trans {
        let name = read_name(reader)?;
        if let Some((base, _variant)) = split_shade_name(&name) {
            if !use_shade {
                use_shade = true;
                shade_to_color = vec![None; num_trans];
                shade_sequences = vec![false; num_trans];
            }
            shade_sequences[idx] = true;
            let map = name_to_index.get_or_insert_with(|| {
                let mut map = HashMap::with_capacity(names.len());
                for (name_idx, prior_name) in names.iter().enumerate() {
                    if split_shade_name(prior_name).is_none() {
                        map.entry(prior_name.clone()).or_insert(name_idx);
                    }
                }
                map
            });
            if let Some(color_idx) = map.get(base) {
                shade_to_color[idx] = Some(*color_idx as u32);
            }
        } else if let Some(map) = name_to_index.as_mut() {
            map.entry(name.clone()).or_insert(idx);
        }
        names.push(name);
    }

    let onlist = read_bitmap(reader, num_trans)?;
    Ok(OnlistData {
        onlist,
        lengths,
        names,
        shade_to_color,
        shade_sequences,
        use_shade,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;
    #[test]
    fn deferred_names_keep_utf8_ids_and_eager_shade_encounter_order() {
        for names in [
            vec![],
            vec![b"ENST1|full|table".as_slice(), b"invalid\xff".as_slice()],
            vec![
                b"x_shade_1".as_slice(),
                b"x".as_slice(),
                b"x_shade_2".as_slice(),
            ],
        ] {
            let mut bytes = Vec::new();
            bytes.extend_from_slice(&(names.len() as i32).to_le_bytes());
            for _ in &names {
                bytes.extend_from_slice(&100i32.to_le_bytes());
            }
            for name in &names {
                bytes.extend_from_slice(&(name.len() as u64).to_le_bytes());
                bytes.extend_from_slice(name);
            }
            bytes.extend_from_slice(&0u64.to_le_bytes());
            let eager = read_onlist_with_lengths(&mut Cursor::new(&bytes)).unwrap();
            let temp = tempfile::NamedTempFile::new().unwrap();
            std::fs::write(temp.path(), &bytes).unwrap();
            for mask in [None, Some(vec![]), Some(vec![true, false, true])] {
                let expected: Vec<_> = eager
                    .names
                    .iter()
                    .enumerate()
                    .filter(|(id, _)| {
                        mask.as_ref()
                            .is_none_or(|mask| mask.get(*id).copied().unwrap_or(false))
                    })
                    .map(|(_, name)| name)
                    .collect();
                let text: String = expected.iter().map(|name| format!("{name}\n")).collect();
                for source in [
                    TranscriptNames::eager(eager.names.clone()),
                    TranscriptNames::deferred(
                        File::open(temp.path()).unwrap(),
                        4 + 4 * names.len() as u64,
                        names.len(),
                    ),
                ] {
                    if source.deferred.is_some() {
                        assert!(source.eager_slice().is_empty());
                    } else {
                        assert_eq!(source.eager_slice(), eager.names.as_slice());
                    }
                    let mut actual = Vec::new();
                    let count = source.write_lines(&mut actual, mask.as_deref()).unwrap();
                    assert_eq!(actual, text.as_bytes());
                    assert_eq!(count, expected.len());
                }
            }

            for capacity in [7, 1024] {
                let mut reader = BufReader::with_capacity(capacity, Cursor::new(&bytes));
                let (scanned, deferred) = read_onlist_scanned(&mut reader).unwrap();
                assert_eq!(scanned.lengths, eager.lengths);
                assert_eq!(scanned.onlist, eager.onlist);
                assert_eq!(scanned.shade_to_color, eager.shade_to_color);
                assert_eq!(scanned.shade_sequences, eager.shade_sequences);
                if let Some((start, count)) = deferred {
                    reader.seek(SeekFrom::Start(start)).unwrap();
                    let actual: Vec<_> = (0..count)
                        .map(|_| read_name(&mut reader).unwrap())
                        .collect();
                    assert_eq!(actual, eager.names);
                } else {
                    assert_eq!(scanned.names, eager.names);
                    assert_eq!(scanned.shade_to_color, [None, None, Some(1)]);
                }
            }
        }
    }
}
