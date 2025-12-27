use std::collections::HashMap;
use std::io::Read;

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

fn split_shade_name(name: &str) -> Option<(&str, &str)> {
    let needle = "_shade_";
    name.find(needle)
        .map(|idx| (&name[..idx], &name[idx + needle.len()..]))
}

pub(super) fn read_onlist<R: Read>(reader: &mut R) -> Result<Option<Vec<bool>>> {
    read_onlist_with_lengths(reader).map(|data| data.onlist)
}

pub(super) fn read_onlist_with_lengths<R: Read>(reader: &mut R) -> Result<OnlistData> {
    let num_trans = match read_i32_le(reader) {
        Ok(v) => v,
        Err(_) => {
            return Ok(OnlistData {
                onlist: None,
                lengths: Vec::new(),
                names: Vec::new(),
                shade_to_color: Vec::new(),
                shade_sequences: Vec::new(),
                use_shade: false,
            });
        }
    };
    if num_trans < 0 {
        return Err(Error::InvalidFormat("negative transcript count".into()));
    }
    let num_trans = num_trans as usize;

    let mut lengths = Vec::with_capacity(num_trans);
    for _ in 0..num_trans {
        let len = read_i32_le(reader)?;
        if len < 0 {
            return Err(Error::InvalidFormat("negative transcript length".into()));
        }
        lengths.push(len as u32);
    }
    let mut names = Vec::with_capacity(num_trans);
    let mut shade_to_color = vec![None; num_trans];
    let mut shade_sequences = vec![false; num_trans];
    let mut name_to_index: HashMap<String, usize> = HashMap::new();
    let mut use_shade = false;
    for idx in 0..num_trans {
        let name_len = read_u64_le(reader)? as usize;
        let mut buf = vec![0u8; name_len];
        reader.read_exact(&mut buf)?;
        let name = String::from_utf8_lossy(&buf).into_owned();
        if let Some((base, _variant)) = split_shade_name(&name) {
            use_shade = true;
            shade_sequences[idx] = true;
            if let Some(color_idx) = name_to_index.get(base) {
                shade_to_color[idx] = Some(*color_idx as u32);
            }
        } else {
            name_to_index.entry(name.clone()).or_insert(idx);
        }
        names.push(name);
    }

    let onlist_size = read_u64_le(reader)?;
    if onlist_size == 0 {
        return Ok(OnlistData {
            onlist: None,
            lengths,
            names,
            shade_to_color,
            shade_sequences,
            use_shade,
        });
    }
    let mut buf = vec![0u8; onlist_size as usize];
    reader.read_exact(&mut buf)?;
    let vals = unsafe { deserialize_roaring_to_vec(&buf) }
        .ok_or_else(|| Error::InvalidFormat("invalid onlist bitmap".into()))?;
    let mut onlist = vec![false; num_trans];
    for v in vals {
        let idx = v as usize;
        if idx < onlist.len() {
            onlist[idx] = true;
        }
    }
    Ok(OnlistData {
        onlist: Some(onlist),
        lengths,
        names,
        shade_to_color,
        shade_sequences,
        use_shade,
    })
}
