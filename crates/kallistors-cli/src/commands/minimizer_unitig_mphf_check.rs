use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, Read, Seek, SeekFrom, Write};
use std::path::PathBuf;

use anyhow::{Result, anyhow};

use kallistors::index::bifrost::{BooPhf, encode_minimizer_rep, read_bitcontainer_values};

pub fn run(index: PathBuf, unitigs: PathBuf, out: Option<PathBuf>) -> Result<()> {
    let unitig_ids = load_unitig_list(&unitigs)?;
    let file = File::open(&index)?;
    let mut reader = BufReader::new(file);

    let _index_version = read_u64_le(&mut reader)?;
    let dbg_size = read_u64_le(&mut reader)?;
    let dbg_start = reader.stream_position()?;
    let _file_format_version = read_u64_le(&mut reader)?;
    let _k = read_i32_le(&mut reader)? as usize;
    let g = read_i32_le(&mut reader)? as usize;

    let v_unitigs_sz = read_u64_le(&mut reader)? as usize;
    let mut unitigs_seq = Vec::with_capacity(v_unitigs_sz);
    for _ in 0..v_unitigs_sz {
        let len = read_u64_le(&mut reader)? as usize;
        let data_sz = len.div_ceil(4);
        let mut data = vec![0u8; data_sz];
        reader.read_exact(&mut data)?;
        unitigs_seq.push(decode_compressed_sequence(&data, len));
    }

    let km_unitigs_sz = read_u64_le(&mut reader)? as usize;
    for _ in 0..km_unitigs_sz {
        let mut buf = vec![0u8; 8];
        reader.read_exact(&mut buf)?;
    }

    let h_kmers_ccov_sz = read_u64_le(&mut reader)? as usize;
    for _ in 0..h_kmers_ccov_sz {
        let mut buf = vec![0u8; 8];
        reader.read_exact(&mut buf)?;
    }

    let _bfg_header = read_u64_le(&mut reader)?;
    let _checksum = read_u64_le(&mut reader)?;
    let _v_unitigs_sz2 = read_u64_le(&mut reader)?;
    let _km_unitigs_sz2 = read_u64_le(&mut reader)?;
    let _h_kmers_sz2 = read_u64_le(&mut reader)?;
    let _hmap_min_unitigs_sz = read_u64_le(&mut reader)?;

    let nb_bmp_unitigs = read_u64_le(&mut reader)? as usize;
    let mut bmp_unitigs = Vec::with_capacity(nb_bmp_unitigs);
    for _ in 0..nb_bmp_unitigs {
        bmp_unitigs.push(read_bitcontainer_values(&mut reader)?);
    }

    let nb_bmp_km = read_u64_le(&mut reader)? as usize;
    for _ in 0..nb_bmp_km {
        let _ = read_bitcontainer_values(&mut reader)?;
    }

    reader.seek(SeekFrom::Start(dbg_start + dbg_size))?;
    let mphf_size = read_u64_le(&mut reader)? as usize;
    let mut mphf_buf = vec![0u8; mphf_size];
    reader.read_exact(&mut mphf_buf)?;
    let mut mphf_cur = std::io::Cursor::new(mphf_buf.as_slice());
    let mphf = BooPhf::load(&mut mphf_cur)?;

    let mut prefix = Vec::with_capacity(unitigs_seq.len());
    let mut acc = 0u64;
    for seq in &unitigs_seq {
        let count = seq.len().saturating_sub(g) + 1;
        acc += count as u64;
        prefix.push(acc);
    }

    let mut counts: HashMap<[u8; 8], usize> = HashMap::new();
    let mut total_positions = 0usize;
    for (i, bitmap) in bmp_unitigs.iter().enumerate() {
        let base = (i as u64) << 32;
        for &pos_bmp in bitmap {
            let pos = base + pos_bmp as u64;
            let unitig_id = match prefix.binary_search(&(pos + 1)) {
                Ok(idx) => idx,
                Err(idx) => idx,
            };
            if !unitig_ids.contains(&unitig_id) {
                continue;
            }
            let prev = if unitig_id == 0 {
                0
            } else {
                prefix[unitig_id - 1]
            };
            let rel = (pos - prev) as usize;
            let unitig = &unitigs_seq[unitig_id];
            if rel + g > unitig.len() {
                continue;
            }
            let slice = &unitig[rel..rel + g];
            let Some(rep) = encode_minimizer_rep(slice) else {
                continue;
            };
            *counts.entry(rep).or_insert(0) += 1;
            total_positions += 1;
        }
    }

    let mut rows: Vec<([u8; 8], usize, bool)> = Vec::new();
    for (rep, count) in counts {
        let hit = mphf.lookup(&rep).is_some();
        rows.push((rep, count, hit));
    }
    rows.sort_by(|a, b| b.1.cmp(&a.1));

    let mut out_writer: Box<dyn Write> = match out {
        Some(path) => Box::new(std::io::BufWriter::new(File::create(path)?)),
        None => Box::new(std::io::BufWriter::new(std::io::stdout())),
    };
    writeln!(out_writer, "minimizer_hex\tcount\tmphf_hit\tminimizer")?;
    for (rep, count, hit) in rows {
        let minimizer_hex = rep.iter().map(|b| format!("{b:02x}")).collect::<String>();
        let minimizer = decode_minimizer_rep(rep, g);
        writeln!(
            out_writer,
            "{}\t{}\t{}\t{}",
            minimizer_hex,
            count,
            if hit { 1 } else { 0 },
            minimizer
        )?;
    }

    if total_positions == 0 {
        return Err(anyhow!("no bitmap positions found for specified unitigs"));
    }
    Ok(())
}

fn load_unitig_list(path: &PathBuf) -> Result<Vec<usize>> {
    let text = std::fs::read_to_string(path)?;
    let mut out = Vec::new();
    for line in text.lines() {
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        out.push(trimmed.parse::<usize>()?);
    }
    out.sort_unstable();
    out.dedup();
    Ok(out)
}

fn decode_minimizer_rep(rep: [u8; 8], g: usize) -> String {
    let val = u64::from_le_bytes(rep);
    let mut out = vec![b'A'; g];
    for (i, slot) in out.iter_mut().enumerate().take(g) {
        let shift = 62 - ((i & 0x1F) << 1);
        let v = ((val >> shift) & 0x3) as u8;
        *slot = match v {
            0 => b'A',
            1 => b'C',
            2 => b'G',
            _ => b'T',
        };
    }
    String::from_utf8_lossy(&out).into_owned()
}

fn read_u64_le<R: Read>(reader: &mut R) -> Result<u64> {
    let mut buf = [0u8; 8];
    reader.read_exact(&mut buf)?;
    Ok(u64::from_le_bytes(buf))
}

fn read_i32_le<R: Read>(reader: &mut R) -> Result<i32> {
    let mut buf = [0u8; 4];
    reader.read_exact(&mut buf)?;
    Ok(i32::from_le_bytes(buf))
}

fn decode_compressed_sequence(data: &[u8], len: usize) -> Vec<u8> {
    let mut out = Vec::with_capacity(len);
    for i in 0..len {
        let byte = data[i >> 2];
        let shift = (i & 0x3) << 1;
        let v = (byte >> shift) & 0x03;
        out.push(match v {
            0 => b'A',
            1 => b'C',
            2 => b'G',
            _ => b'T',
        });
    }
    out
}
