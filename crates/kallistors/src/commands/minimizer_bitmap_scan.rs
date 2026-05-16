use std::fs::File;
use std::io::{BufReader, Read, Seek, SeekFrom, Write};
use std::path::PathBuf;

use anyhow::{Result, anyhow};

use kallistors::index::bifrost::{BooPhf, encode_minimizer_rep, read_bitcontainer_values};

pub fn run(index: PathBuf, minimizer: String, out: Option<PathBuf>, limit: usize) -> Result<()> {
    let file = File::open(&index)?;
    let mut reader = BufReader::new(file);

    let _index_version = read_u64_le(&mut reader)?;
    let dbg_size = read_u64_le(&mut reader)?;
    let dbg_start = reader.stream_position()?;
    let _file_format_version = read_u64_le(&mut reader)?;
    let k = read_i32_le(&mut reader)? as usize;
    let g = read_i32_le(&mut reader)? as usize;

    if minimizer.len() != g {
        return Err(anyhow!(
            "minimizer length {} does not match g {}",
            minimizer.len(),
            g
        ));
    }
    let target_bytes = encode_minimizer_rep(minimizer.as_bytes())
        .ok_or_else(|| anyhow!("invalid minimizer sequence"))?;

    let v_unitigs_sz = read_u64_le(&mut reader)? as usize;
    let mut unitigs = Vec::with_capacity(v_unitigs_sz);
    for _ in 0..v_unitigs_sz {
        let len = read_u64_le(&mut reader)? as usize;
        let data_sz = len.div_ceil(4);
        let mut data = vec![0u8; data_sz];
        reader.read_exact(&mut data)?;
        unitigs.push(decode_compressed_sequence(&data, len));
    }

    let km_unitigs_sz = read_u64_le(&mut reader)? as usize;
    let mut km_unitigs = Vec::with_capacity(km_unitigs_sz);
    for _ in 0..km_unitigs_sz {
        let mut buf = vec![0u8; 8];
        reader.read_exact(&mut buf)?;
        km_unitigs.push(decode_kmer(&buf, k));
    }

    let h_kmers_ccov_sz = read_u64_le(&mut reader)? as usize;
    for _ in 0..h_kmers_ccov_sz {
        let mut buf = vec![0u8; 8];
        reader.read_exact(&mut buf)?;
        let _ = decode_kmer(&buf, k);
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
    let mut bmp_km = Vec::with_capacity(nb_bmp_km);
    for _ in 0..nb_bmp_km {
        bmp_km.push(read_bitcontainer_values(&mut reader)?);
    }

    reader.seek(SeekFrom::Start(dbg_start + dbg_size))?;
    let mphf_size = read_u64_le(&mut reader)? as usize;
    let mut mphf_buf = vec![0u8; mphf_size];
    reader.read_exact(&mut mphf_buf)?;
    let mut mphf_cur = std::io::Cursor::new(mphf_buf.as_slice());
    let mphf = BooPhf::load(&mut mphf_cur)?;
    let (lookup, mphf_debug) = mphf.debug_lookup(&target_bytes);
    let bytes_hex = target_bytes
        .iter()
        .map(|b| format!("{b:02x}"))
        .collect::<Vec<_>>()
        .join("");

    let mut prefix = Vec::with_capacity(unitigs.len());
    let mut acc = 0u64;
    for seq in &unitigs {
        let count = seq.len().saturating_sub(g) + 1;
        acc += count as u64;
        prefix.push(acc);
    }

    let mut matches = 0usize;
    let mut written = 0usize;
    let mut unitig_scanned = 0usize;
    let mut km_scanned = 0usize;

    let mut out_writer: Option<Box<dyn Write>> = match out {
        Some(path) => Some(Box::new(std::io::BufWriter::new(File::create(path)?))),
        None => None,
    };
    if let Some(writer) = out_writer.as_mut() {
        writeln!(
            writer,
            "unitig_id\tunitig_pos\tis_km\tmphf_lookup\tmphf_key_hex"
        )?;
    }

    for (i, bitmap) in bmp_unitigs.iter().enumerate() {
        let base = (i as u64) << 32;
        for &pos_bmp in bitmap {
            unitig_scanned += 1;
            let pos = base + pos_bmp as u64;
            let unitig_id = match prefix.binary_search(&(pos + 1)) {
                Ok(idx) => idx,
                Err(idx) => idx,
            };
            let prev = if unitig_id == 0 {
                0
            } else {
                prefix[unitig_id - 1]
            };
            let rel = (pos - prev) as usize;
            let unitig = &unitigs[unitig_id];
            if rel + g > unitig.len() {
                continue;
            }
            let slice = &unitig[rel..rel + g];
            let Some(rep) = encode_minimizer_rep(slice) else {
                continue;
            };
            if rep != target_bytes {
                continue;
            }
            matches += 1;
            if let Some(writer) = out_writer.as_mut()
                && written < limit
            {
                writeln!(
                    writer,
                    "{}\t{}\t{}\t{}\t{}",
                    unitig_id,
                    rel,
                    0,
                    lookup
                        .map(|v| v.to_string())
                        .unwrap_or_else(|| "-".to_string()),
                    bytes_hex
                )?;
                written += 1;
            }
        }
    }

    let km_glen = k.saturating_sub(g) + 1;
    if km_glen > 0 {
        for (i, bitmap) in bmp_km.iter().enumerate() {
            let base = (i as u64) << 32;
            for &pos_bmp in bitmap {
                km_scanned += 1;
                let pos = base + pos_bmp as u64;
                let km_id = (pos / km_glen as u64) as usize;
                let km_pos = (pos % km_glen as u64) as usize;
                if km_id >= km_unitigs.len() || km_pos + g > k {
                    continue;
                }
                let slice = &km_unitigs[km_id][km_pos..km_pos + g];
                let Some(rep) = encode_minimizer_rep(slice) else {
                    continue;
                };
                if rep != target_bytes {
                    continue;
                }
                matches += 1;
                if let Some(writer) = out_writer.as_mut()
                    && written < limit
                {
                    writeln!(
                        writer,
                        "{}\t{}\t{}\t{}\t{}",
                        km_id,
                        km_pos,
                        1,
                        lookup
                            .map(|v| v.to_string())
                            .unwrap_or_else(|| "-".to_string()),
                        bytes_hex
                    )?;
                    written += 1;
                }
            }
        }
    }

    println!("minimizer\t{minimizer}");
    println!("mphf_key_hex\t{bytes_hex}");
    println!(
        "mphf_lookup\t{}",
        lookup
            .map(|v| v.to_string())
            .unwrap_or_else(|| "-".to_string())
    );
    println!("mphf_level\t{}", mphf_debug.level);
    println!("mphf_hash_raw\t{}", mphf_debug.hash_raw);
    println!("mphf_hash_domain\t{}", mphf_debug.hash_domain);
    println!(
        "mphf_bucket\t{}",
        mphf_debug
            .bucket
            .map(|v| v.to_string())
            .unwrap_or_else(|| "-".to_string())
    );
    println!(
        "mphf_rank\t{}",
        mphf_debug
            .rank
            .map(|v| v.to_string())
            .unwrap_or_else(|| "-".to_string())
    );
    println!(
        "mphf_final_hash_hit\t{}",
        if mphf_debug.final_hash_hit { 1 } else { 0 }
    );
    println!("matches_in_bitmaps\t{matches}");
    println!("bmp_unitig_positions_scanned\t{unitig_scanned}");
    println!("bmp_km_positions_scanned\t{km_scanned}");
    if out_writer.is_some() {
        println!("positions_dumped\t{written}");
    }

    Ok(())
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

fn decode_kmer(bytes: &[u8], k: usize) -> Vec<u8> {
    let mut out = Vec::with_capacity(k);
    let mut longs: Vec<u64> = Vec::with_capacity(bytes.len() / 8);
    for chunk in bytes.chunks_exact(8) {
        let mut arr = [0u8; 8];
        arr.copy_from_slice(chunk);
        longs.push(u64::from_le_bytes(arr));
    }
    for i in 0..k {
        let idx = i >> 5;
        let shift = 62 - ((i & 0x1F) << 1);
        let v = (longs[idx] >> shift) & 0x3;
        out.push(match v {
            0 => b'A',
            1 => b'C',
            2 => b'G',
            _ => b'T',
        });
    }
    out
}
