use std::fs::File;
use std::io::{BufReader, Read, Seek, SeekFrom, Write};
use std::path::PathBuf;

use anyhow::Result;

use kallistors::index::bifrost::{BooPhf, read_bitcontainer_values};

pub fn run(index: PathBuf, out: PathBuf) -> Result<()> {
    let file = File::open(&index)?;
    let mut reader = BufReader::new(file);

    let _index_version = read_u64_le(&mut reader)?;
    let dbg_size = read_u64_le(&mut reader)?;
    let dbg_start = reader.stream_position()?;
    let _file_format_version = read_u64_le(&mut reader)?;
    let _k = read_i32_le(&mut reader)?;
    let _g = read_i32_le(&mut reader)?;

    let v_unitigs_sz = read_u64_le(&mut reader)? as usize;
    for _ in 0..v_unitigs_sz {
        let len = read_u64_le(&mut reader)? as usize;
        let data_sz = len.div_ceil(4);
        reader.seek(SeekFrom::Current(data_sz as i64))?;
    }

    let km_unitigs_sz = read_u64_le(&mut reader)? as usize;
    for _ in 0..km_unitigs_sz {
        reader.seek(SeekFrom::Current(8))?;
    }

    let h_kmers_ccov_sz = read_u64_le(&mut reader)? as usize;
    for _ in 0..h_kmers_ccov_sz {
        reader.seek(SeekFrom::Current(8))?;
    }

    let _bfg_header = read_u64_le(&mut reader)?;
    let _checksum = read_u64_le(&mut reader)?;
    let _v_unitigs_sz2 = read_u64_le(&mut reader)?;
    let _km_unitigs_sz2 = read_u64_le(&mut reader)?;
    let _h_kmers_sz2 = read_u64_le(&mut reader)?;
    let _hmap_min_unitigs_sz = read_u64_le(&mut reader)?;

    let nb_bmp_unitigs = read_u64_le(&mut reader)? as usize;
    for _ in 0..nb_bmp_unitigs {
        let _ = read_bitcontainer_values(&mut reader)?;
    }

    let nb_bmp_km = read_u64_le(&mut reader)? as usize;
    for _ in 0..nb_bmp_km {
        let _ = read_bitcontainer_values(&mut reader)?;
    }

    let nb_special_minz = read_u64_le(&mut reader)? as usize;
    let nb_bmp_special = (nb_special_minz >> 32) + 1;
    for _ in 0..nb_bmp_special {
        let _ = read_bitcontainer_values(&mut reader)?;
    }
    for _ in 0..nb_bmp_special {
        let _ = read_bitcontainer_values(&mut reader)?;
    }
    for _ in 0..nb_special_minz {
        let _minz = read_u64_le(&mut reader)?;
        let _ = read_u64_le(&mut reader)?;
    }

    reader.seek(SeekFrom::Start(dbg_start + dbg_size))?;
    let mphf_size = read_u64_le(&mut reader)? as usize;
    let mut mphf_buf = vec![0u8; mphf_size];
    reader.read_exact(&mut mphf_buf)?;
    let mut mphf_cur = std::io::Cursor::new(mphf_buf.as_slice());
    let mphf = BooPhf::load(&mut mphf_cur)?;

    let entries = mphf.final_hash_entries();
    let mut writer = std::io::BufWriter::new(File::create(out)?);
    writeln!(writer, "mphf_key_hex\tvalue")?;
    for (key, value) in entries {
        let hex = key.iter().map(|b| format!("{b:02x}")).collect::<String>();
        writeln!(writer, "{hex}\t{value}")?;
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
