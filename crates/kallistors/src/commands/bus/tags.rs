//! Barcode and UMI extraction, tag matching and binary encoding.

use super::*;

pub(super) fn normalize_corrected_barcode(mut barcode: Vec<u8>) -> Vec<u8> {
    if let Some(pos) = barcode.iter().position(|base| *base == b'-') {
        barcode.truncate(pos);
    }
    barcode
}

pub(super) fn barcode_sequence(
    batch: &[kallistors::io::FastqRecord],
    spec: &TechnologySpec,
    batch_index: usize,
    batch_barcodes: bool,
) -> Result<Option<Vec<u8>>> {
    let extracted = if spec.bc.first().is_some_and(|v| v.file.is_none()) {
        fake_barcode(batch_index as u64, 16)
    } else {
        let Some(barcode) = concat_slices_optional(batch, &spec.bc, false)? else {
            return Ok(None);
        };
        barcode
    };
    if batch_barcodes && spec.bc.first().is_none_or(|v| v.file.is_some()) {
        if extracted.len() >= 32 {
            bail!("--batch-barcodes requires barcode length shorter than 32 bases");
        }
        let prefix_len = 32usize.saturating_sub(extracted.len()).min(32);
        let mut prefixed = fake_barcode(batch_index as u64, prefix_len);
        prefixed.extend_from_slice(&extracted);
        Ok(Some(prefixed))
    } else {
        Ok(Some(extracted))
    }
}

pub(super) fn configure_tag(
    technology: &str,
    tag_sequence: Option<&str>,
    spec: &mut TechnologySpec,
) -> Result<Option<ConfiguredTag>> {
    let inferred_default = tag_sequence.is_none() && technology.eq_ignore_ascii_case("SMARTSEQ3");
    let tag = tag_sequence
        .map(str::to_string)
        .or_else(|| inferred_default.then(|| "ATTGCGCAATG".to_string()));
    let Some(tag) = tag else {
        return Ok(None);
    };
    validate_tag_sequence(&tag)?;
    let Some(first_umi) = spec.umi.first_mut() else {
        bail!("tag sequence requires a UMI slice");
    };
    let Some(stop) = first_umi.stop else {
        bail!("tag sequence requires a bounded UMI slice");
    };
    first_umi.start = first_umi.start.saturating_add(tag.len());
    if first_umi.start >= stop {
        bail!("tag sequence must be shorter than UMI sequence");
    }
    spec.tag_len = Some(tag.len());
    let mut flag = 0;
    Ok(Some(ConfiguredTag {
        config: TagConfig {
            binary: string_to_binary(tag.as_bytes(), &mut flag, "tag sequence")?,
            sequence: tag.into_bytes(),
        },
        inferred_default,
    }))
}

pub(super) fn validate_tag_sequence(tag: &str) -> Result<()> {
    if tag.is_empty() {
        bail!("tag sequence cannot be empty");
    }
    for (idx, base) in tag.bytes().enumerate() {
        if !matches!(base, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't') {
            bail!(
                "tag sequence contains invalid base {} at position {}",
                char::from(base),
                idx + 1
            );
        }
    }
    Ok(())
}

pub(super) fn umi_value(
    batch: &[kallistors::io::FastqRecord],
    spec: &TechnologySpec,
    tag: Option<&TagConfig>,
) -> Result<Option<UmiValue>> {
    if spec.bulk_like() {
        let mut flag = 0;
        return Ok(Some(UmiValue {
            binary: string_to_binary(b"A", &mut flag, "UMI")?,
            flag,
            len: 1,
            seq: b"A".to_vec(),
            tag_present: false,
            disable_strand_specific: false,
        }));
    }
    if spec.keep_fastq_comments {
        let Some(umi_seq) = rx_umi_from_header(
            batch
                .first()
                .ok_or_else(|| anyhow!("RX UMI extraction requires at least one FASTQ record"))?,
        ) else {
            return Ok(None);
        };
        let mut flag = 0;
        return Ok(Some(UmiValue {
            binary: string_to_binary(&umi_seq, &mut flag, "UMI")?,
            flag,
            len: umi_seq.len(),
            seq: umi_seq,
            tag_present: false,
            disable_strand_specific: false,
        }));
    }

    let Some(tag) = tag else {
        let Some(umi_seq) = concat_slices_optional(batch, &spec.umi, false)? else {
            return Ok(None);
        };
        let mut flag = 0;
        return Ok(Some(UmiValue {
            binary: string_to_binary(&umi_seq, &mut flag, "UMI")?,
            flag,
            len: umi_seq.len(),
            seq: umi_seq,
            tag_present: false,
            disable_strand_specific: false,
        }));
    };

    let Some(first_umi) = spec.umi.first() else {
        bail!("tag sequence requires a UMI slice");
    };
    let Some(mut tagged_umi) = extract_tagged_umi_optional(batch, *first_umi, tag.sequence.len())?
    else {
        return Ok(None);
    };
    for spec in spec.umi.iter().skip(1) {
        let Some(slice) = extract_slice_optional(batch, *spec)? else {
            return Ok(None);
        };
        tagged_umi.extend_from_slice(&slice);
    }
    let mut flag = 0;
    let tagged_binary = string_to_binary(&tagged_umi, &mut flag, "tagged UMI")?;
    let tag_len = tag.sequence.len();
    let hamming_threshold = usize::from(tag_len > 5);
    let tag_present = tagged_umi.len() >= tag_len
        && hamming(
            tag.binary,
            tagged_binary >> (2 * (tagged_umi.len() - tag_len)),
            tag_len,
        ) <= hamming_threshold;
    if tag_present {
        let umi_len = tagged_umi.len() - tag_len;
        let binary = if umi_len == 0 {
            0
        } else {
            tagged_binary & !(u64::MAX << (2 * umi_len))
        };
        Ok(Some(UmiValue {
            binary,
            flag,
            len: umi_len,
            seq: tagged_umi[tag_len..].to_vec(),
            tag_present: true,
            disable_strand_specific: false,
        }))
    } else {
        Ok(Some(UmiValue {
            binary: u64::MAX,
            flag: 0,
            len: 0,
            seq: Vec::new(),
            tag_present: false,
            disable_strand_specific: true,
        }))
    }
}

pub(super) fn extract_tagged_umi_optional(
    records: &[kallistors::io::FastqRecord],
    spec: SliceSpec,
    tag_len: usize,
) -> Result<Option<Vec<u8>>> {
    let mut tagged_spec = spec;
    tagged_spec.start = tagged_spec.start.saturating_sub(tag_len);
    extract_slice_optional(records, tagged_spec)
}

pub(super) fn rx_umi_from_header(record: &kallistors::io::FastqRecord) -> Option<Vec<u8>> {
    let header = record.header.strip_prefix(b"@").unwrap_or(&record.header);
    let pos = header.windows(5).position(|window| window == b"RX:Z:")?;
    let start = pos + 5;
    let end = header[start..]
        .iter()
        .position(|base| matches!(base, b' ' | b'\t'))
        .map_or(header.len(), |offset| start + offset);
    (end > start).then(|| header[start..end].to_vec())
}

pub(super) fn string_to_binary(seq: &[u8], flag: &mut u32, label: &str) -> Result<u64> {
    if seq.len() > 32 {
        bail!("{label} length {} exceeds BUS limit of 32 bases", seq.len());
    }
    *flag = 0;
    let mut result = 0u64;
    let mut num_n = 0u32;
    let mut pos_n = 0u32;
    let (chunks, tail) = seq.as_chunks::<8>();
    for (chunk, bytes) in chunks.iter().enumerate() {
        let word = u64::from_be_bytes(*bytes);
        // These bit classes deliberately preserve kallisto's behavior for
        // every byte, including lowercase and ambiguous/non-DNA bytes.
        let classes = (word & 0x0303_0303_0303_0303) ^ 0x0202_0202_0202_0202;
        let flagged = !((classes & 0x7f7f_7f7f_7f7f_7f7f).wrapping_add(0x7f7f_7f7f_7f7f_7f7f)
            | classes)
            & 0x8080_8080_8080_8080;
        if flagged != 0 {
            if num_n == 0 {
                pos_n = (chunk * 8) as u32 + flagged.leading_zeros() / 8;
            }
            num_n += flagged.count_ones();
        }
        let x = (word & 0x0404_0404_0404_0404) >> 1;
        let mut codes = x | ((x ^ (word & 0x0202_0202_0202_0202)) >> 1);
        codes = (codes | (codes >> 6)) & 0x000f_000f_000f_000f;
        codes = (codes | (codes >> 12)) & 0x0000_00ff_0000_00ff;
        codes = (codes | (codes >> 24)) & 0xffff;
        result = (result << 16) | codes;
    }
    let tail_start = chunks.len() * 8;
    for (offset, base) in tail.iter().enumerate() {
        let idx = tail_start + offset;
        let x = (base & 4) >> 1;
        if (base & 3) == 2 {
            if num_n == 0 {
                pos_n = u32::try_from(idx).unwrap_or(0);
            }
            num_n += 1;
        }
        result <<= 2;
        result |= u64::from(x + ((x ^ (base & 2)) >> 1));
    }
    if num_n > 0 {
        *flag = (num_n.min(3) & 3) | ((pos_n & 31) << 2);
    }
    Ok(result)
}

pub(super) fn hamming(left: u64, right: u64, len: usize) -> usize {
    let diff = left ^ right;
    let mut distance = 0usize;
    for idx in 0..len {
        let shift = 2 * (len - idx - 1);
        if ((diff >> shift) & 0x03) != 0 {
            distance += 1;
        }
    }
    distance
}

#[cfg(test)]
mod binary_word_tests {
    use super::*;
    #[test]
    fn word_encoding_and_ambiguous_flags_match_original_byte_rules() {
        fn original(bytes: &[u8]) -> (u64, u32) {
            let mut word = 0;
            let mut count = 0u32;
            let mut first = 0;
            for (i, &base) in bytes.iter().enumerate() {
                let x = (base & 4) >> 1;
                word = (word << 2) | u64::from(x + ((x ^ (base & 2)) >> 1));
                if base & 3 == 2 {
                    if count == 0 {
                        first = i as u32;
                    }
                    count += 1;
                }
            }
            (
                word,
                if count == 0 {
                    0
                } else {
                    count.min(3) | ((first & 31) << 2)
                },
            )
        }
        for length in [0, 1, 7, 8, 9, 12, 16, 31, 32] {
            for byte in 0..=255u8 {
                for offset in 0..length {
                    let mut bytes = vec![b'A'; length];
                    bytes[offset] = byte;
                    if length > 3 {
                        bytes[length - 2] = b'N';
                    }
                    let mut flag = 99;
                    let actual = string_to_binary(&bytes, &mut flag, "test").unwrap();
                    assert_eq!(
                        (actual, flag),
                        original(&bytes),
                        "length={length} byte={byte} offset={offset}"
                    );
                }
                let bytes = vec![byte; length];
                let mut flag = 99;
                let actual = string_to_binary(&bytes, &mut flag, "test").unwrap();
                assert_eq!((actual, flag), original(&bytes));
            }
        }
        let mut flag = 99;
        assert!(string_to_binary(&[b'A'; 33], &mut flag, "test").is_err());
        assert_eq!(flag, 99);
    }
}
