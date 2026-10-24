use std::collections::HashSet;
use std::io::Cursor;

use super::*;
use crate::index::bifrost::BooPhf;

#[test]
fn compressed_sequence_append_preserves_native_bits_and_tail() {
    for len in 0..128 {
        let sequence = b"ACGTNacgu-"
            .iter()
            .copied()
            .cycle()
            .take(len)
            .collect::<Vec<_>>();
        let mut expected = vec![0u8; len.div_ceil(4)];
        for (i, &base) in sequence.iter().enumerate() {
            let code = match base {
                b'A' => 0u8,
                b'C' => 1,
                b'G' => 2,
                _ => 3,
            };
            expected[i / 4] |= code << (2 * (i % 4));
        }
        let mut actual = vec![0x12, 0x34];
        let range = append_compressed_sequence(&mut actual, &sequence);
        assert_eq!(&actual[..2], &[0x12, 0x34]);
        assert_eq!(&actual[range], expected);
    }
}

#[test]
fn small_bitmap_records_match_croaring() {
    let mut cases = vec![
        vec![],
        vec![0],
        vec![u32::MAX],
        vec![0, u32::MAX],
        vec![1, 1],
        vec![2, 1],
        vec![1, 2, 3],
        (0..6).collect(),
        (0..7).collect(),
        vec![0, 65536, 131072, 196608, 262144, 327680],
    ];
    for len in [3, 4, 5, 6, 7, 8, 16, 31, 255, 1024, 4096, 4097] {
        cases.push((0..len).map(|i| 65536 + 2 * i as u32).collect());
        cases.push((0..len).map(|i| i as u32).collect());
        cases.push((0..len).map(|i| (i as u32) << 16).collect());
    }
    cases.push(vec![65534, 65535, 65536]);
    cases.push(vec![65533, 65535, 65537, u32::MAX]);
    cases.push(vec![u32::MAX - 4, u32::MAX - 2, u32::MAX]);
    for size in [3, 4, 5, 6, 7, 8, 31, 4096, 4097, 65536] {
        cases.push((0..size).collect());
        cases.push((0..size).map(|i| i + i / 3).collect());
        cases.push((0..size.min(32000)).map(|i| 65536 + i + i / 7).collect());
    }
    for values in cases {
        for optimize in [false, true] {
            let mut actual = Vec::new();
            if optimize {
                write_prefixed_bitmap::<true>(&mut actual, &values).unwrap();
            } else {
                write_prefixed_bitmap::<false>(&mut actual, &values).unwrap();
            }
            let native = serialize_roaring(&values, optimize, false).unwrap();
            let mut expected = (native.len() as u64).to_le_bytes().to_vec();
            expected.extend(native);
            assert_eq!(actual, expected, "values={values:?}, optimize={optimize}");
        }
    }
}

#[test]
fn boo_mphf_roundtrips_all_keys_without_all_final_hash() {
    let keys = (0..2000u64)
        .map(|i| u64::from_be_bytes(i.wrapping_mul(0x9e37_79b9_7f4a_7c15).to_le_bytes()))
        .collect::<Vec<_>>();

    let mut buf = Vec::new();
    write_boo_mphf(&mut buf, &keys, 1).expect("write mphf");
    let mut parallel = Vec::new();
    write_boo_mphf(&mut parallel, &keys, 8).expect("parallel mphf");
    assert_eq!(parallel, buf);
    for threads in [1, 8] {
        let mut implicit = Vec::new();
        write_boo_mphf_ids_first::<u32, u32>(&mut implicit, &keys, threads, true).unwrap();
        assert_eq!(implicit, buf);
        let mut wide = Vec::new();
        write_boo_mphf_ids::<usize, u64>(&mut wide, &keys, threads).expect("wide survivor IDs");
        assert_eq!(wide, buf);
        let mut wide_buckets = Vec::new();
        write_boo_mphf_ids::<u32, u64>(&mut wide_buckets, &keys, threads).unwrap();
        assert_eq!(wide_buckets, buf);
    }
    assert_eq!(
        <u32 as MphfBucket>::from_bucket(u32::MAX as u64).bucket(),
        u32::MAX as u64
    );
    assert_eq!(
        <u64 as MphfBucket>::from_bucket(u64::MAX).bucket(),
        u64::MAX
    );
    assert_eq!(
        <u32 as MphfIndex>::from_index(u32::MAX as usize).index(),
        u32::MAX as usize
    );
    assert_eq!(
        <usize as MphfIndex>::from_index(usize::MAX).index(),
        usize::MAX
    );

    let mphf = BooPhf::load(&mut Cursor::new(buf)).expect("load mphf");
    assert_eq!(mphf.size(), keys.len() as u64);
    assert!(mphf.final_hash_entries().len() < keys.len());

    let mut ranks = HashSet::new();
    for key in &keys {
        let rank = mphf.lookup(&key.to_be_bytes()).expect("lookup key");
        assert!(rank < keys.len() as u64);
        assert!(ranks.insert(rank), "duplicate rank {rank}");
    }
    assert_eq!(ranks.len(), keys.len());
}
