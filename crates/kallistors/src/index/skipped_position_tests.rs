use croaring_sys::{
    roaring_bitmap_deserialize, roaring_bitmap_free, roaring_bitmap_get_cardinality,
    roaring_bitmap_to_uint32_array,
};

use super::*;

#[test]
fn direct_array_decode_equals_c_bitmap_set_semantics() {
    for ids in [
        vec![],
        vec![7u32],
        vec![u32::MAX],
        vec![0x8000_0000],
        vec![7, 7],
        vec![1, 7],
        vec![7, 1],
        vec![1, 7, u32::MAX],
        vec![7, 1, 7],
        vec![9, 3, 9, 1],
        vec![0, u32::MAX, 100_000, 65_535, 65_536],
    ] {
        let mut bytes = vec![1u8];
        bytes.extend_from_slice(&(ids.len() as u32).to_le_bytes());
        for id in ids {
            bytes.extend_from_slice(&id.to_le_bytes());
        }
        let expected = unsafe {
            let bitmap = roaring_bitmap_deserialize(bytes.as_ptr().cast());
            assert!(!bitmap.is_null());
            let mut expected = vec![0; roaring_bitmap_get_cardinality(bitmap) as usize];
            roaring_bitmap_to_uint32_array(bitmap, expected.as_mut_ptr());
            roaring_bitmap_free(bitmap);
            expected
        };
        assert_eq!(
            unsafe { deserialize_roaring_to_vec(&bytes) }.unwrap(),
            expected
        );
        assert_eq!(
            unsafe { deserialize_roaring_to_box(&bytes) }
                .unwrap()
                .as_ref(),
            expected.as_slice()
        );
        bytes.extend_from_slice(&[19, 7, 3]);
        assert_eq!(
            unsafe { deserialize_roaring_to_vec(&bytes) }.unwrap(),
            expected
        );
        assert_eq!(
            unsafe { deserialize_roaring_to_box(&bytes) }
                .unwrap()
                .as_ref(),
            expected.as_slice()
        );
        let positions = deserialize_position_set(&bytes).unwrap();
        match positions {
            PositionSet::Single(value) => assert_eq!(expected, [value]),
            PositionSet::Many(values) => assert_eq!(expected, values.as_ref()),
        }
    }
    for bytes in [
        vec![],
        vec![1],
        vec![1, 1, 0, 0, 0],
        vec![1, 2, 0, 0, 0, 5, 0, 0, 0],
    ] {
        assert!(unsafe { deserialize_roaring_to_vec(&bytes) }.is_none());
        assert!(deserialize_position_set(&bytes).is_none());
        assert!(unsafe { deserialize_roaring_to_box(&bytes) }.is_none());
    }
}

fn sparse_bytes(sizes: &[usize]) -> Vec<u8> {
    let mut bitmap = vec![1u8];
    bitmap.extend_from_slice(&3u32.to_le_bytes());
    for id in [1u32, 7, 33] {
        bitmap.extend_from_slice(&id.to_le_bytes());
    }
    let mut bytes = Vec::new();
    bytes.extend_from_slice(&(bitmap.len() as u64).to_le_bytes());
    bytes.extend_from_slice(&bitmap);
    bytes.extend_from_slice(&(sizes.len() as u64).to_le_bytes());
    for &size in sizes {
        bytes.extend_from_slice(&(size as u64).to_le_bytes());
        bytes.extend(std::iter::repeat_n(91u8, size));
    }
    bytes
}

#[test]
fn slice_blocks_match_generic_blocks_and_preserve_existing_prefix() {
    for count in 0..=3u64 {
        let mut bytes = if count == 0 {
            vec![0]
        } else if count == 1 {
            vec![1]
        } else {
            let mut bytes = vec![2];
            bytes.extend_from_slice(&count.to_le_bytes());
            bytes
        };
        for i in 0..count {
            bytes.extend_from_slice(&(i as u32).to_le_bytes());
            bytes.extend_from_slice(&((i + 1) as u32).to_le_bytes());
            bytes.extend_from_slice(&sparse_bytes(&[1, 1025, 7]));
        }
        let expected: Vec<_> = read_block_array_blocks(&mut Cursor::new(&bytes))
            .unwrap()
            .into_iter()
            .map(|b| (b.lb, b.ub, b.ec.into_vec()))
            .collect();
        let mut actual = Vec::new();
        read_block_array_slice(&bytes, |lb, ub, bitmap| {
            let mut ids = vec![9, 9];
            unsafe { append_roaring_to_vec(bitmap, &mut ids) }.unwrap();
            assert_eq!(&ids[..2], &[9, 9]);
            actual.push((lb, ub, ids[2..].to_vec()));
            Ok(())
        })
        .unwrap();
        assert_eq!(actual, expected);
        for cut in [0, bytes.len() - 1] {
            let old = read_block_array_blocks(&mut Cursor::new(&bytes[..cut]));
            let new = read_block_array_slice(&bytes[..cut], |_, _, bitmap| {
                unsafe { deserialize_roaring_to_vec(bitmap) }.ok_or_else(|| {
                    Error::InvalidFormat("failed to deserialize roaring bitmap".into())
                })?;
                Ok(())
            });
            assert_eq!(old.is_err(), new.is_err());
        }
    }
}

#[test]
fn skipped_positions_consume_exact_bytes_and_reject_truncated_payloads() {
    let groups = [
        vec![7u32],
        vec![0x8000_0001, 0x8000_0002],
        (0..16).collect(),
    ];
    let mut bytes = Vec::new();
    let append_bitmap = |bytes: &mut Vec<u8>, values: &[u32]| {
        let size = 5 + values.len() * 4;
        bytes.extend_from_slice(&(size as u64).to_le_bytes());
        bytes.push(1);
        bytes.extend_from_slice(&(values.len() as u32).to_le_bytes());
        for value in values {
            bytes.extend_from_slice(&value.to_le_bytes());
        }
    };
    append_bitmap(&mut bytes, &[1, 7, 33]);
    // The sparse vector carries its position-group count before bitmaps.
    bytes.extend_from_slice(&(groups.len() as u64).to_le_bytes());
    for group in &groups {
        append_bitmap(&mut bytes, group);
    }
    let length = bytes.len();
    bytes.push(173);
    let mut reader = Cursor::new(bytes.as_slice());
    let (ids, positions) = read_sparse_vector_with_positions(&mut reader).unwrap();
    assert_eq!(ids.as_ref(), [1, 7, 33]);
    assert_eq!(reader.position() as usize, length);
    assert_eq!(read_u8(&mut reader).unwrap(), 173);
    for (position, expected) in positions.iter().zip(&groups) {
        let actual = match position {
            PositionSetRef::Single(v) => vec![v],
            PositionSetRef::Many(v) => v.to_vec(),
        };
        assert_eq!(&actual, expected);
    }
    for cut in [0, 7, length - 1] {
        assert!(read_sparse_vector_with_positions(&mut Cursor::new(&bytes[..cut])).is_err());
    }
    // Four singleton groups use the fused header/maximum pass. Put the
    // largest coordinate last so ignoring any group changes the oracle.
    for words in [
        [0, 0x8000_0001, 0x7ffe, 0x8000_7fff],
        [0, 1, 2, 0x8000_8000],
        [0, 1, 2, 0x807f_ffff],
        [0, 1, 2, 0xffff_ffff],
    ] {
        let mut bytes = Vec::new();
        append_bitmap(&mut bytes, &[1, 7, 33, 100]);
        bytes.extend_from_slice(&4u64.to_le_bytes());
        for word in words {
            append_bitmap(&mut bytes, &[word]);
        }
        let (_, positions) =
            read_sparse_vector_with_positions(&mut Cursor::new(bytes.as_slice())).unwrap();
        for (position, expected) in positions.iter().zip(words) {
            assert!(position.contains(expected));
            assert_eq!(
                position.strand_summary(),
                Some(u8::from(expected & 0x8000_0000 == 0))
            );
        }
        assert!(
            read_sparse_vector_with_positions(&mut Cursor::new(&bytes[..bytes.len() - 1])).is_err()
        );
    }
    for sizes in [[0usize, 0, 0], [1, 1024, 1025], [4099, 3, 5]] {
        let bytes = sparse_bytes(&sizes);
        let mut input = bytes.clone();
        input.push(173);
        let mut reader = std::io::Cursor::new(input);
        assert_eq!(read_sparse_vector_ec(&mut reader).unwrap(), [1, 7, 33]);
        assert_eq!(reader.position() as usize, bytes.len());
        assert_eq!(read_u8(&mut reader).unwrap(), 173);
        for cut in [0, 7, bytes.len() - 1] {
            assert!(read_sparse_vector_ec(&mut std::io::Cursor::new(&bytes[..cut])).is_err());
        }
    }
}
