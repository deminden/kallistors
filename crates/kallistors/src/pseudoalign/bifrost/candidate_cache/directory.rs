//! Lossless directory widths; tags only select candidates for full verification.
mod dense19;
use super::{OFFSET_MASK, TAG_MASK, WordSpan};

#[cfg_attr(test, derive(Debug))]
pub(super) enum WordDirectory {
    Wide(Box<[u64]>),
    HashedWide(Box<[u64]>),
    Dense19 {
        entries: Box<[u32]>,
        spans: Box<[u64]>,
    },
    HashedDense19 {
        entries: Box<[u32]>,
        spans: Box<[u64]>,
    },
    TinyMasks {
        entries: Box<[u16]>,
        masks: Box<[u16]>,
    },
    Tiny11 {
        entries: Box<[u16]>,
        spans: Box<[u64]>,
    },
    Narrow {
        entries: Box<[u32]>,
        spans: Box<[u64]>,
    },
}
impl Default for WordDirectory {
    fn default() -> Self {
        Self::Wide(Box::default())
    }
}
impl WordDirectory {
    pub(super) fn narrow_layout(keys: &[WordSpan], source_len: usize) -> bool {
        let multiple = keys.iter().filter(|key| key.end != 0).count();
        source_len <= 128 && multiple <= 128 && multiple * 8 < keys.len() * 4
    }
    pub(super) fn payload_bytes(keys: &[WordSpan], narrow: bool, hashed: bool) -> usize {
        let count = if hashed {
            hash_slots(keys.len())
        } else {
            keys.len()
        };
        if narrow {
            count * 4 + keys.iter().filter(|key| key.end != 0).count() * 8
        } else {
            count * 8
        }
    }
    pub(super) fn build(keys: &[WordSpan], narrow: bool) -> Self {
        if narrow {
            let mut entries = vec![0u32; keys.len()];
            let mut spans = Vec::with_capacity(keys.iter().filter(|key| key.end != 0).count());
            fn arrange(
                keys: &[WordSpan],
                entries: &mut [u32],
                spans: &mut Vec<u64>,
                node: usize,
                rank: &mut usize,
            ) {
                if node >= keys.len() {
                    return;
                }
                arrange(keys, entries, spans, node * 2 + 1, rank);
                let key = &keys[*rank];
                let descriptor = if key.end == 0 {
                    debug_assert!(key.first < 128);
                    key.first
                } else {
                    debug_assert!(spans.len() < 128);
                    let slot = spans.len() as u32;
                    spans.push((u64::from(key.first) << 24) | (u64::from(key.end) << 44));
                    128 | slot
                };
                entries[node] = key.word | (descriptor << 24);
                *rank += 1;
                arrange(keys, entries, spans, node * 2 + 2, rank);
            }
            arrange(keys, &mut entries, &mut spans, 0, &mut 0);
            Self::Narrow {
                entries: entries.into_boxed_slice(),
                spans: spans.into_boxed_slice(),
            }
        } else {
            let mut entries = vec![0u64; keys.len()];
            fn arrange(keys: &[WordSpan], entries: &mut [u64], node: usize, rank: &mut usize) {
                if node >= keys.len() {
                    return;
                }
                arrange(keys, entries, node * 2 + 1, rank);
                let key = &keys[*rank];
                debug_assert!(u64::from(key.first) <= OFFSET_MASK);
                entries[node] =
                    u64::from(key.word) | (u64::from(key.first) << 24) | (u64::from(key.end) << 44);
                *rank += 1;
                arrange(keys, entries, node * 2 + 2, rank);
            }
            arrange(keys, &mut entries, 0, &mut 0);
            match hashed_entries::<_, TAG_MASK>(&entries, u64::MAX) {
                Some(entries) => Self::HashedWide(entries),
                None => Self::Wide(entries.into_boxed_slice()),
            }
        }
    }

    #[cfg(test)]
    pub(super) fn retained_payload_bytes(&self) -> usize {
        match self {
            Self::Wide(entries) | Self::HashedWide(entries) => entries.len() * 8,
            Self::Tiny11 { entries, spans } => entries.len() * 2 + spans.len() * 8,
            Self::TinyMasks { entries, masks } => entries.len() * 2 + masks.len() * 2,
            Self::Narrow { entries, spans }
            | Self::Dense19 { entries, spans }
            | Self::HashedDense19 { entries, spans } => entries.len() * 4 + spans.len() * 8,
        }
    }

    pub(super) fn tiny_masks_supported(keys: &[WordSpan], ordinals: &[u16]) -> bool {
        if !Self::tiny_layout(keys) {
            return false;
        }
        let mut previous = None;
        keys.iter().all(|key| {
            if key.word > 0x7ff || previous.is_some_and(|value| value >= key.word) {
                return false;
            }
            previous = Some(key.word);
            if key.end == 0 {
                return key.first < 16;
            }
            let first = key.first as usize;
            let Some(last) = first.checked_add(key.end as usize) else {
                return false;
            };
            let Some(list) = ordinals.get(first..last) else {
                return false;
            };
            !list.is_empty()
                && list.iter().all(|&value| value < 16)
                && list.windows(2).all(|pair| pair[0] < pair[1])
        })
    }

    pub(super) fn tiny_mask_payload_bytes(keys: &[WordSpan]) -> usize {
        keys.len() * 2 + keys.iter().filter(|key| key.end != 0).count() * 2
    }

    pub(super) fn build_tiny_masks(keys: &[WordSpan], ordinals: &[u16]) -> Option<Self> {
        if !Self::tiny_masks_supported(keys, ordinals) {
            return None;
        }
        let mut entries = vec![0u16; keys.len()];
        let mut masks = Vec::with_capacity(keys.iter().filter(|key| key.end != 0).count());
        fn arrange(
            keys: &[WordSpan],
            ordinals: &[u16],
            entries: &mut [u16],
            masks: &mut Vec<u16>,
            node: usize,
            rank: &mut usize,
        ) {
            if node >= keys.len() {
                return;
            }
            arrange(keys, ordinals, entries, masks, node * 2 + 1, rank);
            let key = &keys[*rank];
            let descriptor = if key.end == 0 {
                key.first as u16
            } else {
                let slot = masks.len() as u16;
                let first = key.first as usize;
                let mask = ordinals[first..first + key.end as usize]
                    .iter()
                    .fold(0u16, |mask, &ordinal| mask | (1u16 << ordinal));
                masks.push(mask);
                16 | slot
            };
            entries[node] = key.word as u16 | (descriptor << 11);
            *rank += 1;
            arrange(keys, ordinals, entries, masks, node * 2 + 2, rank);
        }
        arrange(keys, ordinals, &mut entries, &mut masks, 0, &mut 0);
        Some(Self::TinyMasks {
            entries: entries.into_boxed_slice(),
            masks: masks.into_boxed_slice(),
        })
    }

    #[inline]
    pub(super) fn tiny_mask(
        &self,
        forward: u64,
        reverse: u64,
        allow_forward: bool,
        allow_reverse: bool,
    ) -> Option<u16> {
        let Self::TinyMasks { entries, masks } = self else {
            return None;
        };
        let get = |tag: u64| -> u16 {
            let Some(entry) = search::<_, 0x7ff>(entries, tag >> 13) else {
                return 0;
            };
            let descriptor = (entry >> 11) as usize;
            if descriptor < 16 {
                1u16 << descriptor
            } else {
                masks[descriptor & 15]
            }
        };
        Some(
            if allow_forward { get(forward) } else { 0 }
                | if allow_reverse { get(reverse) } else { 0 },
        )
    }

    pub(super) fn tiny_layout(keys: &[WordSpan]) -> bool {
        keys.iter().filter(|key| key.end != 0).count() <= 16
    }
    pub(super) fn tiny_payload_bytes(keys: &[WordSpan]) -> usize {
        keys.len() * 2 + keys.iter().filter(|key| key.end != 0).count() * 8
    }
    pub(super) fn build_tiny(keys: &[WordSpan]) -> Self {
        let mut entries = vec![0u16; keys.len()];
        let mut spans = Vec::new();
        fn arrange(
            keys: &[WordSpan],
            entries: &mut [u16],
            spans: &mut Vec<u64>,
            node: usize,
            rank: &mut usize,
        ) {
            if node >= keys.len() {
                return;
            }
            arrange(keys, entries, spans, node * 2 + 1, rank);
            let key = &keys[*rank];
            assert!(key.word <= 0x7ff);
            let descriptor = if key.end == 0 {
                assert!(key.first < 16);
                key.first as u16
            } else {
                assert!(spans.len() < 16);
                let slot = spans.len() as u16;
                spans.push((u64::from(key.first) << 24) | (u64::from(key.end) << 44));
                16 | slot
            };
            entries[node] = (key.word as u16) | (descriptor << 11);
            *rank += 1;
            arrange(keys, entries, spans, node * 2 + 2, rank);
        }
        arrange(keys, &mut entries, &mut spans, 0, &mut 0);
        Self::Tiny11 {
            entries: entries.into_boxed_slice(),
            spans: spans.into_boxed_slice(),
        }
    }
    fn expand(&self, entry: u64) -> u64 {
        match self {
            Self::TinyMasks { .. } => unreachable!("mask directory has no ordinal spans"),
            Self::Wide(_) | Self::HashedWide(_) => entry,
            Self::Tiny11 { spans, .. } => {
                let descriptor = (entry >> 11) as usize;
                let tag = (entry & 0x7ff) << 13;
                if descriptor < 16 {
                    tag | ((descriptor as u64) << 24)
                } else {
                    tag | spans[descriptor & 15]
                }
            }
            Self::Dense19 { spans, .. } | Self::HashedDense19 { spans, .. } => {
                let descriptor = (entry >> 19) as usize;
                let tag = (entry & 0x7ffff) << 5;
                if descriptor < 4096 {
                    tag | ((descriptor as u64) << 24)
                } else {
                    tag | spans[descriptor & 4095]
                }
            }
            Self::Narrow { spans, .. } => {
                let descriptor = (entry >> 24) as usize;
                let tag = entry & TAG_MASK;
                if descriptor < 128 {
                    tag | ((descriptor as u64) << 24)
                } else {
                    tag | spans[descriptor & 127]
                }
            }
        }
    }
    pub(super) fn entry(&self, tag: u64) -> Option<u64> {
        let result = match self {
            Self::TinyMasks { .. } => unreachable!("mask directory requires mask traversal"),
            Self::Wide(entries) => search::<_, TAG_MASK>(entries, tag),
            Self::HashedWide(entries) => hash_search::<_, TAG_MASK>(entries, tag, u64::MAX),
            Self::HashedDense19 { entries, .. } => {
                hash_search::<_, 0x7ffff>(entries, tag >> 5, u32::MAX)
            }
            Self::Dense19 { entries, .. } => search::<_, 0x7ffff>(entries, tag >> 5),
            Self::Tiny11 { entries, .. } => search::<_, 0x7ff>(entries, tag >> 13),
            Self::Narrow { entries, .. } => search::<_, TAG_MASK>(entries, tag),
        };
        result.map(|entry| self.expand(entry))
    }
}

// A table has at least one empty slot and at most 3/4 nominal occupancy.
// Complete selector values still decide membership; a hash collision merely
// probes another slot. Real all-ones entries retain the tree representation.
fn hash_slots(entries: usize) -> usize {
    if entries == 0 {
        return 0;
    }
    entries
        .checked_add(entries.div_ceil(3))
        .and_then(usize::checked_next_power_of_two)
        .expect("bounded candidate directory")
}

fn hashed_entries<T: Copy + Into<u64> + PartialEq, const MASK: u64>(
    entries: &[T],
    empty: T,
) -> Option<Box<[T]>> {
    if entries.contains(&empty) {
        return None;
    }
    let slots = hash_slots(entries.len());
    let mut table = vec![empty; slots];
    for &entry in entries {
        let mut slot = (entry.into() & MASK) as usize & (slots - 1);
        while table[slot] != empty {
            slot = (slot + 1) & (slots - 1);
        }
        table[slot] = entry;
    }
    Some(table.into_boxed_slice())
}

#[inline]
fn hash_search<T: Copy + Into<u64> + PartialEq, const MASK: u64>(
    entries: &[T],
    tag: u64,
    empty: T,
) -> Option<u64> {
    let mask = entries.len().checked_sub(1)?;
    let mut slot = tag as usize & mask;
    loop {
        let value = entries[slot];
        if value == empty {
            return None;
        }
        let value = value.into();
        if value & MASK == tag {
            return Some(value);
        }
        slot = (slot + 1) & mask;
    }
}

// One-based heap paths encode left/right turns as zero/one. After
// descent, remove trailing right turns and the preceding left turn to
// recover the last lower-bound candidate, then compare its complete tag.
fn search<T: Copy + Into<u64>, const MASK: u64>(entries: &[T], tag: u64) -> Option<u64> {
    let mut position = 1usize;
    while let Some(&entry) = entries.get(position - 1) {
        let entry: u64 = entry.into();
        position = position * 2 + usize::from((entry & MASK) < tag);
    }
    recover_leaf_entry::<_, MASK>(entries, position, tag)
}

#[inline]
fn recover_leaf_entry<T: Copy + Into<u64>, const MASK: u64>(
    entries: &[T],
    position: usize,
    tag: u64,
) -> Option<u64> {
    let node = position
        .checked_shr(position.trailing_ones() + 1)?
        .checked_sub(1)?;
    let entry: u64 = (*entries.get(node)?).into();
    ((entry & MASK) == tag).then_some(entry)
}

#[cfg(test)]
mod tests {
    #[test]
    fn leaf_search_matches_sorted_membership_for_incomplete_heap_levels() {
        for len in [0, 1, 2, 3, 7, 8, 9, 63, 64, 65, 257] {
            let sorted: Vec<_> = (0..len).map(|i| (i * 5 + 1) as u64).collect();
            let mut directory = vec![0; len];
            fn fill(words: &mut [u64], node: usize, rank: &mut usize, sorted: &[u64]) {
                if node >= words.len() {
                    return;
                }
                fill(words, node * 2 + 1, rank, sorted);
                words[node] = sorted[*rank] | ((*rank as u64) << 24);
                *rank += 1;
                fill(words, node * 2 + 2, rank, sorted);
            }
            fill(&mut directory, 0, &mut 0, &sorted);
            let expected = |query: u64| {
                sorted
                    .binary_search(&query)
                    .ok()
                    .map(|i| query | ((i as u64) << 24))
            };
            for query in (0..=(len * 5 + 2) as u64).chain([super::TAG_MASK]) {
                assert_eq!(
                    super::search::<_, { super::TAG_MASK }>(&directory, query),
                    expected(query)
                );
            }
        }

        // Independent sorted membership through wraparound probe clusters,
        // absent tags, empty tables and an actual reserved-sentinel value.
        for len in [0, 1, 2, 3, 7, 8, 9, 63, 64, 65, 257] {
            let values: Vec<u64> = (0..len)
                .map(|i| ((i * 1024 + 1023) as u64) | ((i as u64) << 24))
                .collect();
            let hashed =
                super::hashed_entries::<_, { super::TAG_MASK }>(&values, u64::MAX).unwrap();
            for query in (0..=(len * 1024 + 1024) as u64)
                .step_by(31)
                .chain(values.iter().map(|v| v & super::TAG_MASK))
            {
                let expected = values
                    .iter()
                    .copied()
                    .find(|v| v & super::TAG_MASK == query);
                assert_eq!(
                    super::hash_search::<_, { super::TAG_MASK }>(&hashed, query, u64::MAX,),
                    expected
                );
            }
            assert_eq!(
                super::hash_search::<_, { super::TAG_MASK }>(
                    &hashed,
                    super::TAG_MASK + 1,
                    u64::MAX,
                ),
                None
            );
        }
        assert!(super::hashed_entries::<_, { super::TAG_MASK }>(&[u64::MAX], u64::MAX,).is_none());
        // The tree fallback must retain the complete all-ones descriptor.
        let extreme = super::WordDirectory::build(
            &[super::WordSpan {
                word: super::TAG_MASK as u32,
                first: super::OFFSET_MASK as u32,
                end: super::OFFSET_MASK as u32,
            }],
            false,
        );
        assert_eq!(extreme.entry(super::TAG_MASK), Some(u64::MAX));
    }
}

#[cfg(test)]
mod tiny_tests {
    use super::*;
    #[test]
    fn tiny_selectors_preserve_inline_heads_lists_and_incomplete_trees() {
        for len in [0, 1, 2, 3, 7, 8, 9, 16, 31, 65, 129] {
            let keys: Vec<_> = (0..len)
                .map(|i| WordSpan {
                    word: if i + 1 == len { 0x7ff } else { i * 7 },
                    first: if i == 15 { OFFSET_MASK as u32 } else { i % 16 },
                    end: if i < 16 {
                        if i == 15 { OFFSET_MASK as u32 } else { 2 }
                    } else {
                        0
                    },
                })
                .collect();
            assert!(WordDirectory::tiny_layout(&keys));
            let tiny = WordDirectory::build_tiny(&keys);
            for word in (0..=1000).chain([0x7ff]) {
                let expected = keys.iter().find(|key| key.word == word).map(|key| {
                    (u64::from(key.word) << 13)
                        | (u64::from(key.first) << 24)
                        | (u64::from(key.end) << 44)
                });
                for tail in [0, 1, 8191] {
                    assert_eq!(tiny.entry((u64::from(word) << 13) | tail), expected);
                }
            }
        }
        let unsupported: Vec<_> = (0..17)
            .map(|word| WordSpan {
                word,
                first: 0,
                end: 2,
            })
            .collect();
        assert!(!WordDirectory::tiny_layout(&unsupported));
        let inline = WordDirectory::build_tiny(&[WordSpan {
            word: 0x7ff,
            first: 15,
            end: 0,
        }]);
        assert_eq!(inline.entry(TAG_MASK), Some(0xffe000 | (15 << 24)));
    }
    #[test]
    fn tiny_masks_preserve_order_full_unions_and_reject_multiset_or_wide_inputs() {
        use super::super::{CachedBucket, candidate_tag};
        // All original positions are distinct ordinals even if their encoded
        // values repeat; merging orientations must emit each ordinal once.
        let positions: Vec<u64> = (0..16).map(|id| (id / 2) as u64).collect();
        let ordinals = vec![0, 3, 8, 15, 1, 3, 14];
        let keys = [
            super::WordSpan {
                word: 0,
                first: 0,
                end: 4,
            },
            super::WordSpan {
                word: 13,
                first: 15,
                end: 0,
            },
            super::WordSpan {
                word: 2047,
                first: 4,
                end: 3,
            },
        ];
        let directory = super::WordDirectory::build_tiny_masks(&keys, &ordinals).unwrap();
        assert_eq!(super::WordDirectory::tiny_mask_payload_bytes(&keys), 10);
        assert_eq!(directory.retained_payload_bytes(), 10);
        let bucket = CachedBucket {
            source_addr: positions.as_ptr() as usize,
            source_len: positions.len(),
            directory,
            ordinals: Box::default(),
        };
        // Invert the existing odd tag mix, independently selecting every
        // boundary tag; encoded values remain only conservative selectors.
        let query = |tag: u32| u64::from((tag << 21).wrapping_mul(0x144c_bc89));
        for a in [0u32, 1, 12, 13, 14, 2046, 2047] {
            for b in [0u32, 1, 13, 14, 2047] {
                let forward = query(a);
                let reverse = query(b);
                assert_eq!(candidate_tag(forward) >> 13, a);
                assert_eq!(candidate_tag(reverse) >> 13, b);
                for flags in 0..4 {
                    let mut expected = Vec::new();
                    for (ordinal, &position) in positions.iter().enumerate() {
                        let member = |tag| match tag {
                            0 => [0, 3, 8, 15].contains(&ordinal),
                            13 => ordinal == 15,
                            2047 => [1, 3, 14].contains(&ordinal),
                            _ => false,
                        };
                        if flags & 1 != 0 && member(a) || flags & 2 != 0 && member(b) {
                            expected.push(position);
                        }
                    }
                    let actual: Vec<_> = bucket
                        .matching_mask_candidates(
                            &positions,
                            forward,
                            reverse,
                            flags & 1 != 0,
                            flags & 2 != 0,
                        )
                        .unwrap()
                        .collect();
                    assert_eq!(actual, expected);
                }
            }
        }
        for invalid in [
            vec![0, 0, 8, 15, 1, 3, 14],
            vec![0, 3, 16, 17, 1, 3, 14],
            vec![0, 8, 3, 15, 1, 3, 14],
            vec![0, 3],
        ] {
            assert!(super::WordDirectory::build_tiny_masks(&keys, &invalid).is_none());
        }
        let inline_wide = [super::WordSpan {
            word: 0,
            first: 16,
            end: 0,
        }];
        assert!(super::WordDirectory::build_tiny_masks(&inline_wide, &[]).is_none());
        let boundary_keys: Vec<_> = (0..16)
            .map(|id| super::WordSpan {
                word: if id == 15 { 2047 } else { id },
                first: id,
                end: 1,
            })
            .collect();
        let boundary_ordinals: Vec<u16> = (0..16).collect();
        let boundary =
            super::WordDirectory::build_tiny_masks(&boundary_keys, &boundary_ordinals).unwrap();
        // Descriptor31 plus tag2047 is a real all-ones u16 entry, not empty.
        assert_eq!(
            boundary.tiny_mask(2047 << 13, 0, true, false),
            Some(1 << 15)
        );
        let mut too_many = boundary_keys;
        too_many.insert(
            15,
            super::WordSpan {
                word: 2046,
                first: 0,
                end: 1,
            },
        );
        assert!(super::WordDirectory::build_tiny_masks(&too_many, &boundary_ordinals).is_none());
        let empty = super::WordDirectory::build_tiny_masks(&[], &[]).unwrap();
        assert_eq!(empty.tiny_mask(0, 0, true, true), Some(0));
        // Every possible u16 mask emits exactly its increasing ordinal list.
        for mask in 0..=u16::MAX {
            let actual: Vec<_> = super::super::MaskMatchingCandidates {
                entries: &positions,
                mask,
            }
            .collect();
            let expected: Vec<_> = positions
                .iter()
                .enumerate()
                .filter_map(|(bit, &value)| (mask & (1 << bit) != 0).then_some(value))
                .collect();
            assert_eq!(actual, expected);
        }
    }
}
