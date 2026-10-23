//! Lossless directory widths; tags only select candidates for full verification.
use super::{OFFSET_MASK, TAG_MASK, WordSpan};

#[cfg_attr(test, derive(Debug))]
pub(super) enum WordDirectory {
    Wide(Box<[u64]>),
    Tiny11 {
        entries: Box<[u16]>,
        spans: Box<[u64]>,
    },
    Dense16 {
        entries: Box<[u32]>,
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
    pub(super) fn payload_bytes(keys: &[WordSpan], narrow: bool) -> usize {
        if narrow {
            keys.len() * 4 + keys.iter().filter(|key| key.end != 0).count() * 8
        } else {
            keys.len() * 8
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
            Self::Wide(entries.into_boxed_slice())
        }
    }
    pub(super) fn build_dense(keys: &[WordSpan]) -> Self {
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
            assert!(key.word <= 0xffff);
            let descriptor = if key.end == 0 {
                assert!(key.first < 32768);
                key.first
            } else {
                assert!(spans.len() < 32768);
                let slot = spans.len() as u32;
                spans.push((u64::from(key.first) << 24) | (u64::from(key.end) << 44));
                32768 | slot
            };
            entries[node] = key.word | (descriptor << 16);
            *rank += 1;
            arrange(keys, entries, spans, node * 2 + 2, rank);
        }
        arrange(keys, &mut entries, &mut spans, 0, &mut 0);
        Self::Dense16 {
            entries: entries.into_boxed_slice(),
            spans: spans.into_boxed_slice(),
        }
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
            Self::Wide(_) => entry,
            Self::Tiny11 { spans, .. } => {
                let descriptor = (entry >> 11) as usize;
                let tag = (entry & 0x7ff) << 13;
                if descriptor < 16 {
                    tag | ((descriptor as u64) << 24)
                } else {
                    tag | spans[descriptor & 15]
                }
            }
            Self::Dense16 { spans, .. } => {
                let descriptor = (entry >> 16) as usize;
                let tag = (entry & 0xffff) << 8;
                if descriptor < 32768 {
                    tag | ((descriptor as u64) << 24)
                } else {
                    tag | spans[descriptor & 32767]
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
            Self::Wide(entries) => search::<_, TAG_MASK>(entries, tag),
            Self::Tiny11 { entries, .. } => search::<_, 0x7ff>(entries, tag >> 13),
            Self::Narrow { entries, .. } => search::<_, TAG_MASK>(entries, tag),
            Self::Dense16 { entries, .. } => search::<_, 0xffff>(entries, tag >> 8),
        };
        result.map(|entry| self.expand(entry))
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
    }
}

#[cfg(test)]
mod dense_tests {
    use super::*;
    #[test]
    fn dense_inline_and_list_descriptors_preserve_all_postings() {
        let keys: Vec<_> = (0..32768u32)
            .map(|i| WordSpan {
                word: i * 2,
                first: if i == 32767 { OFFSET_MASK as u32 } else { i },
                end: if i == 32767 { OFFSET_MASK as u32 } else { 2 },
            })
            .collect();
        let dense = WordDirectory::build_dense(&keys);
        for key in &keys {
            let expected = (u64::from(key.word) << 8)
                | (u64::from(key.first) << 24)
                | (u64::from(key.end) << 44);
            for low in [0, 1, 255] {
                assert_eq!(
                    dense.entry((u64::from(key.word) << 8) | low),
                    Some(expected)
                );
                assert_eq!(dense.entry((u64::from(key.word + 1) << 8) | low), None);
            }
        }
        for ordinal in [0, 127, 128, 4095, 32767] {
            let dense = WordDirectory::build_dense(&[WordSpan {
                word: 65535,
                first: ordinal,
                end: 0,
            }]);
            assert_eq!(
                dense.entry(TAG_MASK),
                Some(0xffff00 | (u64::from(ordinal) << 24))
            );
        }
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
}
