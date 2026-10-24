//! Exact19-bit conservative selector in the existing32-bit entry width.

use super::{OFFSET_MASK, WordDirectory, WordSpan, hashed_entries};

impl WordDirectory {
    pub(in crate::pseudoalign::bifrost::candidate_cache) fn dense19_layout(
        keys: &[WordSpan],
    ) -> bool {
        keys.iter().all(|key| {
            key.word <= 0x7ffff
                && if key.end == 0 {
                    key.first < 4096
                } else {
                    u64::from(key.first) <= OFFSET_MASK && u64::from(key.end) <= OFFSET_MASK
                }
        }) && keys.iter().filter(|key| key.end != 0).count() <= 4096
    }

    pub(in crate::pseudoalign::bifrost::candidate_cache) fn build_dense19(
        keys: &[WordSpan],
    ) -> Option<Self> {
        if !Self::dense19_layout(keys) {
            return None;
        }
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
                key.first
            } else {
                let slot = spans.len() as u32;
                spans.push((u64::from(key.first) << 24) | (u64::from(key.end) << 44));
                4096 | slot
            };
            entries[node] = key.word | (descriptor << 19);
            *rank += 1;
            arrange(keys, entries, spans, node * 2 + 2, rank);
        }
        arrange(keys, &mut entries, &mut spans, 0, &mut 0);
        let spans = spans.into_boxed_slice();
        Some(match hashed_entries::<_, 0x7ffff>(&entries, u32::MAX) {
            Some(entries) => Self::HashedDense19 { entries, spans },
            None => Self::Dense19 {
                entries: entries.into_boxed_slice(),
                spans,
            },
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dense19_matches_independent_sorted_descriptors_and_checked_fallback() {
        let keys: Vec<_> = (0..4096u32)
            .map(|i| WordSpan {
                word: if i == 4095 { 0x7ffff } else { i * 3 },
                first: if i == 4095 { OFFSET_MASK as u32 } else { i },
                end: if i == 4095 { OFFSET_MASK as u32 } else { 2 },
            })
            .collect();
        let directory = WordDirectory::build_dense19(&keys).unwrap();
        // Last entry is all ones, forcing the tree/sentinel fallback.
        assert!(matches!(directory, WordDirectory::Dense19 { .. }));
        for word in (0..=12286u32).chain([0x7fffe, 0x7ffff]) {
            let expected = keys.iter().find(|key| key.word == word).map(|key| {
                (u64::from(key.word) << 5)
                    | (u64::from(key.first) << 24)
                    | (u64::from(key.end) << 44)
            });
            for low in [0, 1, 31] {
                assert_eq!(directory.entry((u64::from(word) << 5) | low), expected);
            }
        }
        for ordinal in [0, 127, 128, 2047, 2048, 4095] {
            let directory = WordDirectory::build_dense19(&[WordSpan {
                word: 0x7ffff,
                first: ordinal,
                end: 0,
            }])
            .unwrap();
            assert!(matches!(directory, WordDirectory::HashedDense19 { .. }));
            assert_eq!(
                directory.entry(super::super::TAG_MASK),
                Some(0xffffe0 | (u64::from(ordinal) << 24))
            );
        }
        for bad in [
            WordSpan {
                word: 0x80000,
                first: 0,
                end: 0,
            },
            WordSpan {
                word: 0,
                first: 4096,
                end: 0,
            },
            WordSpan {
                word: 0,
                first: (OFFSET_MASK + 1) as u32,
                end: 2,
            },
            WordSpan {
                word: 0,
                first: 0,
                end: (OFFSET_MASK + 1) as u32,
            },
        ] {
            assert!(WordDirectory::build_dense19(&[bad]).is_none());
        }
        let too_many: Vec<_> = (0..4097)
            .map(|word| WordSpan {
                word,
                first: 0,
                end: 2,
            })
            .collect();
        assert!(WordDirectory::build_dense19(&too_many).is_none());
    }
}
