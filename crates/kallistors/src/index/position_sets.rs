//! Full-width grouped positions in one allocation, with compact and wide offsets.
use super::PositionSet;

pub(crate) struct PositionSets {
    groups: usize,
    data: Box<[u32]>,
}
// Adaptive descriptors retain the complete coordinate and strand value.
// Full-u32 descriptors remain available whenever narrower layouts cannot fit. A separate bit per group marks
// offsets into length-prefixed repeated coordinates. The wide fallback retains
// full-width offsets and lengths when the arena exceeds the compact domain.
const WIDE: usize = 1usize << (usize::BITS - 1);
const NARROW: usize = 1usize << (usize::BITS - 2);
const MEDIUM: usize = 1usize << (usize::BITS - 3);
// On 64-bit targets one complete u32 coordinate fits below a fourth layout
// bit. An empty box then needs no allocation; narrower targets retain the
// existing allocated representation.
const INLINE: usize = if usize::BITS > 32 {
    1usize << (usize::BITS - 4)
} else {
    0
};
const INLINE_ARITY_SHIFT: u32 = usize::BITS - 8;
// Proven all-singleton arenas need neither indirect-entry flags nor reads
// from those flags. Other arenas and 32-bit targets keep their old layout.
const ALL_SINGLE: usize = if usize::BITS > 32 {
    1usize << (usize::BITS - 5)
} else {
    0
};
const LAYOUT: usize = WIDE | NARROW | MEDIUM | INLINE | ALL_SINGLE;
const INDIRECT: u64 = 1u64 << 63;

pub(crate) enum PositionSetRef<'a> {
    Single(u32),
    Many(&'a [u32]),
}

impl PositionSets {
    pub(super) fn from_values(sets: Vec<PositionSet>) -> Self {
        Self::pack(sets, false)
    }

    // A validated all-singleton payload can provide exact words directly.
    // Select the same descriptors as pack(), without a transient PositionSet Vec.
    #[cfg(test)]
    pub(super) fn from_singletons(values: impl ExactSizeIterator<Item = u32> + Clone) -> Self {
        let maximum = values
            .clone()
            .map(|raw| raw & 0x7fff_ffff)
            .max()
            .unwrap_or(0);
        Self::from_singletons_with_max(values, maximum)
    }

    pub(super) fn from_singletons_with_max(
        mut values: impl ExactSizeIterator<Item = u32> + Clone,
        max_coordinate: u32,
    ) -> Self {
        let groups = values.len();
        if INLINE != 0 && (1..=3).contains(&groups) {
            let mut words = [0; 3];
            for (slot, value) in words.iter_mut().zip(values.clone()) {
                *slot = value;
            }
            if let Some(packed) = Self::inline_values(&words[..groups]) {
                return packed;
            }
        }
        // Fewer than four groups always use complete u32 descriptors. Build
        // those exact small boxes without width scans or a temporary Vec header.
        if groups < 4 {
            let mut next = || values.next().expect("validated singleton count");
            let data: Box<[u32]> = match groups {
                0 => Box::new([]),
                1 => Box::new([next(), 0]),
                2 => Box::new([next(), next(), 0]),
                3 => Box::new([next(), next(), next(), 0]),
                _ => unreachable!(),
            };
            return Self { groups, data };
        }
        let flags = if ALL_SINGLE != 0 {
            0
        } else {
            groups.div_ceil(32)
        };
        let narrow_prefix = groups
            .div_ceil(2)
            .checked_add(flags)
            .expect("position descriptor capacity");
        let medium_prefix = groups
            .checked_mul(3)
            .expect("position descriptor capacity")
            .div_ceil(4)
            .checked_add(flags)
            .expect("position descriptor capacity");
        let compact_prefix = groups
            .checked_add(flags)
            .expect("position descriptor capacity");
        let narrow = groups >= 4 && max_coordinate <= 0x7fff && narrow_prefix <= u16::MAX as usize;
        let medium =
            !narrow && groups >= 8 && max_coordinate <= 0x7f_ffff && medium_prefix <= 0xff_ffff;
        let wide = !narrow && !medium && u32::try_from(compact_prefix).is_err();
        debug_assert_eq!(groups & LAYOUT, 0);
        let (layout, data) = if wide {
            let size = groups.checked_mul(2).expect("position descriptor capacity");
            (
                WIDE,
                Self::pack_words::<64>(values.map(PositionSet::Single), groups, size, size),
            )
        } else if narrow {
            (
                NARROW,
                Self::pack_words::<16>(
                    values.map(PositionSet::Single),
                    groups,
                    narrow_prefix,
                    narrow_prefix,
                ),
            )
        } else if medium {
            (
                MEDIUM,
                Self::pack_words::<24>(
                    values.map(PositionSet::Single),
                    groups,
                    medium_prefix,
                    medium_prefix,
                ),
            )
        } else {
            (
                0,
                Self::pack_words::<32>(
                    values.map(PositionSet::Single),
                    groups,
                    compact_prefix,
                    compact_prefix,
                ),
            )
        };
        Self {
            groups: groups | layout | ALL_SINGLE,
            data,
        }
    }

    // The payload occupies at most48 bits; the arity tag is above it and
    // below every layout flag. Reject any coordinate that cannot be retained
    // in full, then keep its sign in the highest bit of its own descriptor.
    fn inline_values(values: &[u32]) -> Option<Self> {
        let (bits, coordinate_mask) = match values.len() {
            1 => (32, 0x7fff_ffff),
            2 => (24, 0x7f_ffff),
            3 => (16, 0x7fff),
            _ => return None,
        };
        let mut payload = 0usize;
        for (index, &raw) in values.iter().enumerate() {
            if raw & 0x7fff_ffff > coordinate_mask {
                return None;
            }
            let word = (raw & coordinate_mask) | ((raw >> (32 - bits)) & (1 << (bits - 1)));
            payload |= (word as usize) << (index * bits);
        }
        Some(Self {
            groups: INLINE | ((values.len() - 1) << INLINE_ARITY_SHIFT) | payload,
            data: Box::new([]),
        })
    }

    #[inline]
    fn inline_raw(&self, index: usize, groups: usize) -> u32 {
        let bits = match groups {
            1 => return self.groups as u32,
            2 => 24,
            _ => 16,
        };
        let word = (self.groups >> (index * bits)) as u32;
        let sign = 1u32 << (bits - 1);
        (word & (sign - 1)) | ((word & sign) << (32 - bits))
    }

    fn pack(sets: Vec<PositionSet>, force_wide: bool) -> Self {
        let groups = sets.len();
        if INLINE != 0 && !force_wide && (1..=3).contains(&groups) {
            let mut words = [0; 3];
            let mut all_single = true;
            for (word, set) in words.iter_mut().zip(&sets) {
                if let PositionSet::Single(raw) = set {
                    *word = *raw;
                } else {
                    all_single = false;
                    break;
                }
            }
            if all_single && let Some(packed) = Self::inline_values(&words[..groups]) {
                return packed;
            }
        }
        let mut repeated = 0usize;
        let mut many = 0usize;
        let mut max_coordinate = 0;
        for set in &sets {
            match set {
                PositionSet::Single(raw) => max_coordinate = max_coordinate.max(raw & 0x7fff_ffff),
                PositionSet::Many(values) => {
                    repeated = repeated
                        .checked_add(values.len())
                        .expect("position arena capacity");
                    many += 1;
                }
            }
        }
        let size_with_prefix = |prefix: usize, lengths: usize| {
            prefix
                .checked_add(repeated)
                .and_then(|n| n.checked_add(lengths))
                .expect("position arena capacity")
        };
        let flags = groups.div_ceil(32);
        let narrow_prefix = groups
            .div_ceil(2)
            .checked_add(flags)
            .expect("position descriptor capacity");
        let narrow_size = size_with_prefix(narrow_prefix, many);
        let medium_words = groups
            .checked_mul(3)
            .expect("position descriptor capacity")
            .div_ceil(4);
        let medium_prefix = medium_words
            .checked_add(flags)
            .expect("position descriptor capacity");
        let medium_size = size_with_prefix(medium_prefix, many);
        // Narrow only when at least two descriptor words disappear; coordinate
        // and arena bounds prove every encoded value fits without truncation.
        let narrow = !force_wide
            && groups >= 4
            && max_coordinate <= 0x7fff
            && narrow_size <= u16::MAX as usize;
        let medium = !force_wide
            && !narrow
            && groups >= 8
            && max_coordinate <= 0x7f_ffff
            && medium_size <= 0xff_ffff;
        let compact_prefix = groups
            .checked_add(flags)
            .expect("position descriptor capacity");
        let compact_size = size_with_prefix(compact_prefix, many);
        let wide = force_wide || (!narrow && !medium && u32::try_from(compact_size).is_err());
        let (prefix, size) = if wide {
            let prefix = groups.checked_mul(2).expect("position descriptor capacity");
            (
                prefix,
                size_with_prefix(
                    prefix,
                    many.checked_mul(2).expect("position lengths capacity"),
                ),
            )
        } else if narrow {
            (narrow_prefix, narrow_size)
        } else if medium {
            (medium_prefix, medium_size)
        } else {
            (compact_prefix, compact_size)
        };
        // The input Vec bounds groups below all three layout flag bits.
        debug_assert_eq!(groups & LAYOUT, 0);
        let data = if wide {
            Self::pack_words::<64>(sets, groups, prefix, size)
        } else if narrow {
            Self::pack_words::<16>(sets, groups, prefix, size)
        } else if medium {
            Self::pack_words::<24>(sets, groups, prefix, size)
        } else {
            Self::pack_words::<32>(sets, groups, prefix, size)
        };
        Self {
            groups: groups
                | if wide {
                    WIDE
                } else if narrow {
                    NARROW
                } else if medium {
                    MEDIUM
                } else {
                    0
                },
            data,
        }
    }

    // Dispatch once per block; each monomorph keeps width branches outside the
    // per-position loop, including large blocks which retain 24/32-bit storage.
    #[inline]
    fn pack_words<const BITS: usize>(
        sets: impl IntoIterator<Item = PositionSet>,
        groups: usize,
        prefix: usize,
        size: usize,
    ) -> Box<[u32]> {
        let flags_start = if BITS == 16 {
            groups.div_ceil(2)
        } else if BITS == 24 {
            (groups * 3).div_ceil(4)
        } else {
            groups
        };
        let mut data = Vec::with_capacity(size);
        data.resize(prefix, 0);
        for (index, set) in sets.into_iter().enumerate() {
            if BITS == 64 {
                let entry = match set {
                    PositionSet::Single(raw) => u64::from(raw),
                    PositionSet::Many(values) => {
                        let entry = INDIRECT | data.len() as u64;
                        let len = values.len() as u64;
                        data.push(len as u32);
                        data.push((len >> 32) as u32);
                        data.extend_from_slice(&values);
                        entry
                    }
                };
                data[index * 2] = entry as u32;
                data[index * 2 + 1] = (entry >> 32) as u32;
            } else {
                let entry = match set {
                    PositionSet::Single(raw) if BITS == 16 => {
                        (raw & 0x7fff) | ((raw >> 16) & 0x8000)
                    }
                    PositionSet::Single(raw) if BITS == 24 => {
                        (raw & 0x7f_ffff) | ((raw >> 8) & 0x80_0000)
                    }
                    PositionSet::Single(raw) => raw,
                    PositionSet::Many(values) => {
                        let offset = data.len() as u32;
                        data[flags_start + index / 32] |= 1 << (index % 32);
                        data.push(values.len() as u32);
                        data.extend_from_slice(&values);
                        offset
                    }
                };
                if BITS == 16 {
                    debug_assert!(entry <= u16::MAX as u32);
                    data[index / 2] |= entry << ((index % 2) * 16);
                } else if BITS == 24 {
                    debug_assert!(entry <= 0xff_ffff);
                    let bit = index * 24;
                    let word = bit / 32;
                    let shift = bit % 32;
                    data[word] |= entry << shift;
                    if shift > 8 {
                        data[word + 1] |= entry >> (32 - shift);
                    }
                } else {
                    data[index] = entry;
                }
            }
        }
        debug_assert_eq!(data.len(), size);
        data.into_boxed_slice()
    }

    pub(crate) fn len(&self) -> usize {
        if self.groups & INLINE != 0 {
            ((self.groups >> INLINE_ARITY_SHIFT) & 3) + 1
        } else {
            self.groups & !LAYOUT
        }
    }

    #[inline]
    pub(crate) fn get(&self, index: usize) -> Option<PositionSetRef<'_>> {
        let groups = self.len();
        if index >= groups {
            return None;
        }
        if self.groups & INLINE != 0 {
            return Some(PositionSetRef::Single(self.inline_raw(index, groups)));
        }
        self.get_allocated(index, groups)
    }

    // Keep the complete allocated layouts shared; only the bounded inline
    // prefix belongs in callers' per-transcript loops.
    #[inline(never)]
    fn get_allocated(&self, index: usize, groups: usize) -> Option<PositionSetRef<'_>> {
        let (start, len) = if self.groups & WIDE == 0 {
            let narrow = self.groups & NARROW != 0;
            let medium = self.groups & MEDIUM != 0;
            let (raw, flags_start) = if narrow {
                (
                    (self.data[index / 2] >> ((index % 2) * 16)) & 0xffff,
                    groups.div_ceil(2),
                )
            } else if medium {
                let bit = index * 24;
                let word = bit / 32;
                let shift = bit % 32;
                let mut raw = self.data[word] >> shift;
                if shift > 8 {
                    raw |= self.data[word + 1] << (32 - shift);
                }
                (raw & 0xff_ffff, (groups * 3).div_ceil(4))
            } else {
                (self.data[index], groups)
            };
            if self.groups & ALL_SINGLE != 0
                || self.data[flags_start + index / 32] & (1 << (index % 32)) == 0
            {
                let raw = if narrow {
                    (raw & 0x7fff) | ((raw & 0x8000) << 16)
                } else if medium {
                    (raw & 0x7f_ffff) | ((raw & 0x80_0000) << 8)
                } else {
                    raw
                };
                return Some(PositionSetRef::Single(raw));
            }
            let start = raw as usize;
            (start.checked_add(1)?, *self.data.get(start)? as usize)
        } else {
            let entry =
                u64::from(self.data[index * 2]) | (u64::from(self.data[index * 2 + 1]) << 32);
            if entry & INDIRECT == 0 {
                return Some(PositionSetRef::Single(entry as u32));
            }
            let start = usize::try_from(entry & !INDIRECT).ok()?;
            let len = u64::from(*self.data.get(start)?)
                | (u64::from(*self.data.get(start.checked_add(1)?)?) << 32);
            (start.checked_add(2)?, usize::try_from(len).ok()?)
        };
        Some(PositionSetRef::Many(
            self.data.get(start..start.checked_add(len)?)?,
        ))
    }

    // Inline descriptors already hold every sign bit. Prove a uniform
    // orientation directly from those bits, without a cached summary or any
    // additional work when constructing the index.
    #[inline]
    pub(crate) fn inline_uniform_strand(&self) -> Option<u8> {
        if self.groups & INLINE == 0 {
            return None;
        }
        let mask: u64 = match self.len() {
            1 => 1 << 31,
            2 => (1 << 23) | (1 << 47),
            3 => (1 << 15) | (1 << 31) | (1 << 47),
            _ => return None,
        };
        match self.groups as u64 & mask {
            0 => Some(1),
            signs if signs == mask => Some(0),
            _ => None,
        }
    }

    // Strand filtering consumes only signs. Decode singleton descriptors or
    // repeated-position endpoints without constructing coordinate views.
    #[inline]
    pub(crate) fn strand_summary_at(&self, index: usize) -> Option<u8> {
        let groups = self.len();
        if index >= groups {
            return None;
        }
        if self.groups & INLINE != 0 {
            return Some(u8::from(self.inline_raw(index, groups) & 0x8000_0000 == 0));
        }
        self.strand_summary_allocated(index, groups)
    }

    #[inline(never)]
    fn strand_summary_allocated(&self, index: usize, groups: usize) -> Option<u8> {
        let (start, len) = if self.groups & WIDE == 0 {
            let (raw, flags_start, sign_bit) = if self.groups & NARROW != 0 {
                (
                    (self.data[index / 2] >> ((index % 2) * 16)) & 0xffff,
                    groups.div_ceil(2),
                    0x8000,
                )
            } else if self.groups & MEDIUM != 0 {
                let bit = index * 24;
                let word = bit / 32;
                let shift = bit % 32;
                let mut raw = self.data[word] >> shift;
                if shift > 8 {
                    raw |= self.data[word + 1] << (32 - shift);
                }
                (raw & 0xff_ffff, (groups * 3).div_ceil(4), 0x80_0000)
            } else {
                (self.data[index], groups, 0x8000_0000)
            };
            if self.groups & ALL_SINGLE != 0
                || self.data[flags_start + index / 32] & (1 << (index % 32)) == 0
            {
                return Some(u8::from(raw & sign_bit == 0));
            }
            let start = raw as usize;
            (start.checked_add(1)?, *self.data.get(start)? as usize)
        } else {
            let entry =
                u64::from(self.data[index * 2]) | (u64::from(self.data[index * 2 + 1]) << 32);
            if entry & INDIRECT == 0 {
                return Some(u8::from(entry & 0x8000_0000 == 0));
            }
            let start = usize::try_from(entry & !INDIRECT).ok()?;
            let len = u64::from(*self.data.get(start)?)
                | (u64::from(*self.data.get(start.checked_add(1)?)?) << 32);
            (start.checked_add(2)?, usize::try_from(len).ok()?)
        };
        let first = *self.data.get(start)?;
        let last = *self.data.get(start.checked_add(len.checked_sub(1)?)?)?;
        Some(if (first ^ last) & 0x8000_0000 != 0 {
            2
        } else {
            u8::from(first & 0x8000_0000 == 0)
        })
    }

    #[cfg(test)]
    pub(crate) fn iter(&self) -> impl ExactSizeIterator<Item = PositionSetRef<'_>> {
        (0..self.len()).map(|i| self.get(i).expect("validated position entry"))
    }
}

impl PositionSetRef<'_> {
    pub(crate) fn strand_summary(&self) -> Option<u8> {
        let (min, max) = match self {
            Self::Single(raw) => (*raw, *raw),
            Self::Many(values) => (*values.first()?, *values.last()?),
        };
        let min_sense = min >> 31 == 0;
        let max_sense = max >> 31 == 0;
        Some(if min_sense == max_sense {
            u8::from(min_sense)
        } else {
            2
        })
    }
    pub(crate) fn contains(&self, value: u32) -> bool {
        match self {
            Self::Single(raw) => *raw == value,
            Self::Many(values) => values.binary_search(&value).is_ok(),
        }
    }
    pub(crate) fn any(&self, mut f: impl FnMut(u32) -> bool) -> bool {
        match self {
            Self::Single(raw) => f(*raw),
            Self::Many(values) => values.iter().copied().any(f),
        }
    }
    #[cfg(test)]
    pub(crate) fn for_each(&self, mut f: impl FnMut(u32)) {
        match self {
            Self::Single(raw) => f(*raw),
            Self::Many(values) => {
                for &raw in *values {
                    f(raw);
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn compact_and_wide_position_groups_preserve_every_coordinate_and_strand_bit() {
        for groups in [
            vec![],
            vec![vec![]],
            vec![vec![u32::MAX]],
            vec![vec![0]],
            vec![vec![0x7fff_ffff]],
            vec![vec![0x8000_0000], vec![0x7fff_ffff]],
            vec![vec![0x007f_ffff], vec![0x807f_ffff]],
            vec![vec![0x0080_0000], vec![0x8080_0000]],
            vec![vec![0x0000_7fff], vec![0x8000_7fff], vec![0x8000_0000]],
            vec![vec![0x0000_8000], vec![0x8000_8000], vec![0]],
            vec![vec![0x8000_0001], vec![1], vec![0x100_0000]],
            vec![vec![0], vec![0x7f_ffff]],
            vec![vec![0x8000_0000], vec![0x807f_ffff]],
            vec![vec![0], vec![1], vec![0x7fff]],
            vec![vec![0x8000_0000], vec![0x8000_0001], vec![0x8000_7fff]],
            // Width boundaries retain the strand bit and all larger coordinates.
            vec![vec![0x7fff], vec![0x8000_7fff], vec![0], vec![0x8000_0000]],
            vec![vec![0x8000], vec![0x8000_8000], vec![0], vec![1]],
            vec![
                vec![0x7f_ffff],
                vec![0x807f_ffff],
                vec![0],
                vec![1],
                vec![2],
                vec![3],
                vec![4],
                vec![5],
            ],
            vec![
                vec![0x80_0000],
                vec![0x8080_0000],
                vec![0],
                vec![1],
                vec![2],
                vec![3],
                vec![4],
                vec![5],
            ],
            vec![
                (0..65_536).collect(),
                vec![0],
                vec![1],
                vec![2],
                vec![3],
                vec![4],
                vec![5],
                vec![6],
            ],
            (0..65)
                .map(|i| {
                    if i % 3 == 0 || matches!(i, 31 | 32 | 63 | 64) {
                        vec![i, 0x8000_0000 | i]
                    } else {
                        vec![i]
                    }
                })
                .collect(),
            vec![vec![], vec![u32::MAX], vec![0, 0x8000_0000, u32::MAX]],
            vec![
                vec![0],
                vec![u32::MAX],
                vec![0x8000_0000],
                vec![0x7fff_ffff],
            ],
            vec![
                vec![u32::MAX, 0, 0x8000_0000, 0, u32::MAX],
                vec![17, 17],
                vec![123],
            ],
        ] {
            let expected: Vec<Vec<u32>> = groups
                .iter()
                .map(|values| {
                    let mut values = values.clone();
                    values.sort_unstable();
                    values.dedup();
                    values
                })
                .collect();
            for variant in 0..3 {
                if variant == 2 && !expected.iter().all(|values| values.len() == 1) {
                    continue;
                }
                let packed = if variant == 2 {
                    PositionSets::from_singletons(expected.iter().map(|values| values[0]))
                } else {
                    PositionSets::pack(
                        groups
                            .iter()
                            .cloned()
                            .map(PositionSet::from_values)
                            .collect(),
                        variant == 1,
                    )
                };
                if INLINE != 0 && variant != 1 && expected.len() == 1 && groups[0].len() == 1 {
                    assert!(packed.data.is_empty(), "inline singleton allocated storage");
                }
                if variant == 2 && expected.len() >= 4 && ALL_SINGLE != 0 {
                    assert_ne!(packed.groups & ALL_SINGLE, 0);
                    if packed.groups & WIDE == 0 {
                        let words = if packed.groups & NARROW != 0 {
                            expected.len().div_ceil(2)
                        } else if packed.groups & MEDIUM != 0 {
                            (expected.len() * 3).div_ceil(4)
                        } else {
                            expected.len()
                        };
                        assert_eq!(
                            packed.data.len(),
                            words,
                            "singleton arena retained flag words"
                        );
                    }
                }
                let uniform = if packed.groups & INLINE != 0 {
                    let first = expected[0][0] >> 31;
                    expected
                        .iter()
                        .all(|values| values.len() == 1 && values[0] >> 31 == first)
                        .then_some(u8::from(first == 0))
                } else {
                    None
                };
                assert_eq!(packed.inline_uniform_strand(), uniform);
                assert_eq!(packed.len(), expected.len());
                assert!(packed.get(expected.len()).is_none());
                assert!(packed.strand_summary_at(expected.len()).is_none());
                for (index, (set, expected)) in packed.iter().zip(&expected).enumerate() {
                    let mut actual = Vec::new();
                    set.for_each(|value| actual.push(value));
                    assert_eq!(&actual, expected);
                    let summary = expected.first().zip(expected.last()).map(|(&min, &max)| {
                        if min >> 31 == max >> 31 {
                            u8::from(min >> 31 == 0)
                        } else {
                            2
                        }
                    });
                    assert_eq!(set.strand_summary(), summary);
                    assert_eq!(packed.strand_summary_at(index), summary);
                    for query in [0, 1, 17, 123, 0x7fff_ffff, 0x8000_0000, u32::MAX] {
                        assert_eq!(set.contains(query), expected.binary_search(&query).is_ok());
                        assert_eq!(set.any(|raw| raw == query), expected.contains(&query));
                    }
                }
            }
        }
    }
}
