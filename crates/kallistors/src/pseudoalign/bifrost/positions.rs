//! Lossless delta coding of temporary minimizer positions during loading.

pub(super) enum BitmapPositions {
    Raw(Vec<u32>),
    Deltas {
        values: Vec<u8>,
        len: usize,
        last: Option<u32>,
    },
}
impl BitmapPositions {
    pub(super) fn new(values: Vec<u32>, compact: bool) -> Self {
        if !compact {
            return Self::Raw(values);
        }
        let first_escape = usize::from(values.first().is_some_and(|&value| value >= 255));
        let (escapes, unsorted) =
            values
                .windows(2)
                .fold((first_escape, false), |(escapes, unsorted), pair| {
                    (
                        escapes + usize::from(pair[1].wrapping_sub(pair[0]) >= 255),
                        unsorted | (pair[0] > pair[1]),
                    )
                });
        let encoded_len = values.len() + escapes * 4;
        if unsorted || encoded_len > values.len() * std::mem::size_of::<u32>() {
            return Self::Raw(values);
        }
        let mut previous = 0;
        let mut deltas = Vec::<u8>::with_capacity(encoded_len);
        let mut written = 0;
        for &value in &values {
            let gap = value - previous;
            // SAFETY: the preceding pass counted every escape and reserved
            // exactly one or five bytes per value. Writes cover disjoint
            // slots inside that capacity; the input allocation is distinct.
            unsafe {
                let destination = deltas.as_mut_ptr().add(written);
                if gap < 255 {
                    destination.write(gap as u8);
                    written += 1;
                } else {
                    destination.write(255);
                    std::ptr::copy_nonoverlapping(
                        gap.to_le_bytes().as_ptr(),
                        destination.add(1),
                        4,
                    );
                    written += 5;
                }
            }
            previous = value;
        }
        debug_assert_eq!(written, encoded_len);
        // SAFETY: each byte of this prefix was initialized above.
        unsafe {
            deltas.set_len(written);
        }
        let len = values.len();
        let last = values.last().copied();
        crate::index::scratch_pages::retire(values);
        Self::Deltas {
            values: deltas,
            len,
            last,
        }
    }
    pub(super) fn retire(self) {
        match self {
            Self::Raw(values) => crate::index::scratch_pages::retire(values),
            Self::Deltas { values, .. } => crate::index::scratch_pages::retire(values),
        }
    }
    pub(super) fn len(&self) -> usize {
        match self {
            Self::Raw(v) => v.len(),
            Self::Deltas { len, .. } => *len,
        }
    }
    pub(super) fn last(&self) -> Option<u32> {
        match self {
            Self::Raw(v) => v.last().copied(),
            Self::Deltas { last, .. } => *last,
        }
    }
    pub(super) fn iter(&self) -> PositionIter<'_> {
        match self {
            Self::Raw(v) => PositionIter::Raw(v.iter()),
            Self::Deltas { values, .. } => PositionIter::Deltas {
                iter: values.iter(),
                previous: 0,
            },
        }
    }
}
pub(super) enum PositionIter<'a> {
    Raw(std::slice::Iter<'a, u32>),
    Deltas {
        iter: std::slice::Iter<'a, u8>,
        previous: u32,
    },
}
impl Iterator for PositionIter<'_> {
    type Item = u32;
    #[inline]
    fn next(&mut self) -> Option<u32> {
        match self {
            Self::Raw(iter) => iter.next().copied(),
            Self::Deltas { iter, previous } => {
                let byte = *iter.next()?;
                let gap = if byte == 255 {
                    let mut bytes = [0; 4];
                    for byte in &mut bytes {
                        *byte = *iter.next().expect("validated position escape");
                    }
                    u32::from_le_bytes(bytes)
                } else {
                    u32::from(byte)
                };
                *previous += gap;
                Some(*previous)
            }
        }
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn delta_positions_preserve_zero_bounds_large_gaps_and_unsorted_fallback() {
        for original in [
            vec![],
            vec![0],
            vec![0, 1, 65535, 65536, 131071],
            vec![0, 65536],
            vec![u32::MAX],
            vec![2, 1],
            vec![0, 0, 1],
            vec![0, 254, 509, 765, 66300, u32::MAX],
        ] {
            for compact in [false, true] {
                let got = BitmapPositions::new(original.clone(), compact);
                assert_eq!(got.len(), original.len());
                assert_eq!(got.last(), original.last().copied());
                assert_eq!(got.iter().collect::<Vec<_>>(), original);
            }
        }
    }
}
