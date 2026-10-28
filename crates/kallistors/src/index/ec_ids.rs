//! Immutable transcript IDs with allocation-free native one/two-ID classes.
//!
//! On 64-bit targets the Box slice's unused null-pointer niche holds the inline
//! variant without enlarging the 16-byte value. Narrow targets keep boxed
//! storage. All ownership, moves and drops use ordinary safe Rust enum fields.
use std::ops::Deref;

#[derive(Clone)]
pub(crate) enum EcIds {
    #[cfg(target_pointer_width = "64")]
    // Equal entries encode a singleton. A real duplicate pair stays on the
    // heap, preserving its multiplicity instead of collapsing it to one ID.
    Inline([u32; 2]),
    Heap(Box<[u32]>),
}
impl EcIds {
    // Nonpositional readers consume their boxes into Vecs; preserve the
    // original allocation there instead of creating a needless round trip.
    pub(super) fn heap(values: Box<[u32]>) -> Self {
        Self::Heap(values)
    }
    pub(super) fn into_vec(self) -> Vec<u32> {
        match self {
            #[cfg(target_pointer_width = "64")]
            Self::Inline(values) if values[0] == values[1] => vec![values[0]],
            #[cfg(target_pointer_width = "64")]
            Self::Inline(values) => values.to_vec(),
            Self::Heap(values) => values.into_vec(),
        }
    }
    pub(super) fn sort_unstable(&mut self) {
        match self {
            #[cfg(target_pointer_width = "64")]
            Self::Inline(values) if values[0] != values[1] => values.sort_unstable(),
            #[cfg(target_pointer_width = "64")]
            Self::Inline(_) => {}
            Self::Heap(values) => values.sort_unstable(),
        }
    }
}
impl From<Box<[u32]>> for EcIds {
    fn from(values: Box<[u32]>) -> Self {
        #[cfg(target_pointer_width = "64")]
        match *values.as_ref() {
            [value] => return Self::Inline([value, value]),
            [a, b] if a != b => return Self::Inline([a, b]),
            _ => {}
        }
        Self::Heap(values)
    }
}
impl From<Vec<u32>> for EcIds {
    fn from(values: Vec<u32>) -> Self {
        #[cfg(target_pointer_width = "64")]
        match *values.as_slice() {
            [value] => return Self::Inline([value, value]),
            [a, b] if a != b => return Self::Inline([a, b]),
            _ => {}
        }
        Self::Heap(values.into_boxed_slice())
    }
}
impl Deref for EcIds {
    type Target = [u32];
    fn deref(&self) -> &[u32] {
        match self {
            #[cfg(target_pointer_width = "64")]
            Self::Inline(values) if values[0] == values[1] => &values[..1],
            #[cfg(target_pointer_width = "64")]
            Self::Inline(values) => values,
            Self::Heap(values) => values,
        }
    }
}
impl AsRef<[u32]> for EcIds {
    fn as_ref(&self) -> &[u32] {
        self
    }
}
impl PartialEq for EcIds {
    fn eq(&self, other: &Self) -> bool {
        self.as_ref() == other.as_ref()
    }
}
impl Eq for EcIds {}
impl std::fmt::Debug for EcIds {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_list().entries(self.iter()).finish()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn inline_ec_ids_preserve_width_multisets_moves_sort_and_box_footprint() {
        assert_eq!(
            std::mem::size_of::<EcIds>(),
            std::mem::size_of::<Box<[u32]>>()
        );
        assert_eq!(
            std::mem::align_of::<EcIds>(),
            std::mem::align_of::<Box<[u32]>>()
        );
        type OriginalShape = (u32, u32, Box<[u32]>, Option<super::super::PositionSets>);
        assert_eq!(
            std::mem::size_of::<super::super::EcBlock>(),
            std::mem::size_of::<OriginalShape>()
        );
        let mut inputs = vec![
            vec![],
            vec![0],
            vec![u32::MAX],
            vec![0, u32::MAX],
            vec![17, 3],
            vec![7, 7],
            vec![1, 3, 5],
        ];
        for len in [0usize, 1, 2, 3, 4, 7, 31, 32, 33, 257] {
            for salt in [0u32, 7, u32::MAX] {
                inputs.push(
                    (0..len)
                        .map(|i| (i as u32).wrapping_mul(17).wrapping_add(salt))
                        .collect(),
                );
            }
        }
        let mut moved = Vec::new();
        for values in &inputs {
            let ids = EcIds::from(values.clone());
            assert_eq!(ids.as_ref(), values);
            // In particular duplicate pairs remain two IDs, not a singleton.
            assert_eq!(ids.len(), values.len());
            let boxed = EcIds::heap(values.clone().into_boxed_slice());
            assert_eq!(ids, boxed);
            assert_eq!(boxed.into_vec(), *values);
            assert_eq!(ids.clone().into_vec(), *values);
            let mut sorted = ids.clone();
            let mut expected = values.clone();
            sorted.sort_unstable();
            expected.sort_unstable();
            assert_eq!(sorted.as_ref(), expected);
            moved.push(ids);
        }
        moved.shrink_to_fit();
        let moved = moved.into_boxed_slice();
        for (ids, values) in moved.iter().zip(&inputs) {
            assert_eq!(ids.as_ref(), values);
        }
        fn require_send_sync<T: Send + Sync>() {}
        require_send_sync::<EcIds>();
    }
}
