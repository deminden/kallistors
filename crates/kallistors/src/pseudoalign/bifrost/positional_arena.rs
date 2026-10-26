//! Ordered positional blocks in bounded contiguous pages, without one box per node.
use crate::Result;
use crate::index::EcBlock;
use std::io::Cursor;
use std::ops::Index;
const PAGE: usize = 1 << 16;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Span {
    chunk: usize,
    start: usize,
    len: usize,
}
impl Span {
    fn packed(self) -> Option<u64> {
        let chunk = u16::try_from(self.chunk).ok()?;
        let start = u16::try_from(self.start).ok()?;
        let len = u32::try_from(self.len).ok()?;
        Some((u64::from(chunk) << 48) | (u64::from(start) << 32) | u64::from(len))
    }
    fn unpack(value: u64) -> Self {
        Self {
            chunk: (value >> 48) as usize,
            start: ((value >> 32) & 0xffff) as usize,
            len: (value & 0xffff_ffff) as usize,
        }
    }
}
const EXTRA: u32 = 1 << 31;
enum Directory {
    Small { words: Vec<u32>, extras: Vec<u64> },
    Compact(Vec<u64>),
    Wide(Vec<Span>),
}
impl Directory {
    fn new(capacity: usize) -> Self {
        if capacity >= 4096 {
            Self::small(capacity)
        } else {
            Self::Compact(Vec::with_capacity(capacity))
        }
    }
    fn small(capacity: usize) -> Self {
        Self::Small {
            words: Vec::with_capacity(capacity),
            extras: Vec::new(),
        }
    }
    fn len(&self) -> usize {
        match self {
            Self::Small { words, .. } => words.len(),
            Self::Compact(v) => v.len(),
            Self::Wide(v) => v.len(),
        }
    }
    #[inline]
    fn decode_small(word: u32, extras: &[u64]) -> Span {
        if word == 0 {
            Span {
                chunk: 0,
                start: 0,
                len: 0,
            }
        } else if word & EXTRA == 0 {
            let global = (word - 1) as usize;
            Span {
                chunk: global >> 16,
                start: global & 0xffff,
                len: 1,
            }
        } else {
            Span::unpack(extras[(word & !EXTRA) as usize])
        }
    }
    fn push(&mut self, span: Span) {
        match self {
            Self::Small { words, extras } => {
                let packed = span.packed();
                let direct = if span
                    == (Span {
                        chunk: 0,
                        start: 0,
                        len: 0,
                    }) {
                    Some(0)
                } else if span.len == 1 && packed.is_some() {
                    let global = ((span.chunk as u64) << 16) | span.start as u64;
                    u32::try_from(global + 1).ok().filter(|&word| word < EXTRA)
                } else {
                    None
                };
                if let Some(word) = direct {
                    words.push(word);
                    return;
                }
                if let Some(value) = packed {
                    if extras.len() < (EXTRA as usize) {
                        words.push(EXTRA | extras.len() as u32);
                        extras.push(value);
                        // Sparse exceptions save space. If they become dense,
                        // restore the old8-byte representation before filling
                        // the reference-sized directory. Preserve its capacity.
                        if words.len() >= 4096 && extras.len() > words.len() / 2 {
                            let mut full = Vec::with_capacity(words.capacity());
                            full.extend(words.iter().map(|&word| {
                                Self::decode_small(word, extras)
                                    .packed()
                                    .expect("representable stored span")
                            }));
                            crate::index::scratch_pages::retire(std::mem::take(words));
                            crate::index::scratch_pages::retire(std::mem::take(extras));
                            *self = Self::Compact(full);
                        }
                        return;
                    }
                    let mut full = Vec::with_capacity(words.capacity().max(words.len() + 1));
                    full.extend(words.iter().map(|&word| {
                        Self::decode_small(word, extras)
                            .packed()
                            .expect("representable stored span")
                    }));
                    full.push(value);
                    crate::index::scratch_pages::retire(std::mem::take(words));
                    crate::index::scratch_pages::retire(std::mem::take(extras));
                    *self = Self::Compact(full);
                } else {
                    let mut wide = Vec::with_capacity(words.capacity().max(words.len() + 1));
                    wide.extend(words.iter().map(|&word| Self::decode_small(word, extras)));
                    wide.push(span);
                    crate::index::scratch_pages::retire(std::mem::take(words));
                    crate::index::scratch_pages::retire(std::mem::take(extras));
                    *self = Self::Wide(wide);
                }
            }
            Self::Compact(v) => {
                if let Some(word) = span.packed() {
                    v.push(word);
                    return;
                }
                let mut wide: Vec<_> = v.iter().copied().map(Span::unpack).collect();
                wide.push(span);
                crate::index::scratch_pages::retire(std::mem::take(v));
                *self = Self::Wide(wide);
            }
            Self::Wide(v) => v.push(span),
        }
    }
    #[inline]
    fn get(&self, index: usize) -> Option<Span> {
        match self {
            Self::Small { words, extras } => words
                .get(index)
                .map(|&word| Self::decode_small(word, extras)),
            Self::Compact(v) => v.get(index).copied().map(Span::unpack),
            Self::Wide(v) => v.get(index).copied(),
        }
    }
    fn retire(self) {
        match self {
            Self::Small { words, extras } => {
                crate::index::scratch_pages::retire(words);
                crate::index::scratch_pages::retire(extras);
            }
            Self::Compact(v) => crate::index::scratch_pages::retire(v),
            Self::Wide(v) => crate::index::scratch_pages::retire(v),
        }
    }
}

pub(crate) struct OwnedEcBlocks {
    rows: Directory,
    chunks: Vec<Vec<EcBlock>>,
    len: usize,
}
impl OwnedEcBlocks {
    pub(super) fn empty() -> Self {
        Self {
            rows: Directory::new(0),
            chunks: Vec::new(),
            len: 0,
        }
    }
    #[inline]
    pub(crate) fn len(&self) -> usize {
        self.len
    }
    #[inline]
    pub(crate) fn is_empty(&self) -> bool {
        self.len == 0
    }
    #[inline]
    pub(crate) fn get(&self, index: usize) -> Option<&[EcBlock]> {
        let span = self.rows.get(index)?;
        if span.len == 0 {
            return Some(&[]);
        }
        self.chunks
            .get(span.chunk)?
            .get(span.start..span.start.checked_add(span.len)?)
    }
    pub(crate) fn iter(&self) -> OwnedIter<'_> {
        OwnedIter {
            owner: self,
            front: 0,
            back: self.len,
        }
    }
    pub(super) fn retire(self) {
        self.rows.retire();
        // Clearing this outer Vec drops every block/page first. Retire once,
        // rather than requesting page reclamation separately for every page.
        crate::index::scratch_pages::retire(self.chunks);
    }
}
impl Index<usize> for OwnedEcBlocks {
    type Output = [EcBlock];
    #[inline]
    fn index(&self, index: usize) -> &Self::Output {
        self.get(index).expect("valid unitig EC blocks")
    }
}
pub(crate) struct OwnedIter<'a> {
    owner: &'a OwnedEcBlocks,
    front: usize,
    back: usize,
}
impl<'a> Iterator for OwnedIter<'a> {
    type Item = &'a [EcBlock];
    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        if self.front == self.back {
            return None;
        }
        let index = self.front;
        self.front += 1;
        self.owner.get(index)
    }
    fn size_hint(&self) -> (usize, Option<usize>) {
        let n = self.back - self.front;
        (n, Some(n))
    }
}
impl DoubleEndedIterator for OwnedIter<'_> {
    fn next_back(&mut self) -> Option<Self::Item> {
        if self.front == self.back {
            return None;
        }
        self.back -= 1;
        self.owner.get(self.back)
    }
}
impl ExactSizeIterator for OwnedIter<'_> {}
impl<'a> IntoIterator for &'a OwnedEcBlocks {
    type Item = &'a [EcBlock];
    type IntoIter = OwnedIter<'a>;
    fn into_iter(self) -> Self::IntoIter {
        self.iter()
    }
}
impl FromIterator<Box<[EcBlock]>> for OwnedEcBlocks {
    fn from_iter<T: IntoIterator<Item = Box<[EcBlock]>>>(iter: T) -> Self {
        let iter = iter.into_iter();
        let mut builder = EcBlocksBuilder::new(iter.size_hint().0);
        for node in iter {
            builder.push(node.into_vec());
        }
        builder.finish()
    }
}

pub(super) struct EcBlocksBuilder {
    rows: Directory,
    chunks: Vec<Vec<EcBlock>>,
    current: Vec<EcBlock>,
    page: usize,
}
impl EcBlocksBuilder {
    pub(super) fn new(nodes: usize) -> Self {
        Self::with_page(nodes, nodes.clamp(1, PAGE))
    }
    fn with_page(nodes: usize, page: usize) -> Self {
        assert!((1..=PAGE).contains(&page));
        Self {
            rows: Directory::new(nodes),
            chunks: Vec::new(),
            current: Vec::new(),
            page,
        }
    }
    fn flush(&mut self) {
        if !self.current.is_empty() {
            self.chunks.push(std::mem::take(&mut self.current));
        }
    }
    fn begin(&mut self, count: usize) -> Span {
        if count == 0 {
            return Span {
                chunk: 0,
                start: 0,
                len: 0,
            };
        }
        if count > self.current.capacity() - self.current.len() {
            self.flush();
        }
        if self.current.capacity() == 0 {
            self.current = Vec::with_capacity(self.page.max(count));
        }
        Span {
            chunk: self.chunks.len(),
            start: self.current.len(),
            len: count,
        }
    }
    fn commit(&mut self, span: Span) {
        debug_assert!(span.len == 0 || self.current.len() == span.start + span.len);
        self.rows.push(span);
        if self.current.capacity() > self.page {
            self.flush();
        }
    }
    pub(super) fn push(&mut self, node: Vec<EcBlock>) {
        let span = self.begin(node.len());
        self.current.extend(node);
        self.commit(span);
    }
    pub(super) fn append_payload(&mut self, bytes: &[u8], positional: bool) -> Result<()> {
        let mut reader = Cursor::new(bytes);
        let _node_id = super::read_u32_le(&mut reader)?;
        if !positional {
            self.push(crate::index::read_block_array_blocks(&mut reader)?);
            return Ok(());
        }
        let count = crate::index::read_positional_block_count(&mut reader)?;
        let Ok(count_usize) = usize::try_from(count) else {
            // Preserve the original Vec/cursor decoder on narrow targets if
            // the serialized count is not addressable; it retains its errors.
            let mut blocks = Vec::with_capacity(count as usize);
            crate::index::read_positional_blocks_into(&mut reader, count, |block| {
                blocks.push(block)
            })?;
            self.push(blocks);
            return Ok(());
        };
        let span = self.begin(count_usize);
        crate::index::read_positional_blocks_into(&mut reader, count, |block| {
            self.current.push(block)
        })?;
        self.commit(span);
        Ok(())
    }

    pub(super) fn finish_for_transfer(mut self) -> OwnedEcBlocks {
        self.flush();
        // Task-sized pages must not accumulate their spare tails across the
        // complete reference. Shrink every page before its ownership moves.
        for page in &mut self.chunks {
            *page = std::mem::take(page).into_boxed_slice().into_vec();
        }
        self.finish()
    }
    pub(super) fn append_finished(&mut self, other: OwnedEcBlocks) -> Result<()> {
        self.flush();
        let OwnedEcBlocks {
            rows,
            mut chunks,
            len,
        } = other;
        let base = self.chunks.len();
        for row in 0..len {
            let mut span = rows.get(row).expect("complete transferred node directory");
            if span.len != 0 {
                span.chunk = span.chunk.checked_add(base).ok_or_else(|| {
                    crate::Error::InvalidFormat("positional page directory overflow".into())
                })?;
            }
            self.rows.push(span);
        }
        self.chunks.append(&mut chunks);
        rows.retire();
        Ok(())
    }
    pub(super) fn finish(mut self) -> OwnedEcBlocks {
        self.flush();
        let len = self.rows.len();
        OwnedEcBlocks {
            rows: self.rows,
            chunks: self.chunks,
            len,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn node(id: u32, count: usize) -> Vec<EcBlock> {
        (0..count)
            .map(|i| EcBlock {
                lb: i as u32,
                ub: i as u32 + 7,
                ec: vec![id, u32::MAX].into_boxed_slice().into(),
                positions: None,
            })
            .collect()
    }
    #[test]
    fn positional_pages_preserve_empty_crossings_oversized_nodes_and_wide_directory() {
        for page in [1usize, 2, 3, 7, 16] {
            let lengths = [0, 1, 2, 0, 7, 3, 1, 18, 0, 2, 4, 1];
            let mut builder = EcBlocksBuilder::with_page(lengths.len(), page);
            for (id, &len) in lengths.iter().enumerate() {
                builder.push(node(id as u32, len));
            }
            let blocks = builder.finish();
            assert_eq!(blocks.len(), lengths.len());
            assert!(!blocks.is_empty());
            assert!(blocks.get(lengths.len()).is_none());
            for (id, (actual, &len)) in blocks.iter().zip(&lengths).enumerate() {
                assert_eq!(actual.len(), len);
                for (i, block) in actual.iter().enumerate() {
                    assert_eq!((block.lb, block.ub), (i as u32, i as u32 + 7));
                    assert_eq!(block.ec.as_ref(), [id as u32, u32::MAX]);
                    assert!(block.positions.is_none());
                }
            }
            assert_eq!(
                blocks
                    .iter()
                    .rev()
                    .map(<[EcBlock]>::len)
                    .collect::<Vec<_>>(),
                lengths.iter().rev().copied().collect::<Vec<_>>()
            );
        }
        let mut cases = vec![
            Span {
                chunk: 0,
                start: 0,
                len: 0,
            },
            Span {
                chunk: 65535,
                start: 65535,
                len: u32::MAX as usize,
            },
            Span {
                chunk: 65536,
                start: 1,
                len: 2,
            },
            Span {
                chunk: 1,
                start: 65536,
                len: 2,
            },
        ];
        if let Some(len) = (u32::MAX as usize).checked_add(1) {
            cases.push(Span {
                chunk: 1,
                start: 1,
                len,
            });
        }
        let mut dir = Directory::new(0);
        for &span in &cases {
            dir.push(span);
        }
        assert!(matches!(dir, Directory::Wide(_)));
        for (i, &span) in cases.iter().enumerate() {
            assert_eq!(dir.get(i), Some(span));
        }
    }
    #[test]
    fn positional_page_stream_preserves_positions_and_truncation_errors() {
        let bitmap = |words: &[u32]| {
            let mut bytes = vec![1];
            bytes.extend_from_slice(&(words.len() as u32).to_le_bytes());
            for word in words {
                bytes.extend_from_slice(&word.to_le_bytes());
            }
            bytes
        };
        let ids = bitmap(&[1, 7]);
        let mut bytes = 9u32.to_le_bytes().to_vec();
        bytes.push(1);
        bytes.extend_from_slice(&2u32.to_le_bytes());
        bytes.extend_from_slice(&31u32.to_le_bytes());
        bytes.extend_from_slice(&(ids.len() as u64).to_le_bytes());
        bytes.extend_from_slice(&ids);
        bytes.extend_from_slice(&2u64.to_le_bytes());
        for set in [&[0x8000_0003][..], &[4, 9][..]] {
            let values = bitmap(set);
            bytes.extend_from_slice(&(values.len() as u64).to_le_bytes());
            bytes.extend(values);
        }
        let mut multiple = 9u32.to_le_bytes().to_vec();
        multiple.push(2);
        multiple.extend_from_slice(&3u64.to_le_bytes());
        for _ in 0..3 {
            multiple.extend_from_slice(&bytes[5..]);
        }
        let mut empty = 9u32.to_le_bytes().to_vec();
        empty.push(0);
        for (frame, count) in [(empty, 0usize), (bytes, 1), (multiple, 3)] {
            for positional in [false, true] {
                let mut builder = EcBlocksBuilder::with_page(1, 2);
                builder.append_payload(&frame, positional).unwrap();
                let data = builder.finish();
                assert_eq!(data.len(), 1);
                assert_eq!(data[0].len(), count);
                for block in &data[0] {
                    assert_eq!((block.lb, block.ub), (2, 31));
                    assert_eq!(block.ec.as_ref(), [1, 7]);
                    assert_eq!(block.positions.is_some(), positional);
                    if let Some(positions) = &block.positions {
                        let values: Vec<Vec<u32>> = positions
                            .iter()
                            .map(|p| match p {
                                crate::index::PositionSetRef::Single(v) => vec![v],
                                crate::index::PositionSetRef::Many(v) => v.to_vec(),
                            })
                            .collect();
                        assert_eq!(values, [vec![0x8000_0003], vec![4, 9]]);
                    }
                }
                let mut padded = frame.clone();
                padded.extend_from_slice(&[9, 2, 8]);
                let mut builder = EcBlocksBuilder::with_page(1, 2);
                builder.append_payload(&padded, positional).unwrap();
                assert_eq!(builder.finish()[0].len(), count);
                for end in 0..frame.len() {
                    let mut builder = EcBlocksBuilder::with_page(1, 2);
                    assert!(builder.append_payload(&frame[..end], positional).is_err());
                }
            }
        }
    }
    #[test]
    fn narrow_positional_directory_preserves_tag_edges_and_adapts_dense_exceptions() {
        let cases = [
            Span {
                chunk: 0,
                start: 0,
                len: 0,
            },
            Span {
                chunk: 0,
                start: 0,
                len: 1,
            },
            Span {
                chunk: 32767,
                start: 65534,
                len: 1,
            },
            Span {
                chunk: 32767,
                start: 65535,
                len: 1,
            },
            Span {
                chunk: 65535,
                start: 65535,
                len: 1,
            },
            Span {
                chunk: 7,
                start: 9,
                len: 2,
            },
            Span {
                chunk: 65535,
                start: 65535,
                len: u32::MAX as usize,
            },
            Span {
                chunk: 1,
                start: 2,
                len: 0,
            },
        ];
        let mut dir = Directory::small(cases.len());
        for &span in &cases {
            dir.push(span);
        }
        assert!(matches!(dir, Directory::Small { .. }));
        for (i, &span) in cases.iter().enumerate() {
            assert_eq!(dir.get(i), Some(span));
        }
        assert!(dir.get(cases.len()).is_none());
        dir.push(Span {
            chunk: 65536,
            start: 1,
            len: 3,
        });
        assert!(matches!(dir, Directory::Wide(_)));
        for (i, &span) in cases.iter().enumerate() {
            assert_eq!(dir.get(i), Some(span));
        }
        let mut dense = Directory::new(5000);
        for i in 0..5000 {
            dense.push(Span {
                chunk: 0,
                start: i,
                len: 2,
            });
        }
        assert!(matches!(dense, Directory::Compact(_)));
        for i in 0..5000 {
            assert_eq!(
                dense.get(i),
                Some(Span {
                    chunk: 0,
                    start: i,
                    len: 2
                })
            );
        }
        let mut sparse = Directory::new(5000);
        for i in 0..5000 {
            sparse.push(Span {
                chunk: 0,
                start: i,
                len: if i % 100 == 0 { 2 } else { 1 },
            });
        }
        if let Directory::Small { words, extras } = &sparse {
            assert_eq!(words.len(), 5000);
            assert_eq!(extras.len(), 50);
            assert!(words.capacity() * 4 + extras.capacity() * 8 <= 5000 * 8);
        } else {
            panic!("sparse directory should retain compact words");
        }
        for i in 0..5000 {
            assert_eq!(
                sparse.get(i),
                Some(Span {
                    chunk: 0,
                    start: i,
                    len: if i % 100 == 0 { 2 } else { 1 }
                })
            );
        }
    }
    #[test]
    fn positional_page_transfer_preserves_wide_rebased_spans_and_tight_capacity() {
        let mut source = EcBlocksBuilder::with_page(4, 2);
        for (id, count) in [0, 1, 7, 2].into_iter().enumerate() {
            source.push(node(id as u32, count));
        }
        let source = source.finish_for_transfer();
        for page in &source.chunks {
            assert_eq!(page.len(), page.capacity());
        }
        let mut target = EcBlocksBuilder::new(4);
        // Force the complete usize-wide chunk fallback without a huge graph.
        target.chunks.resize_with(65536, Vec::new);
        target.append_finished(source).unwrap();
        let owned = target.finish();
        assert!(matches!(owned.rows, Directory::Wide(_)));
        assert_eq!(owned.len(), 4);
        for (id, count) in [0, 1, 7, 2].into_iter().enumerate() {
            let blocks = &owned[id];
            assert_eq!(blocks.len(), count);
            for (b, block) in blocks.iter().enumerate() {
                assert_eq!(&*block.ec, &[id as u32, u32::MAX]);
                assert_eq!((block.lb, block.ub), (b as u32, b as u32 + 7));
            }
        }
        assert!(owned.get(4).is_none());
    }
}
