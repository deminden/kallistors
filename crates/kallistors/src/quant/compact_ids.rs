//! Narrow dense transcript IDs in their existing allocation.

pub(super) struct CompactTranscriptIds {
    words: Vec<usize>,
    len: usize,
}
impl CompactTranscriptIds {
    pub(super) fn new(mut words: Vec<usize>) -> Self {
        let len = words.len();
        for i in 0..len {
            let value = u32::try_from(words[i]).expect("validated dense transcript ID");
            // SAFETY: usize is at least 32 bits on supported targets. Each
            // four-byte write stays in initialized storage and behind the
            // next unread usize. Alignment is inherited from Vec<usize>.
            unsafe {
                words.as_mut_ptr().cast::<u32>().add(i).write(value);
            }
        }
        words.truncate((len * std::mem::size_of::<u32>()).div_ceil(std::mem::size_of::<usize>()));
        crate::index::scratch_pages::discard_unused(&mut words, len);
        Self { words, len }
    }
}
impl std::ops::Deref for CompactTranscriptIds {
    type Target = [u32];
    fn deref(&self) -> &Self::Target {
        // SAFETY: construction initialized exactly len native u32 values in
        // this aligned prefix. The original allocation still owns the bytes;
        // immutable borrows cannot outlive it or mutate its backing words.
        unsafe { std::slice::from_raw_parts(self.words.as_ptr().cast::<u32>(), self.len) }
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn narrowing_preserves_empty_odd_even_and_page_boundaries() {
        for len in [0, 1, 2, 3, 4095, 4096, 4097, 8195] {
            let source: Vec<usize> = (0..len)
                .map(|i| if i % 7 == 0 { u32::MAX as usize } else { i })
                .collect();
            let pointer = source.as_ptr();
            let expected: Vec<u32> = source.iter().map(|&id| id as u32).collect();
            let compact = CompactTranscriptIds::new(source);
            assert_eq!(&*compact, expected);
            assert_eq!(compact.words.as_ptr(), pointer);
        }
    }
}
