//! Release fully unused pages from a deduplicated scratch allocation.
#[cfg(any(target_os = "macos", target_os = "linux"))]
pub(crate) fn discard_unused<T>(words: &mut Vec<T>, written: usize) {
    let written = written.min(words.capacity());
    if written <= words.len() {
        return;
    }
    // SAFETY: sysconf has no pointer arguments or memory side effects.
    let Ok(page) = usize::try_from(unsafe { libc::sysconf(libc::_SC_PAGESIZE) }) else {
        return;
    };
    if page == 0 {
        return;
    }
    let base = words.as_mut_ptr() as usize;
    let Some(live_end) = words
        .len()
        .checked_mul(std::mem::size_of::<T>())
        .and_then(|len| base.checked_add(len))
    else {
        return;
    };
    let Some(written_end) = written
        .checked_mul(std::mem::size_of::<T>())
        .and_then(|len| base.checked_add(len))
    else {
        return;
    };
    let Some(first) = live_end.checked_add(page - 1).map(|end| end / page * page) else {
        return;
    };
    let last = written_end / page * page;
    if first >= last {
        return;
    }
    // SAFETY: all collector/sort workers have joined. This page-aligned range
    // lies wholly inside this allocation's spare capacity, beyond every live
    // element and before its original initialized end. No allocator headers,
    // live bytes or shared references occupy it. Mapping it read/write again
    // keeps the original address reservation accessible for reuse/deallocation.
    // Failure leaves the optimization unapplied; no live contents are modified.
    unsafe {
        #[cfg(target_os = "macos")]
        let _ = libc::mmap(
            first as *mut libc::c_void,
            last - first,
            libc::PROT_READ | libc::PROT_WRITE,
            libc::MAP_PRIVATE | libc::MAP_ANON | libc::MAP_FIXED,
            -1,
            0,
        );
        #[cfg(target_os = "linux")]
        let _ = libc::madvise(
            first as *mut libc::c_void,
            last - first,
            libc::MADV_DONTNEED,
        );
    }
}

#[cfg(not(any(target_os = "macos", target_os = "linux")))]
pub(crate) fn discard_unused<T>(_words: &mut Vec<T>, _written: usize) {}

// Clear first so destructors finish before releasing the allocation's wholly
// unused pages. The Vec retains ownership through normal allocator deallocation.
pub(crate) fn retire<T>(mut buffer: Vec<T>) {
    let written = buffer.len();
    buffer.clear();
    discard_unused(&mut buffer, written);
}

#[cfg(all(test, any(target_os = "macos", target_os = "linux")))]
mod tests {
    use super::*;
    #[test]
    fn discarded_scratch_preserves_partial_live_pages_and_capacity_reuse() {
        let mut words = (0..131072u64).collect::<Vec<_>>();
        let written = words.len();
        let retained = 4099;
        words.truncate(retained);
        discard_unused(&mut words, written);
        assert!(
            words
                .iter()
                .enumerate()
                .all(|(i, &value)| value == i as u64)
        );
        words.resize(written, 42);
        assert!(
            words[..retained]
                .iter()
                .enumerate()
                .all(|(i, &value)| value == i as u64)
        );
        assert!(words[retained..].iter().all(|&value| value == 42));
        let mut bytes = vec![7u8; 131072];
        let written = bytes.len();
        bytes.truncate(16387);
        discard_unused(&mut bytes, written);
        assert!(bytes.iter().all(|&value| value == 7));
        bytes.resize(written, 11);
        assert!(bytes[16387..].iter().all(|&value| value == 11));
        retire(bytes);
        let mut indices = vec![7u32; 32768];
        let written = indices.len();
        indices.truncate(4099);
        discard_unused(&mut indices, written);
        assert!(indices.iter().all(|&value| value == 7));
        indices.resize(written, 11);
        assert!(indices[4099..].iter().all(|&value| value == 11));
        retire(indices);
        let mut empty = Vec::<()>::with_capacity(4096);
        empty.resize(4096, ());
        empty.clear();
        discard_unused(&mut empty, 4096);
    }
}
