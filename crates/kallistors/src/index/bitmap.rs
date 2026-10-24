//! Bitmap decoding for native node payloads and portable-compatible containers.
//!
//! The native array path keeps exact reservation and set normalisation. The
//! compatibility path tries portable decoding first, then native decoding.
//! These format policies intentionally remain distinct.

use croaring_sys::{
    roaring_bitmap_deserialize, roaring_bitmap_free, roaring_bitmap_get_cardinality,
    roaring_bitmap_portable_deserialize_safe, roaring_bitmap_t, roaring_bitmap_to_uint32_array,
};

#[inline]
unsafe fn copy_and_free(bitmap: *mut roaring_bitmap_t, destination: *mut u32) {
    unsafe {
        roaring_bitmap_to_uint32_array(bitmap, destination);
        roaring_bitmap_free(bitmap);
    }
}

pub(super) unsafe fn deserialize_native_to_vec(buf: &[u8]) -> Option<Vec<u32>> {
    // Immutable full-node metadata becomes boxed after decoding. Reserve the
    // exact small-array length to avoid a subsequent shrink reallocation.
    // Validate its complete payload before allocating, as the decoder does.
    let capacity = if buf.first() == Some(&1) {
        let count = u32::from_le_bytes(buf.get(1..5)?.try_into().ok()?) as usize;
        let end = count.checked_mul(4)?.checked_add(5)?;
        buf.get(5..end)?;
        count
    } else {
        0
    };
    let mut out = Vec::with_capacity(capacity);
    unsafe { append_native_to_vec(buf, &mut out) }?;
    Some(out)
}

pub(crate) unsafe fn append_native_to_vec(buf: &[u8], out: &mut Vec<u32>) -> Option<()> {
    if buf.first() == Some(&1) {
        let count = u32::from_le_bytes(buf.get(1..5)?.try_into().ok()?) as usize;
        let end = count.checked_mul(4)?.checked_add(5)?;
        let payload = buf.get(5..end)?;
        let start = out.len();
        out.extend(
            payload
                .as_chunks::<4>()
                .0
                .iter()
                .map(|word| u32::from_le_bytes(*word)),
        );
        if out[start..].windows(2).any(|pair| pair[0] >= pair[1]) {
            out[start..].sort_unstable();
            let mut write = start;
            for read in start..out.len() {
                let value = out[read];
                if write == start || out[write - 1] != value {
                    out[write] = value;
                    write += 1;
                }
            }
            out.truncate(write);
        }
        return Some(());
    }
    if buf.is_empty() {
        return None;
    }
    let ptr = unsafe { roaring_bitmap_deserialize(buf.as_ptr().cast()) };
    if ptr.is_null() {
        return None;
    }
    let card = unsafe { roaring_bitmap_get_cardinality(ptr) } as usize;
    let start = out.len();
    out.resize(start + card, 0);
    unsafe {
        copy_and_free(ptr, out.as_mut_ptr().add(start));
    }
    Some(())
}

pub(crate) unsafe fn deserialize_compatible_to_vec(buf: &[u8]) -> Option<Vec<u32>> {
    let ptr = unsafe { roaring_bitmap_portable_deserialize_safe(buf.as_ptr().cast(), buf.len()) };
    let ptr = if ptr.is_null() {
        unsafe { roaring_bitmap_deserialize(buf.as_ptr().cast()) }
    } else {
        ptr
    };
    if ptr.is_null() {
        return None;
    }
    let card = unsafe { roaring_bitmap_get_cardinality(ptr) };
    let mut out = vec![0u32; card as usize];
    unsafe {
        copy_and_free(ptr, out.as_mut_ptr());
    }
    Some(out)
}
