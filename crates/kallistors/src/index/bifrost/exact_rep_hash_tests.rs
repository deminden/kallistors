#[test]
fn fixed_minimizer_keys_equal_generic_bytes_for_all_bits_and_seeds() {
    let check = |word: u64| {
        let bytes = word.to_le_bytes();
        for seed in [0, 1, u64::MAX] {
            assert_eq!(
                super::minimizer_hash(&bytes, seed),
                super::wyhash(&bytes, 0) ^ seed
            );
        }
    };
    for bit in 0..64 {
        check(1u64 << bit);
        check(!(1u64 << bit));
    }
    for word in [0, u64::MAX, u32::MAX as u64, 1u64 << 32] {
        check(word);
    }
    let mut state = 59u64;
    for _ in 0..4096 {
        state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
        check(state);
    }
}

#[test]
fn fixed_words_equal_generic_byte_hash_in_both_orders() {
    let mut state = 17u64;
    for first in [0, 1, u32::MAX as u64, 1u64 << 32, u64::MAX] {
        for second in [0, 1, u64::MAX] {
            check(first, second);
        }
    }
    for _ in 0..4096 {
        state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
        let first = state;
        state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
        check(first, state);
        check(state, first);
    }
}
fn check(first: u64, second: u64) {
    let mut bytes = [0u8; 16];
    bytes[..8].copy_from_slice(&first.to_le_bytes());
    bytes[8..].copy_from_slice(&second.to_le_bytes());
    assert_eq!(
        super::wyhash16_words(first, second),
        super::wyhash(&bytes, 0)
    );
}
