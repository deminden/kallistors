//! Shared exact grouped REP initialization for index construction and mapping.
pub(crate) const REP_HASH_VALS: [u64; 4] = [
    2053695854357871005,
    5073395517033431291,
    10060236952204337488,
    7783083932390163561,
];
// Hash groups use the same A/C/T/G REP codes as the byte oracle, but the input
// is the validated A/C/G/T two-bit word. XOR/rotation composition is exact.
const REP_HASH_GROUPS: [[u64; 256]; 5] = {
    let mut groups = [[0u64; 256]; 5];
    let mut len = 1;
    while len <= 4 {
        let mut word = 0;
        while word < 256 {
            let mut hash = 0u64;
            let mut i = 0;
            while i < len {
                let base = (word >> (2 * (len - i - 1))) & 3;
                let rep_code = base ^ (base >> 1);
                hash = hash.rotate_left(1) ^ REP_HASH_VALS[rep_code];
                i += 1;
            }
            groups[len][word] = hash;
            word += 1;
        }
        len += 1;
    }
    groups
};
// Three-base groups keep the initializer table within512 bytes.
static REP_GROUP_THREE: [u64; 64] = {
    let mut groups = [0; 64];
    let mut i = 0;
    while i < 64 {
        groups[i] = REP_HASH_GROUPS[3][i];
        i += 1;
    }
    groups
};
#[inline]
fn rep_grouped_23(word: u64) -> u64 {
    REP_GROUP_THREE[(word & 63) as usize]
        ^ REP_GROUP_THREE[((word >> 6) & 63) as usize].rotate_left(3)
        ^ REP_GROUP_THREE[((word >> 12) & 63) as usize].rotate_left(6)
        ^ REP_GROUP_THREE[((word >> 18) & 63) as usize].rotate_left(9)
        ^ REP_GROUP_THREE[((word >> 24) & 63) as usize].rotate_left(12)
        ^ REP_GROUP_THREE[((word >> 30) & 63) as usize].rotate_left(15)
        ^ REP_GROUP_THREE[((word >> 36) & 63) as usize].rotate_left(18)
        ^ REP_HASH_GROUPS[2][((word >> 42) & 15) as usize].rotate_left(21)
}

/// Exact REP state from a validated right-aligned A/C/G/T two-bit word.
#[inline]
pub(crate) fn grouped_state(mut word: u64, g: usize) -> u64 {
    debug_assert!(g <= 32);
    if g == 23 {
        return rep_grouped_23(word);
    }
    let mut state = 0u64;
    for chunk in 0..g / 4 {
        state ^= REP_HASH_GROUPS[4][(word & 255) as usize].rotate_left((chunk * 4) as u32);
        word >>= 8;
    }
    let tail = g % 4;
    if tail != 0 {
        state ^= REP_HASH_GROUPS[tail][word as usize].rotate_left((g - tail) as u32);
    }
    state
}
