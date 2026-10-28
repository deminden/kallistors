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
// Three-base groups keep the initializer table within 512 bytes.
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

// Paired initialization shares each validated forward-word table lookup.
// The second lane contains the exact REP state of that group's reverse
// complement; opposite rotations place it in the reversed complete window.
const fn paired_groups<const N: usize>(bases: usize) -> [(u64, u64); N] {
    let mut out = [(0, 0); N];
    let mut word = 0;
    while word < N {
        let mut reverse = 0;
        let mut remaining = word;
        let mut i = 0;
        while i < bases {
            reverse = (reverse << 2) | ((remaining & 3) ^ 3);
            remaining >>= 2;
            i += 1;
        }
        out[word] = (
            REP_HASH_GROUPS[bases][word],
            REP_HASH_GROUPS[bases][reverse],
        );
        word += 1;
    }
    out
}

static REP_PAIRED_THREE: [(u64, u64); 64] = paired_groups(3);
static REP_PAIRED_TWO: [(u64, u64); 16] = paired_groups(2);

#[inline]
pub(crate) fn paired_state_23(word: u64) -> (u64, u64) {
    let tail = REP_PAIRED_TWO[((word >> 42) & 15) as usize];
    let mut forward = tail.0.rotate_left(21);
    let mut reverse = tail.1;
    macro_rules! add_group {
        ($offset:literal) => {{
            let pair = REP_PAIRED_THREE[((word >> (2 * $offset)) & 63) as usize];
            forward ^= pair.0.rotate_left($offset);
            reverse ^= pair.1.rotate_left(20 - $offset);
        }};
    }
    add_group!(0);
    add_group!(3);
    add_group!(6);
    add_group!(9);
    add_group!(12);
    add_group!(15);
    add_group!(18);
    (forward, reverse)
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
