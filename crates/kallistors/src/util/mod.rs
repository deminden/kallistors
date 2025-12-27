//! Utility helpers for kallistors.

/// Returns true if the sequence contains only A/C/G/T (case-insensitive).
pub fn is_valid_dna(seq: &[u8]) -> bool {
    seq.iter()
        .all(|b| matches!(b, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't'))
}

/// Returns true if the sequence contains an ambiguous base (N/n).
pub fn contains_n(seq: &[u8]) -> bool {
    seq.iter().any(|b| *b == b'N' || *b == b'n')
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dna_validation_accepts_acgt() {
        assert!(is_valid_dna(b"ACGTacgt"));
        assert!(!is_valid_dna(b"ACGTN"));
    }

    #[test]
    fn contains_n_detects_ambiguous() {
        assert!(contains_n(b"ACNT"));
        assert!(!contains_n(b"ACGT"));
    }
}
