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

pub fn aa_to_comma_free(seq: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity(seq.len() * 3);
    for &aa in seq {
        out.extend_from_slice(aa_comma_free_codon(aa));
    }
    out
}

pub fn nucleotide_to_comma_free(seq: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity(seq.len());
    for codon in seq.chunks_exact(3) {
        out.extend_from_slice(nucleotide_comma_free_codon(codon));
    }
    out
}

pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|base| match base.to_ascii_uppercase() {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' | b'U' => b'A',
            _ => b'N',
        })
        .collect()
}

fn aa_comma_free_codon(aa: u8) -> &'static [u8; 3] {
    match aa.to_ascii_uppercase() {
        b'F' => b"ACC",
        b'L' | b'J' => b"ACA",
        b'I' => b"ATA",
        b'M' => b"ATC",
        b'V' => b"ATT",
        b'S' => b"CTA",
        b'P' => b"CTC",
        b'T' => b"CTT",
        b'A' => b"AGA",
        b'Y' => b"AGC",
        b'H' => b"AGT",
        b'Q' => b"AGG",
        b'N' => b"CGA",
        b'K' => b"CGC",
        b'D' | b'B' => b"CGT",
        b'E' | b'Z' => b"CGG",
        b'C' => b"TGA",
        b'W' => b"TGC",
        b'R' => b"TGT",
        b'G' => b"TGG",
        _ => b"NNN",
    }
}

fn nucleotide_comma_free_codon(codon: &[u8]) -> &'static [u8; 3] {
    let b0 = codon[0].to_ascii_uppercase();
    let b1 = codon[1].to_ascii_uppercase();
    let b2 = codon[2].to_ascii_uppercase();
    match (b0, b1, b2) {
        (b'T', b'T', b'T' | b'C') => b"ACC",
        (b'T', b'T', b'A' | b'G') | (b'C', b'T', _) => b"ACA",
        (b'A', b'T', b'T' | b'C' | b'A') => b"ATA",
        (b'A', b'T', b'G') => b"ATC",
        (b'G', b'T', _) => b"ATT",
        (b'T', b'C', _) | (b'A', b'G', b'T' | b'C') => b"CTA",
        (b'C', b'C', _) => b"CTC",
        (b'A', b'C', _) => b"CTT",
        (b'G', b'C', _) => b"AGA",
        (b'T', b'A', b'T' | b'C') => b"AGC",
        (b'C', b'A', b'T' | b'C') => b"AGT",
        (b'C', b'A', b'A' | b'G') => b"AGG",
        (b'A', b'A', b'T' | b'C') => b"CGA",
        (b'A', b'A', b'A' | b'G') => b"CGC",
        (b'G', b'A', b'T' | b'C') => b"CGT",
        (b'G', b'A', b'A' | b'G') => b"CGG",
        (b'T', b'G', b'T' | b'C') => b"TGA",
        (b'T', b'G', b'G') => b"TGC",
        (b'C', b'G', _) | (b'A', b'G', b'A' | b'G') => b"TGT",
        (b'G', b'G', _) => b"TGG",
        _ => b"NNN",
    }
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

    #[test]
    fn aa_to_comma_free_matches_kallisto_map() {
        assert_eq!(aa_to_comma_free(b"FLIMV"), b"ACCACAATAATCATT");
    }

    #[test]
    fn nucleotide_to_comma_free_translates_codons() {
        assert_eq!(nucleotide_to_comma_free(b"TTTCTTATGGGT"), b"ACCACAATCTGG");
    }
}
