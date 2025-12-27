# Kallisto compatibility TODOs

This checklist is ordered. Do not proceed to the next item until the current item is validated against kallisto output.

## Phase A: Index parsing parity (metadata)
- [x] Parse and report index version (should be 13)
- [x] Parse and report k-mer length (match `kallisto inspect`)
- [x] Parse and report minimizer length (match `kallisto inspect`)
- [x] Parse and report number of unitigs (match `kallisto inspect`)
- [x] Parse and report number of k-mers (match `kallisto inspect`)
- [x] Parse and report number of targets and total length (match `kallisto inspect` + transcript FASTA)

Acceptance test:
- `kallisto inspect data/synthetic.idx` matches `kallistors-cli index-info` for all reported fields.

## Phase B: EC data parity
- [x] Parse EC mapping and emit list in kallisto-compatible order (debug-only)
- [x] Support loading kallisto EC text files (`output.ec.txt`) and validate ordering
- [x] Export EC list to text and compare with kallisto EC list (if available)
- [x] Add debug-only EC compare utility (remove later)

Notes:
- kallisto does not expose an EC list export for a plain index; comparison requires an external tool or custom build.

Acceptance test:
- EC list checksum matches a known reference for a synthetic index.

## Phase C: Pseudoalignment parity (single-end)
- [x] Implement naive k-mer lookup against kallisto index structures (debug path)
- [x] Compute EC per read and count (debug path)
- [ ] Validate EC counts for toy reads against kallisto pseudoalignment output

Acceptance test:
- EC counts identical for a fixed synthetic reads set.

Notes:
- kallisto does not expose EC counts directly; validation will require a custom tool or patched build.

## Phase D: Quant parity (basic EM)
- [x] Implement EM on EC counts
- [x] Validate against kallisto quant on toy data (within tolerance)

Acceptance test:
- TPMs/est_counts within tolerance for synthetic data.

Notes:
- `kallisto quant --bias` segfaults on the current simple synthetic dataset; bias parity remains blocked upstream.
