# Compatibility

kallistors supports kallisto v13 indexes, paired/single-end quantification and the BUS modes
described below. See the [README](../README.md#usage) for common commands and the
[command-line guide](usage.md) for additional examples. Unlisted kallisto options require
separate compatibility checks.

## Quantification

- Complete EC parity means identical counts for each transcript set, independent of EC numbering.
  Equal `run_info.json` totals alone do not establish EC parity.
- `--fr-stranded` and `--rf-stranded` retain strand annotations. A mapped mate rejected by
  strand rejects the pair; a genuinely unmapped mate can still allow its matching partner.
- Paired effective lengths use the empirical fragment histogram. The default
  [sampling policy](fragment_length_sampling.md) follows original kallisto's 10,000-observation
  goal and 8 MiB sequence buffers. Concurrent scheduling can change the sampled histogram.
- Single-end quantification requires a fragment mean and standard deviation (`-l` and `-s`).
- Different fragment models can change effective lengths, EM convergence and ambiguous isoform
  estimates even when EC counts match. Matching upstream behavior establishes compatibility;
  biological accuracy requires independent evidence.
- Sequence-specific bias correction is supported with `--bias --transcripts FASTA`.
  General real-data abundance parity with upstream bias correction is not established.
- Bootstrap HDF5 layout is compatible. Bootstrap RNG and sample values can differ from upstream.
- BAM output is supported by `bus`. The `quant --pseudobam` and `quant --genomebam` flags
  are recognized but return an unsupported-feature error.
- General long-read quantification and fusion detection are not implemented.

## Indexes

The Rust builder supports nucleotide transcript FASTA and amino-acid BUS indexes (`index --aa`).
It is checked against standard upstream builds using `MAX_KMER_SIZE=32`.
Distinguish mode, genomic D-list construction and nonstandard upstream compile-time k-mer widths
are outside its supported surface. See [index format](index_format_notes.md) for layout and
platform assumptions. Format-compatible indexes need not be byte-identical.

## BUS

Supported features include fixed presets and custom technology triples, batch/interleaved FASTQ,
BAM input tags, SmartSeq3 tags, long-read BUS filtering, amino-acid BUS mapping and transcriptome
or genome BAM output.

Nucleotide BUS defaults to unstranded processing. `--fr-stranded` retains cDNA reads
in transcript orientation; `--rf-stranded` retains reverse-complement cDNA reads.
For paired technologies these require opposite orientations for the second mate.
Explicit strand flags load transcript orientation annotations for FASTQ and BAM
input; annotations are read from existing v13 indexes without rebuilding them.
Choose the flag from the library protocol; the technology preset alone does not
infer strandedness. Ordinary stranded and unstranded counting use packed graph
storage while preserving the annotations needed by their filters.

Ordinary FASTQ BUS counting uses ordered parallel mapping for single-cDNA,
compound-barcode, synthetic-barcode, tagged and paired layouts. Ordinary BAM
input is decoded in bounded batches for the same counting workers. Read limits,
ambiguity flags, raw/corrected tag precedence and output order are preserved.
Paired BUS records every qualifying fragment length, independently of the
sampling quota used for bulk fragment estimation. Modes requesting long-read
filtering, amino-acid mapping, unmapped-read diagnostics or BAM output retain
the metadata required by those features.

Amino-acid mode supports single-cDNA technologies and rejects paired BUS technologies.
Unpaired multi-sequence custom technologies are rejected for BAM output and strand-specific
processing. Genome BAM output requires matching GTF annotations and chromosome lengths.

## Output

Quantification writes `abundance.tsv` and `run_info.json`; HDF5 is enabled by default.
`--plaintext` writes bootstrap TSVs instead of HDF5. HDF5 uses DEFLATE level 1 and chunks
bounded to 16,384 elements. Dataset schemas and values can be checked independently with
`scripts/check_hdf5_output.py` using native HDF5 through h5py.

Bootstrap scheduling shares the requested worker budget across samples and large within-sample
components. Sampling retains the original rejection checks, draws and final RNG state.
For the same supplied EC matrix, seed and effective lengths, serial and parallel bootstrap
quantification retain the same samples and native numeric results. FASTQ runs can assign
EC numbers in a different order across thread counts, so a fixed seed alone does not imply
identical bootstrap values between such runs or between different quantifiers.
