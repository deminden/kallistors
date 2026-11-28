//! Maps the same paired FASTQ through kallistors' own reader and through a
//! readkit batch source, and checks the two results are identical.
//!
//! cargo run --release --features readkit --example readkit_parity -- \
//!     INDEX R1.fastq R2.fastq [THREADS]

use std::fs::File;
use std::io::BufReader;
use std::path::Path;
use std::time::Instant;

use kallistors::io::open_fastq_reader;
use kallistors::pseudoalign::{
    EcCounts, PseudoalignOptions, Strand, build_bifrost_index_quant_threaded,
    pseudoalign_paired_batches_threaded, pseudoalign_paired_bifrost_with_options_threaded,
};
use readkit::fastq::{FastqInput, FastqOptions, FastqSource};

/// Nonzero classes (sorted), reads processed, reads aligned, length histogram.
type Summary = (Vec<(Vec<u32>, u32)>, u64, u64, Option<Vec<u32>>);

fn summary(counts: &EcCounts) -> Summary {
    let mut classes: Vec<_> = counts
        .ec_list
        .iter()
        .cloned()
        .zip(counts.counts.iter().copied())
        .filter(|(_, count)| *count > 0)
        .collect();
    classes.sort_unstable();
    (
        classes,
        counts.reads_processed,
        counts.reads_aligned,
        counts.fragment_length_hist.clone(),
    )
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let [index, r1, r2] = [&args[0], &args[1], &args[2]];
    let threads: usize = args.get(3).map_or(Ok(8), |t| t.parse())?;
    let index = build_bifrost_index_quant_threaded(Path::new(index), false, threads)?;
    let options = PseudoalignOptions::default();

    let start = Instant::now();
    let mut reader1 = open_fastq_reader(Path::new(r1))?;
    let mut reader2 = open_fastq_reader(Path::new(r2))?;
    let files = pseudoalign_paired_bifrost_with_options_threaded(
        &index,
        &mut reader1,
        &mut reader2,
        Strand::Unstranded,
        options,
        threads,
    )?;
    println!("kallistors reader: {:.3}s", start.elapsed().as_secs_f64());

    let start = Instant::now();
    let open = |path: &str| -> std::io::Result<FastqInput> {
        Ok(FastqInput::named(
            BufReader::with_capacity(1 << 20, File::open(path)?),
            path,
        ))
    };
    let mut source = FastqSource::paired(open(r1)?, open(r2)?, FastqOptions::default());
    let batches = pseudoalign_paired_batches_threaded(
        &index,
        &mut source,
        Strand::Unstranded,
        options,
        threads,
    )?;
    println!("readkit source:    {:.3}s", start.elapsed().as_secs_f64());

    let (a, b) = (summary(&files), summary(&batches));
    println!(
        "classes {} / {}, processed {} / {}, aligned {} / {}",
        a.0.len(),
        b.0.len(),
        a.1,
        b.1,
        a.2,
        b.2
    );
    println!("EC counts identical: {}", a.0 == b.0);
    println!("fragment-length histogram identical: {}", a.3 == b.3);
    if a != b {
        return Err("results differ".into());
    }
    Ok(())
}
