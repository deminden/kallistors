use anyhow::{Result, anyhow};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

#[allow(clippy::too_many_arguments)]
pub fn run(
    index: &Path,
    ec: &Path,
    out: &Path,
    bias: bool,
    bias_counts: Option<&Path>,
    transcripts: Option<&Path>,
    fragment_length: Option<f64>,
    fragment_length_sd: Option<f64>,
    fr_stranded: bool,
    rf_stranded: bool,
) -> Result<()> {
    if fr_stranded && rf_stranded {
        return Err(anyhow!(
            "--fr-stranded and --rf-stranded are mutually exclusive"
        ));
    }
    if bias && transcripts.is_none() {
        return Err(anyhow!("--bias requires --transcripts"));
    }

    let index_meta = kallistors::index::Index::load(index)?;
    let ec_input = kallistors::quant::load_ec_counts(ec)?;

    let bias_counts = if bias {
        let path = bias_counts
            .map(|p| p.to_path_buf())
            .unwrap_or_else(|| ec.with_extension("bias.txt"));
        Some(kallistors::bias::BiasCounts::load_text(&path)?)
    } else {
        None
    };

    let transcript_seqs = if bias {
        Some(kallistors::quant::load_transcript_sequences(
            transcripts.unwrap(),
            &index_meta
                .transcripts
                .iter()
                .map(|t| t.name.clone())
                .collect::<Vec<_>>(),
        )?)
    } else {
        None
    };

    let strand_specific = if fr_stranded {
        Some(kallistors::pseudoalign::StrandSpecific::FR)
    } else if rf_stranded {
        Some(kallistors::pseudoalign::StrandSpecific::RF)
    } else {
        None
    };

    let result = kallistors::quant::em_quantify(
        &ec_input,
        &index_meta
            .transcripts
            .iter()
            .map(|t| t.length)
            .collect::<Vec<_>>(),
        transcript_seqs.as_deref(),
        bias_counts.as_ref(),
        kallistors::quant::QuantOptions {
            mean_fragment_length: fragment_length.unwrap_or(200.0),
            fragment_length_sd: fragment_length_sd.unwrap_or(20.0),
            bias,
            strand_specific,
            ..kallistors::quant::QuantOptions::default()
        },
    )?;

    let out_file = File::create(out)?;
    let mut writer = BufWriter::new(out_file);
    writeln!(writer, "transcript_id\tname\tlength\test_counts\ttpm")?;
    for (i, tr) in index_meta.transcripts.iter().enumerate() {
        let est = result.est_counts.get(i).copied().unwrap_or(0.0);
        let tpm = result.tpm.get(i).copied().unwrap_or(0.0);
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}",
            i, tr.name, tr.length, est, tpm
        )?;
    }
    Ok(())
}
