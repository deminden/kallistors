//! Quantification routines (EM over ECs with optional sequence bias).

use std::borrow::Cow;
use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use crate::bias::{BiasCounts, hexamer_to_int, update_hexamer};
use crate::ec::EcList;
use crate::index::Index;
use crate::pseudoalign::StrandSpecific;
use crate::{Error, Result};
use hdf5_pure::FileBuilder;
use hdf5_pure::datatype::{CharacterSet, Datatype, StringPadding};
use hdf5_pure::type_builders::DatasetBuilder;

const MIN_ALPHA: f64 = 1e-8;
const ALPHA_LIMIT: f64 = 1e-7;
const ALPHA_CHANGE_LIMIT: f64 = 1e-2;
const ALPHA_CHANGE: f64 = 1e-2;
const TOLERANCE: f64 = f64::from_bits(1);
const KALLISTO_BIAS_LEN: usize = 4096;

struct EmMultiEc {
    start: usize,
    len: usize,
    count: f64,
}

struct EmWork {
    singleton_alpha: Vec<f64>,
    multi_ecs: Vec<EmMultiEc>,
    multi_transcripts: Vec<usize>,
    multi_weights: Vec<f64>,
}

#[derive(Debug, Clone, Copy)]
pub struct QuantOptions {
    pub mean_fragment_length: f64,
    pub fragment_length_sd: f64,
    pub max_iter: usize,
    pub min_rounds: usize,
    pub bias: bool,
    pub strand_specific: Option<StrandSpecific>,
}

impl Default for QuantOptions {
    fn default() -> Self {
        Self {
            mean_fragment_length: 200.0,
            fragment_length_sd: 20.0,
            max_iter: 10_000,
            min_rounds: 50,
            bias: false,
            strand_specific: None,
        }
    }
}

pub struct QuantResult {
    pub est_counts: Vec<f64>,
    pub tpm: Vec<f64>,
    pub eff_lengths: Vec<f64>,
    pub post_bias: Option<Vec<f64>>,
}

pub struct H5QuantOutput<'a> {
    pub target_names: &'a [&'a str],
    pub target_lengths: &'a [u32],
    pub result: &'a QuantResult,
    pub bootstrap_est_counts: &'a [&'a [f64]],
    pub reads_processed: u64,
    pub fragment_length_hist: Option<&'a [u32]>,
    pub bias_observed: Option<&'a [u32]>,
    pub index_version: u64,
    pub kallisto_version: &'a str,
    pub call: &'a str,
    pub start_time: &'a str,
}

/// Run metadata matching kallisto's `run_info.json` fields.
pub struct RunInfo {
    pub n_targets: usize,
    pub n_bootstraps: usize,
    pub n_processed: u64,
    pub n_pseudoaligned: u64,
    pub n_unique: u64,
    pub p_pseudoaligned: f64,
    pub p_unique: f64,
    pub kallisto_version: String,
    pub index_version: u64,
    pub kmer_length: u32,
    pub start_time: String,
    pub call: String,
}

impl RunInfo {
    pub fn write_json(&self, path: &Path) -> Result<()> {
        let file = File::create(path)?;
        let mut writer = BufWriter::new(file);
        let entries = vec![
            ("n_targets", self.n_targets.to_json()),
            ("n_bootstraps", self.n_bootstraps.to_json()),
            ("n_processed", self.n_processed.to_json()),
            ("n_pseudoaligned", self.n_pseudoaligned.to_json()),
            ("n_unique", self.n_unique.to_json()),
            ("p_pseudoaligned", format!("{:.1}", self.p_pseudoaligned)),
            ("p_unique", format!("{:.1}", self.p_unique)),
            ("kallisto_version", self.kallisto_version.to_json()),
            ("index_version", self.index_version.to_json()),
            ("k-mer length", self.kmer_length.to_json()),
            ("start_time", self.start_time.to_json()),
            ("call", self.call.to_json()),
        ];
        writeln!(writer, "{{")?;
        for (idx, (key, value)) in entries.iter().enumerate() {
            let sep = if idx + 1 == entries.len() { "" } else { "," };
            writeln!(writer, "\"{}\": {}{}", escape_json(key), value, sep)?;
        }
        writeln!(writer, "}}")?;
        Ok(())
    }
}

/// Write kallisto-compatible `abundance.tsv`.
pub fn write_abundance_tsv(path: &Path, index: &Index, result: &QuantResult) -> Result<()> {
    let names = index
        .transcripts
        .iter()
        .map(|tr| tr.name.as_str())
        .collect::<Vec<_>>();
    let lengths = index
        .transcripts
        .iter()
        .map(|tr| tr.length)
        .collect::<Vec<_>>();
    write_abundance_tsv_from_parts(path, &names, &lengths, result)
}

/// Write kallisto-compatible `abundance.tsv` from already loaded transcript metadata.
pub fn write_abundance_tsv_from_parts(
    path: &Path,
    target_names: &[&str],
    target_lengths: &[u32],
    result: &QuantResult,
) -> Result<()> {
    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);
    writeln!(writer, "target_id\tlength\teff_length\test_counts\ttpm")?;
    for (i, name) in target_names.iter().enumerate() {
        let length = target_lengths.get(i).copied().unwrap_or(0);
        let eff_len = result.eff_lengths.get(i).copied().unwrap_or(0.0);
        let est = result.est_counts.get(i).copied().unwrap_or(0.0);
        let tpm = result.tpm.get(i).copied().unwrap_or(0.0);
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}",
            name, length, eff_len, est, tpm
        )?;
    }
    Ok(())
}

pub fn write_abundance_h5(path: &Path, output: H5QuantOutput<'_>) -> Result<()> {
    let mut file = FileBuilder::new();
    let mut aux = file.create_group("aux");

    write_h5_f64(file.create_dataset("est_counts"), &output.result.est_counts);
    write_h5_i32(
        aux.create_dataset("num_bootstrap"),
        &[
            i32::try_from(output.bootstrap_est_counts.len()).map_err(|_| {
                Error::InvalidFormat("too many bootstrap samples for HDF5 output".into())
            })?,
        ],
    );
    write_h5_i32(
        aux.create_dataset("num_processed"),
        &[i32::try_from(output.reads_processed).map_err(|_| {
            Error::InvalidFormat("n_processed exceeds kallisto-compatible HDF5 range".into())
        })?],
    );
    write_h5_i32(
        aux.create_dataset("fld"),
        &output
            .fragment_length_hist
            .map(u32_slice_to_i32)
            .transpose()?
            .unwrap_or_default(),
    );
    let bias_observed = match output.bias_observed {
        Some(values) => Cow::Owned(u32_slice_to_i32(values)?),
        None => Cow::Owned(vec![1; KALLISTO_BIAS_LEN]),
    };
    write_h5_i32(aux.create_dataset("bias_observed"), &bias_observed);

    let bias_normalized = match output.result.post_bias.as_deref() {
        Some(values) => Cow::Borrowed(values),
        None => Cow::Owned(vec![1.0; KALLISTO_BIAS_LEN]),
    };
    write_h5_f64(aux.create_dataset("bias_normalized"), &bias_normalized);
    write_h5_string_vec(
        aux.create_dataset("kallisto_version"),
        &[output.kallisto_version],
    );
    write_h5_i32(
        aux.create_dataset("index_version"),
        &[i32::try_from(output.index_version).map_err(|_| {
            Error::InvalidFormat("index version exceeds kallisto-compatible HDF5 range".into())
        })?],
    );
    write_h5_string_vec(aux.create_dataset("call"), &[output.call]);
    write_h5_string_vec(aux.create_dataset("start_time"), &[output.start_time]);
    write_h5_string_vec(aux.create_dataset("ids"), output.target_names);
    write_h5_f64(
        aux.create_dataset("eff_lengths"),
        &output.result.eff_lengths,
    );
    write_h5_i32(
        aux.create_dataset("lengths"),
        &u32_slice_to_i32(output.target_lengths)?,
    );
    file.add_group(aux.finish());

    if !output.bootstrap_est_counts.is_empty() {
        let mut bootstrap = file.create_group("bootstrap");
        for (idx, counts) in output.bootstrap_est_counts.iter().enumerate() {
            write_h5_f64(bootstrap.create_dataset(&format!("bs{idx}")), counts);
        }
        file.add_group(bootstrap.finish());
    }

    let bytes = file.finish().map_err(h5_error)?;
    fs::write(path, bytes)?;
    Ok(())
}

fn write_h5_f64(dataset: &mut DatasetBuilder, values: &[f64]) {
    dataset
        .with_f64_data(values)
        .with_shape(&[values.len() as u64]);
    set_h5_chunking(dataset, values.len());
}

fn write_h5_i32(dataset: &mut DatasetBuilder, values: &[i32]) {
    dataset
        .with_i32_data(values)
        .with_shape(&[values.len() as u64]);
    set_h5_chunking(dataset, values.len());
}

fn write_h5_string_vec(dataset: &mut DatasetBuilder, values: &[&str]) {
    let width = values
        .iter()
        .map(|value| value.len())
        .max()
        .unwrap_or(0)
        .saturating_add(1)
        .max(1);
    dataset
        .with_compound_data(
            fixed_ascii_type(width),
            fixed_string_bytes(values, width),
            values.len() as u64,
        )
        .with_shape(&[values.len() as u64]);
    set_h5_chunking(dataset, values.len());
}

fn set_h5_chunking(dataset: &mut DatasetBuilder, len: usize) {
    if len > 0 {
        dataset.with_chunks(&[len as u64]).with_deflate(6);
    }
}

fn fixed_ascii_type(width: usize) -> Datatype {
    Datatype::String {
        size: width as u32,
        padding: StringPadding::NullTerminate,
        charset: CharacterSet::Ascii,
    }
}

fn fixed_string_bytes(values: &[&str], width: usize) -> Vec<u8> {
    let mut out = vec![0u8; values.len() * width];
    for (row, value) in values.iter().enumerate() {
        let start = row * width;
        let bytes = value.as_bytes();
        out[start..start + bytes.len()].copy_from_slice(bytes);
    }
    out
}

fn u32_slice_to_i32(values: &[u32]) -> Result<Vec<i32>> {
    values
        .iter()
        .map(|&value| {
            i32::try_from(value).map_err(|_| {
                Error::InvalidFormat("value exceeds kallisto-compatible HDF5 range".into())
            })
        })
        .collect()
}

fn h5_error(err: hdf5_pure::Error) -> Error {
    Error::InvalidFormat(format!("HDF5 output failed: {err}"))
}

pub struct EcCountsInput {
    pub ec_list: EcList,
    pub counts: Vec<u32>,
}

pub fn load_ec_counts(path: &Path) -> Result<EcCountsInput> {
    let file = File::open(path).map_err(|_| Error::MissingFile(path.to_path_buf()))?;
    let reader = BufReader::new(file);
    let mut classes: Vec<Vec<u32>> = Vec::new();
    let mut counts: Vec<u32> = Vec::new();

    for (line_no, line) in reader.lines().enumerate() {
        let line = line?;
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }
        let mut parts = line.splitn(3, '\t');
        let ec_id: usize = parts
            .next()
            .ok_or_else(|| Error::InvalidFormat(format!("missing EC id at line {}", line_no + 1)))?
            .parse()
            .map_err(|_| Error::InvalidFormat(format!("invalid EC id at line {}", line_no + 1)))?;
        let count: u32 = parts
            .next()
            .ok_or_else(|| {
                Error::InvalidFormat(format!("missing EC count at line {}", line_no + 1))
            })?
            .parse()
            .map_err(|_| {
                Error::InvalidFormat(format!("invalid EC count at line {}", line_no + 1))
            })?;
        let list_str = parts
            .next()
            .ok_or_else(|| {
                Error::InvalidFormat(format!("missing EC list at line {}", line_no + 1))
            })?
            .trim();
        let mut ids = Vec::new();
        if !list_str.is_empty() {
            for token in list_str.split(',') {
                if token.is_empty() {
                    continue;
                }
                let val: u32 = token.parse().map_err(|_| {
                    Error::InvalidFormat(format!("invalid transcript id at line {}", line_no + 1))
                })?;
                ids.push(val);
            }
        }
        if ec_id != classes.len() {
            return Err(Error::InvalidFormat(format!(
                "EC id out of order at line {}",
                line_no + 1
            )));
        }
        classes.push(ids);
        counts.push(count);
    }

    Ok(EcCountsInput {
        ec_list: EcList { classes },
        counts,
    })
}

trait JsonValue {
    fn to_json(&self) -> String;
}

impl JsonValue for usize {
    fn to_json(&self) -> String {
        self.to_string()
    }
}

impl JsonValue for u64 {
    fn to_json(&self) -> String {
        self.to_string()
    }
}

impl JsonValue for u32 {
    fn to_json(&self) -> String {
        self.to_string()
    }
}

impl JsonValue for f64 {
    fn to_json(&self) -> String {
        if self.is_finite() {
            self.to_string()
        } else {
            "0".to_string()
        }
    }
}

impl JsonValue for bool {
    fn to_json(&self) -> String {
        self.to_string()
    }
}

impl JsonValue for String {
    fn to_json(&self) -> String {
        format!("\"{}\"", escape_json(self))
    }
}

fn escape_json(input: &str) -> String {
    let mut out = String::with_capacity(input.len());
    for ch in input.chars() {
        match ch {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c if c.is_control() => {
                out.push_str("\\u");
                out.push_str(&format!("{:04x}", c as u32));
            }
            c => out.push(c),
        }
    }
    out
}

pub fn load_transcript_sequences(fasta: &Path, names: &[String]) -> Result<Vec<Vec<u8>>> {
    let file = File::open(fasta).map_err(|_| Error::MissingFile(fasta.to_path_buf()))?;
    let reader = BufReader::new(file);
    let mut seqs: HashMap<String, Vec<u8>> = HashMap::new();
    let mut current_name: Option<String> = None;
    let mut current_seq: Vec<u8> = Vec::new();
    for line in reader.lines() {
        let line = line?;
        if let Some(rest) = line.strip_prefix('>') {
            if let Some(name) = current_name.take() {
                seqs.insert(name, std::mem::take(&mut current_seq));
            }
            current_name = Some(rest.trim().to_string());
        } else {
            current_seq.extend_from_slice(line.trim().as_bytes());
        }
    }
    if let Some(name) = current_name.take() {
        seqs.insert(name, current_seq);
    }

    let mut out = Vec::with_capacity(names.len());
    for name in names {
        let seq = seqs
            .get(name)
            .cloned()
            .ok_or_else(|| Error::InvalidFormat(format!("missing transcript {}", name)))?;
        out.push(seq);
    }
    Ok(out)
}

pub fn em_quantify(
    input: &EcCountsInput,
    transcript_lengths: &[u32],
    transcript_sequences: Option<&[Vec<u8>]>,
    bias_counts: Option<&BiasCounts>,
    options: QuantOptions,
) -> Result<QuantResult> {
    let num_trans = transcript_lengths.len();
    if num_trans == 0 {
        return Err(Error::InvalidFormat("no transcripts provided".into()));
    }
    if input.ec_list.classes.len() != input.counts.len() {
        return Err(Error::InvalidFormat(
            "EC list and counts length mismatch".into(),
        ));
    }

    let means = mean_frag_lens_by_transcript(
        transcript_lengths,
        options.mean_fragment_length,
        options.fragment_length_sd,
    );
    let mut eff_lens = calc_eff_lens(transcript_lengths, &means);
    let mut post_bias = None;

    let mut alpha = vec![1.0 / num_trans as f64; num_trans];
    let mut next_alpha = vec![0.0f64; num_trans];

    let mut em_work = prepare_em_work(&input.ec_list.classes, &input.counts, &eff_lens);

    let mut final_round = false;
    for iter in 0..options.max_iter {
        if options.bias
            && (iter == options.min_rounds || iter == options.min_rounds + 500)
            && let (Some(bias_counts), Some(transcript_sequences)) =
                (bias_counts, transcript_sequences)
        {
            let (updated_eff, bias_model) = update_eff_lens(
                &means,
                bias_counts,
                transcript_sequences,
                &alpha,
                &eff_lens,
                options.strand_specific,
            );
            eff_lens = updated_eff;
            post_bias = Some(bias_model);
            refresh_multi_weights(&mut em_work, &eff_lens);
        }

        next_alpha.copy_from_slice(&em_work.singleton_alpha);

        for ec in &em_work.multi_ecs {
            let range = ec.start..ec.start + ec.len;
            let mut denom = 0.0;
            for i in range.clone() {
                let tr = em_work.multi_transcripts[i];
                denom += alpha[tr] * em_work.multi_weights[i];
            }
            if denom < TOLERANCE {
                continue;
            }
            let scale = ec.count / denom;
            for i in range {
                let tr = em_work.multi_transcripts[i];
                next_alpha[tr] += em_work.multi_weights[i] * alpha[tr] * scale;
            }
        }

        let mut chcount = 0;
        for i in 0..num_trans {
            if next_alpha[i] > ALPHA_CHANGE_LIMIT {
                let rel = (next_alpha[i] - alpha[i]).abs() / next_alpha[i];
                if rel > ALPHA_CHANGE {
                    chcount += 1;
                }
            }
            alpha[i] = next_alpha[i];
        }

        let stop_em = chcount == 0 && iter > options.min_rounds;
        if final_round {
            break;
        }
        if stop_em {
            final_round = true;
            for val in &mut alpha {
                if *val < ALPHA_LIMIT / 10.0 {
                    *val = 0.0;
                }
            }
        }
    }

    let mut rates = Vec::with_capacity(num_trans);
    let mut rate_sum = 0.0;
    for i in 0..num_trans {
        let rate = if eff_lens[i] > 0.0 {
            alpha[i] / eff_lens[i]
        } else {
            0.0
        };
        rates.push(rate);
        rate_sum += rate;
    }

    let mut tpm = vec![0.0; num_trans];
    if rate_sum > 0.0 {
        for i in 0..num_trans {
            tpm[i] = rates[i] / rate_sum * 1_000_000.0;
        }
    }

    Ok(QuantResult {
        est_counts: alpha,
        tpm,
        eff_lengths: eff_lens,
        post_bias,
    })
}

pub fn bootstrap_quantify(
    input: &EcCountsInput,
    eff_lens: &[f64],
    bootstrap_samples: usize,
    seed: u64,
    options: QuantOptions,
) -> Result<Vec<QuantResult>> {
    if input.ec_list.classes.len() != input.counts.len() {
        return Err(Error::InvalidFormat(
            "EC list and counts length mismatch".into(),
        ));
    }

    let mut seed_rng = SplitMix64::new(seed);
    let mut out = Vec::with_capacity(bootstrap_samples);
    for _ in 0..bootstrap_samples {
        let mut sample_rng = SplitMix64::new(seed_rng.next_u64());
        let counts = sample_ec_counts(&input.counts, &mut sample_rng);
        out.push(em_quantify_fixed_eff_lens(
            &input.ec_list.classes,
            &counts,
            eff_lens,
            options.max_iter,
            options.min_rounds,
        ));
    }
    Ok(out)
}

fn em_quantify_fixed_eff_lens(
    ec_classes: &[Vec<u32>],
    counts: &[u32],
    eff_lens: &[f64],
    max_iter: usize,
    min_rounds: usize,
) -> QuantResult {
    let num_trans = eff_lens.len();
    let mut alpha = vec![1.0 / num_trans as f64; num_trans];
    let mut next_alpha = vec![0.0f64; num_trans];
    let em_work = prepare_em_work(ec_classes, counts, eff_lens);

    let mut final_round = false;
    for iter in 0..max_iter {
        next_alpha.copy_from_slice(&em_work.singleton_alpha);

        for ec in &em_work.multi_ecs {
            let range = ec.start..ec.start + ec.len;
            let mut denom = 0.0;
            for i in range.clone() {
                let tr = em_work.multi_transcripts[i];
                denom += alpha[tr] * em_work.multi_weights[i];
            }
            if denom < TOLERANCE {
                continue;
            }
            let scale = ec.count / denom;
            for i in range {
                let tr = em_work.multi_transcripts[i];
                next_alpha[tr] += em_work.multi_weights[i] * alpha[tr] * scale;
            }
        }

        let mut chcount = 0;
        for i in 0..num_trans {
            if next_alpha[i] > ALPHA_CHANGE_LIMIT {
                let rel = (next_alpha[i] - alpha[i]).abs() / next_alpha[i];
                if rel > ALPHA_CHANGE {
                    chcount += 1;
                }
            }
            alpha[i] = next_alpha[i];
        }

        let stop_em = chcount == 0 && iter > min_rounds;
        if final_round {
            break;
        }
        if stop_em {
            final_round = true;
            for val in &mut alpha {
                if *val < ALPHA_LIMIT / 10.0 {
                    *val = 0.0;
                }
            }
        }
    }

    let tpm = counts_to_tpm(&alpha, eff_lens);
    QuantResult {
        est_counts: alpha,
        tpm,
        eff_lengths: eff_lens.to_vec(),
        post_bias: None,
    }
}

fn prepare_em_work(ec_classes: &[Vec<u32>], counts: &[u32], eff_lens: &[f64]) -> EmWork {
    let mut singleton_alpha = vec![0.0f64; eff_lens.len()];
    let mut multi_ecs = Vec::new();
    let mut multi_transcripts = Vec::new();
    let mut multi_weights = Vec::new();
    for (ec_id, ec) in ec_classes.iter().enumerate() {
        let count = counts.get(ec_id).copied().unwrap_or(0) as f64;
        if ec.len() == 1 {
            let tr = ec[0] as usize;
            if tr < singleton_alpha.len() {
                singleton_alpha[tr] = count;
            }
        } else if ec.len() > 1 && count != 0.0 {
            let start = multi_transcripts.len();
            for &tr in ec {
                let tr = tr as usize;
                let len = eff_lens.get(tr).copied().unwrap_or(1.0);
                multi_transcripts.push(tr);
                multi_weights.push(count / len);
            }
            multi_ecs.push(EmMultiEc {
                start,
                len: ec.len(),
                count,
            });
        }
    }
    EmWork {
        singleton_alpha,
        multi_ecs,
        multi_transcripts,
        multi_weights,
    }
}

fn sample_ec_counts(counts: &[u32], rng: &mut SplitMix64) -> Vec<u32> {
    let total = counts.iter().map(|&v| u64::from(v)).sum::<u64>();
    let mut sampled = vec![0u32; counts.len()];
    if total == 0 {
        return sampled;
    }

    let mut cumulative = Vec::with_capacity(counts.len());
    let mut running = 0u64;
    for &count in counts {
        running = running.saturating_add(u64::from(count));
        cumulative.push(running);
    }

    for _ in 0..total {
        let draw = rng.gen_range(total);
        let idx = cumulative.partition_point(|&bound| bound <= draw);
        if let Some(slot) = sampled.get_mut(idx) {
            *slot = slot.saturating_add(1);
        }
    }
    sampled
}

fn counts_to_tpm(est_counts: &[f64], eff_lens: &[f64]) -> Vec<f64> {
    let mut tpm = vec![0.0; est_counts.len()];
    let mut total_mass = 0.0;
    for i in 0..est_counts.len() {
        let eff_len = eff_lens.get(i).copied().unwrap_or(0.0);
        if eff_len > 0.0 {
            tpm[i] = est_counts[i] / eff_len;
            total_mass += tpm[i];
        }
    }
    if total_mass > 0.0 {
        for val in &mut tpm {
            *val = *val / total_mass * 1_000_000.0;
        }
    }
    tpm
}

struct SplitMix64 {
    state: u64,
}

impl SplitMix64 {
    fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9e37_79b9_7f4a_7c15);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        z ^ (z >> 31)
    }

    fn gen_range(&mut self, upper: u64) -> u64 {
        let zone = u64::MAX - u64::MAX % upper;
        loop {
            let value = self.next_u64();
            if value < zone {
                return value % upper;
            }
        }
    }
}

fn refresh_multi_weights(work: &mut EmWork, eff_lens: &[f64]) {
    for ec in &work.multi_ecs {
        for i in ec.start..ec.start + ec.len {
            let tr = work.multi_transcripts[i];
            let len = eff_lens.get(tr).copied().unwrap_or(1.0);
            work.multi_weights[i] = ec.count / len;
        }
    }
}

fn calc_eff_lens(lengths: &[u32], means: &[f64]) -> Vec<f64> {
    let mut eff = Vec::with_capacity(lengths.len());
    for (i, &len) in lengths.iter().enumerate() {
        let mean = means.get(i).copied().unwrap_or(0.0);
        let mut cur = len as f64 - mean + 1.0;
        if cur < 1.0 {
            cur = len as f64;
        }
        eff.push(cur);
    }
    eff
}

fn mean_frag_lens_by_transcript(lengths: &[u32], mean: f64, sd: f64) -> Vec<f64> {
    const MAX_FRAG_LEN: usize = 1000;
    if sd <= 0.0 {
        return vec![mean; lengths.len()];
    }
    let mean_trunc = trunc_gaussian_fld(0, MAX_FRAG_LEN, mean, sd);
    let marginal = *mean_trunc.last().unwrap_or(&mean);
    lengths
        .iter()
        .map(|&len| {
            let len = len as usize;
            if len >= MAX_FRAG_LEN {
                marginal
            } else {
                mean_trunc[len]
            }
        })
        .collect()
}

fn trunc_gaussian_fld(start: usize, stop: usize, mean: f64, sd: f64) -> Vec<f64> {
    let n = stop.saturating_sub(start);
    let mut mean_fl = vec![0.0; n];
    let mut total_mass = 0.0;
    let mut total_density = 0.0;

    for (i, slot) in mean_fl.iter_mut().enumerate() {
        let mut x = (start + i) as f64;
        x = (x - mean) / sd;
        let cur_density = (-0.5 * x * x).exp() / sd;
        total_mass += cur_density * i as f64;
        total_density += cur_density;
        if total_mass > 0.0 {
            *slot = total_mass / total_density;
        }
    }
    mean_fl
}

fn update_eff_lens(
    means: &[f64],
    bias_counts: &BiasCounts,
    transcript_seqs: &[Vec<u8>],
    alpha: &[f64],
    eff_lens: &[f64],
    strand_specific: Option<StrandSpecific>,
) -> (Vec<f64>, Vec<f64>) {
    let num6 = 4096usize;
    let bias_data_norm: f64 = bias_counts.counts.iter().map(|&v| v as f64).sum();
    if bias_data_norm == 0.0 {
        return (eff_lens.to_vec(), vec![1.0; num6]);
    }

    let mut dbias5 = vec![0.0; num6];
    let mut bias_alpha_norm = 0.0;

    for i in 0..transcript_seqs.len() {
        if i >= means.len() || i >= alpha.len() || i >= eff_lens.len() {
            break;
        }
        if alpha[i] < MIN_ALPHA {
            continue;
        }
        let mean = means[i];
        let seq = &transcript_seqs[i];
        if seq.len() < 6 || (seq.len() as f64) < mean {
            continue;
        }
        let mut contrib = 0.5 * alpha[i] / eff_lens[i];
        if strand_specific.is_some() {
            contrib = alpha[i] / eff_lens[i];
        }
        let seqlen = seq.len() as isize;

        if (strand_specific.is_none() || strand_specific == Some(StrandSpecific::FR))
            && let Some(mut hex) = hexamer_to_int(seq, false)
        {
            let fwlimit = (seqlen as f64 - mean - 6.0).max(0.0) as isize;
            for j in 0..fwlimit {
                dbias5[hex] += contrib;
                if let Some(next) = update_hexamer(hex, seq[(j + 6) as usize], false) {
                    hex = next;
                } else {
                    break;
                }
            }
        }

        if strand_specific.is_none() || strand_specific == Some(StrandSpecific::RF) {
            let bwlimit = (mean - 6.0).max(0.0) as isize;
            let start = bwlimit as usize;
            if start + 6 <= seq.len()
                && let Some(mut hex) = hexamer_to_int(&seq[start..], true)
            {
                for j in bwlimit..(seqlen - 6) {
                    dbias5[hex] += contrib;
                    if (j as usize) < seq.len() - 6 {
                        if let Some(next) = update_hexamer(hex, seq[(j + 6) as usize], true) {
                            hex = next;
                        } else {
                            break;
                        }
                    }
                }
            }
        }
    }

    for val in &dbias5 {
        bias_alpha_norm += *val;
    }
    if bias_alpha_norm == 0.0 {
        return (eff_lens.to_vec(), dbias5);
    }

    let mut bias_lens = vec![0.0; transcript_seqs.len()];
    for i in 0..transcript_seqs.len() {
        if i >= means.len() || i >= alpha.len() || i >= eff_lens.len() {
            break;
        }
        let mean = means[i];
        let seq = &transcript_seqs[i];
        if seq.len() < 6 || (seq.len() as f64) < mean || alpha[i] < MIN_ALPHA {
            bias_lens[i] = eff_lens[i];
            continue;
        }
        let seqlen = seq.len() as isize;
        let mut efflen = 0.0;

        if (strand_specific.is_none() || strand_specific == Some(StrandSpecific::FR))
            && let Some(mut hex) = hexamer_to_int(seq, false)
        {
            let fwlimit = (seqlen as f64 - mean - 6.0).max(0.0) as isize;
            for j in 0..fwlimit {
                let denom = dbias5[hex];
                if denom > 0.0 {
                    efflen += bias_counts.counts[hex] as f64 / denom;
                }
                if let Some(next) = update_hexamer(hex, seq[(j + 6) as usize], false) {
                    hex = next;
                } else {
                    break;
                }
            }
        }

        if strand_specific.is_none() || strand_specific == Some(StrandSpecific::RF) {
            let bwlimit = (mean - 6.0).max(0.0) as isize;
            let start = bwlimit as usize;
            if start + 6 <= seq.len()
                && let Some(mut hex) = hexamer_to_int(&seq[start..], true)
            {
                for j in bwlimit..(seqlen - 6) {
                    let denom = dbias5[hex];
                    if denom > 0.0 {
                        efflen += bias_counts.counts[hex] as f64 / denom;
                    }
                    if (j as usize) < seq.len() - 6 {
                        if let Some(next) = update_hexamer(hex, seq[(j + 6) as usize], true) {
                            hex = next;
                        } else {
                            break;
                        }
                    }
                }
            }
        }

        if strand_specific.is_none() {
            efflen *= 0.5 * bias_alpha_norm / bias_data_norm;
        } else {
            efflen *= bias_alpha_norm / bias_data_norm;
        }

        if efflen > mean {
            bias_lens[i] = efflen;
        } else {
            bias_lens[i] = eff_lens[i];
        }
    }

    (bias_lens, dbias5)
}
