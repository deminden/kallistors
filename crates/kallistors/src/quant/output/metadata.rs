//! Run metadata with upstream-compatible field order and JSON escaping.

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use crate::Result;
use crate::pseudoalign::MatchCacheAllocation;

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
        self.write_json_inner(path, None)
    }

    /// Include requested and effective match-cache allocation for reproducibility.
    pub fn write_json_with_match_cache(
        &self,
        path: &Path,
        cache: MatchCacheAllocation,
    ) -> Result<()> {
        self.write_json_inner(path, Some(cache))
    }

    fn write_json_inner(&self, path: &Path, cache: Option<MatchCacheAllocation>) -> Result<()> {
        let file = File::create(path)?;
        let mut writer = BufWriter::new(file);
        let mut entries = vec![
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
        if let Some(cache) = cache {
            entries.push(("match_cache", format!(
                "{{\"policy\":{},\"requested_bytes\":{},\"allocated_bytes\":{},\"bytes_per_worker\":{},\"workers\":{}}}",
                cache.policy().to_string().to_json(), cache.requested_bytes(),
                cache.total_bytes(), cache.bytes_per_worker(), cache.workers(),
            )));
        }
        writeln!(writer, "{{")?;
        for (idx, (key, value)) in entries.iter().enumerate() {
            let sep = if idx + 1 == entries.len() { "" } else { "," };
            writeln!(writer, "\"{}\": {}{}", escape_json(key), value, sep)?;
        }
        writeln!(writer, "}}")?;
        Ok(())
    }
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
