pub mod ec_compare;
pub mod ec_export;
pub mod ec_from_index;
pub mod ec_stats;
pub mod index_ec_dump;
pub mod index_info;
pub mod minimizer_bitmap_scan;
pub mod minimizer_lookup;
pub mod minimizer_mphf_keys;
pub mod minimizer_unitig_mphf_check;
pub mod pseudoalign;
pub mod quant;
pub mod trace_compare;
pub mod trace_reads;
pub mod transcript_lookup;
pub mod unitig_window;
pub mod version;

fn env_flag(name: &str) -> bool {
    std::env::var(name)
        .ok()
        .map(|v| {
            let v = v.trim();
            v == "1"
                || v.eq_ignore_ascii_case("true")
                || v.eq_ignore_ascii_case("yes")
                || v.eq_ignore_ascii_case("on")
        })
        .unwrap_or(false)
}

pub fn investigation_options_from_env() -> kallistors::pseudoalign::InvestigationOptions {
    kallistors::pseudoalign::InvestigationOptions {
        trace_fast_path: env_flag("KALLISTORS_TRACE_FAST_PATH"),
        bounded_incremental_scan: env_flag("KALLISTORS_FAST_DELTA_BOUNDED_INCREMENTAL_SCAN"),
        faster_probe_state: env_flag("KALLISTORS_FAST_DELTA_PROBE_STATE"),
        altered_hit_bookkeeping: env_flag("KALLISTORS_FAST_DELTA_HIT_BOOKKEEPING"),
        candidate_reuse: env_flag("KALLISTORS_FAST_DELTA_CANDIDATE_REUSE"),
    }
}
