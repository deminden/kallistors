pub mod bus;
pub mod ec_compare;
pub mod ec_export;
pub mod ec_from_index;
pub mod ec_stats;
pub mod index;
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

// Retain nested blocks for positional filters, diagnostics and alternate scanners.
// Ordinary counting uses the same flat EC values and jump bounds without a duplicate graph.
// Packed graph bases retain the same positional annotations for ordinary
// fragment and strand filters. Diagnostics and alternate scanners stay full.
fn use_packed_positions(
    load_positional_info: bool,
    options: kallistors::pseudoalign::PseudoalignOptions,
    debug: bool,
) -> bool {
    (load_positional_info || options.strand_specific.is_some())
        && !debug
        && !options.bias
        && !options.kallisto_direct_kmer
        && !options.kallisto_fallback
        && !options.do_union
        && !options.dfk_onlist
        && !options.no_jump
        && !options.kallisto_enum
        && !options.kallisto_strict
        && !options.kallisto_local_fallback
        && !options.discard_special_only
        && !options.kallisto_bifrost_find
        && !options.investigation.any_enabled()
}
fn build_mapping_index(
    path: &std::path::Path,
    load_positional_info: bool,
    threads: usize,
    options: kallistors::pseudoalign::PseudoalignOptions,
    debug: bool,
) -> kallistors::Result<kallistors::pseudoalign::BifrostIndex> {
    if use_packed_positions(load_positional_info, options, debug) {
        return kallistors::pseudoalign::build_bifrost_index_positional_packed_threaded(
            path, threads, false,
        );
    }
    build_mapping_index_impl(path, load_positional_info, threads, options, debug, false)
}
fn build_quant_index(
    path: &std::path::Path,
    load_positional_info: bool,
    threads: usize,
    options: kallistors::pseudoalign::PseudoalignOptions,
) -> kallistors::Result<kallistors::pseudoalign::BifrostIndex> {
    if use_packed_positions(load_positional_info, options, false) {
        return kallistors::pseudoalign::build_bifrost_index_positional_packed_threaded(
            path,
            threads,
            !options.bias,
        );
    }
    build_mapping_index_impl(path, load_positional_info, threads, options, false, true)
}
fn build_mapping_index_impl(
    path: &std::path::Path,
    load_positional_info: bool,
    threads: usize,
    options: kallistors::pseudoalign::PseudoalignOptions,
    debug: bool,
    defer_names: bool,
) -> kallistors::Result<kallistors::pseudoalign::BifrostIndex> {
    use kallistors::pseudoalign::*;
    let load_positional_info = load_positional_info || options.strand_specific.is_some();
    if options.kallisto_direct_kmer {
        build_bifrost_index_with_kmer_pos(path, load_positional_info)
    } else if options.kallisto_fallback {
        build_bifrost_index_with_positions_and_kmer(path, load_positional_info)
    } else if load_positional_info
        || debug
        || options.strand_specific.is_some()
        || options.do_union
        || options.dfk_onlist
        || options.no_jump
        || options.kallisto_enum
        || options.kallisto_strict
        || options.kallisto_local_fallback
        || options.discard_special_only
        || options.kallisto_bifrost_find
        || options.investigation.any_enabled()
    {
        build_bifrost_index_with_positions_threaded(path, load_positional_info, threads)
    } else if options.bias {
        build_bifrost_index_quant_threaded(path, false, threads)
    } else if defer_names {
        build_bifrost_index_quant_deferred_threaded(path, threads)
    } else {
        build_bifrost_index_quant_packed_threaded(path, threads)
    }
}
