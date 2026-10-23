//! Focused pseudoalignment core modules.
//!
//! The read scanners are still large because they mirror kallisto's control flow,
//! but support code is split by role so changes land in smaller files.

mod count;
mod fast;
mod filters;
mod local;
mod matcher;
mod read;
mod runtime;

pub(super) use filters::{
    apply_strand_filter, block_index_for_position, filter_ec_by_fragment,
    filter_ec_by_fragment_cow, filter_ec_by_fragment_in_place, filter_ec_by_strand_cow,
};
pub use local::local_kmer_hits;
pub(super) use matcher::{special_unitig_for_kmer, special_unitig_for_kmer_raw};
pub(super) use read::{
    ec_for_read_bifrost, ec_for_read_bifrost_counting, ec_for_read_bifrost_counting_borrowed,
    ec_for_read_bifrost_single, ec_for_read_bifrost_single_borrowed,
};
pub(crate) use runtime::reset_thread_local_caches;

pub(crate) use filters::filter_ec_by_fragment_and_strand_cow;
