//! Index inspection and read traces, kept separate from mapping kernels.

mod ec_dump;
mod minimizer;
mod tracing;

pub use ec_dump::{write_index_ec_dump, write_index_ec_dump_filtered};
pub use minimizer::debug_minimizer_lookup;
pub(super) use tracing::placement_for_match;
pub(crate) use tracing::trace_result_from_debug_state;
pub use tracing::{trace_read_bifrost, trace_read_bifrost_with_hits};
