use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Mutex, OnceLock};
use std::time::{Duration, Instant};

const STAGE_COUNT: usize = 8;

#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum Stage {
    IndexHeaderParse = 0,
    GraphDecode = 1,
    MinimizerCountPass = 2,
    MinimizerFillPass = 3,
    FastqReadDecompress = 4,
    Pseudoalign = 5,
    EcMerge = 6,
    Em = 7,
}

impl Stage {
    pub const ALL: [Stage; STAGE_COUNT] = [
        Stage::IndexHeaderParse,
        Stage::GraphDecode,
        Stage::MinimizerCountPass,
        Stage::MinimizerFillPass,
        Stage::FastqReadDecompress,
        Stage::Pseudoalign,
        Stage::EcMerge,
        Stage::Em,
    ];

    pub const fn label(self) -> &'static str {
        match self {
            Stage::IndexHeaderParse => "index_header_parse",
            Stage::GraphDecode => "graph_decode",
            Stage::MinimizerCountPass => "minimizer_count_pass",
            Stage::MinimizerFillPass => "minimizer_fill_pass",
            Stage::FastqReadDecompress => "fastq_read_decompress",
            Stage::Pseudoalign => "pseudoalign",
            Stage::EcMerge => "ec_merge",
            Stage::Em => "em",
        }
    }
}

#[derive(Clone, Debug)]
pub struct StageTiming {
    pub stage: Stage,
    pub duration: Duration,
}

static ENABLED: AtomicBool = AtomicBool::new(false);
static TOTALS: OnceLock<Mutex<[u128; STAGE_COUNT]>> = OnceLock::new();

fn totals() -> &'static Mutex<[u128; STAGE_COUNT]> {
    TOTALS.get_or_init(|| Mutex::new([0; STAGE_COUNT]))
}

pub fn set_enabled(enabled: bool) {
    ENABLED.store(enabled, Ordering::Relaxed);
}

pub fn is_enabled() -> bool {
    ENABLED.load(Ordering::Relaxed)
}

pub fn reset() {
    let mut guard = totals().lock().expect("timing totals poisoned");
    *guard = [0; STAGE_COUNT];
}

pub fn add_duration(stage: Stage, duration: Duration) {
    if !is_enabled() {
        return;
    }
    let mut guard = totals().lock().expect("timing totals poisoned");
    guard[stage as usize] += duration.as_nanos();
}

pub fn snapshot() -> Vec<StageTiming> {
    let guard = totals().lock().expect("timing totals poisoned");
    Stage::ALL
        .into_iter()
        .map(|stage| {
            let nanos = guard[stage as usize].min(u64::MAX as u128) as u64;
            StageTiming {
                stage,
                duration: Duration::from_nanos(nanos),
            }
        })
        .collect()
}

pub struct StageGuard {
    stage: Stage,
    start: Instant,
}

impl StageGuard {
    pub fn new(stage: Stage) -> Self {
        Self {
            stage,
            start: Instant::now(),
        }
    }
}

impl Drop for StageGuard {
    fn drop(&mut self) {
        add_duration(self.stage, self.start.elapsed());
    }
}

pub fn scoped(stage: Stage) -> Option<StageGuard> {
    is_enabled().then(|| StageGuard::new(stage))
}
