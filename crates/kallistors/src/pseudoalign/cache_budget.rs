//! Aggregate budgets for the private exact-match tables used by bulk workers.

use std::fmt;
use std::str::FromStr;

use crate::{Error, Result};

pub(crate) const MATCH_CACHE_ENTRY_BYTES: usize = 32;
const MIB: u64 = 1 << 20;

/// Total exact k-mer match-cache allowance across bulk mapping workers.
///
/// This excludes the index, other caches and quantification storage. Each worker
/// owns a private power-of-two table; rounding can leave some allowance unused.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
#[repr(transparent)]
pub struct MatchCacheSize(u64);

/// Resolved allocation without allocating or touching the cache itself.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MatchCacheAllocation {
    policy: MatchCacheSize,
    requested_bytes: usize,
    workers: usize,
    entries_per_worker: usize,
}

impl MatchCacheSize {
    /// 512 MiB for 1–2 workers, 1 GiB for 3–8, then a doubling aggregate tier:
    /// 2 GiB for 9–16, 4 GiB for 17–32, 8 GiB for 33–64, and so on.
    /// Automatic budgets on 32-bit targets are capped at 1 GiB.
    pub const AUTO: Self = Self(0);

    /// Explicit aggregate byte allowance, independent of worker count.
    pub fn from_bytes(bytes: u64) -> Result<Self> {
        if bytes < MATCH_CACHE_ENTRY_BYTES as u64 {
            return Err(Error::InvalidFormat(
                "match-cache budget must allow at least one 32-byte entry".into(),
            ));
        }
        addressable_bytes(bytes)?;
        Ok(Self(bytes))
    }

    /// Resolve once before mapping. Zero requested threads means one worker.
    pub fn allocation(self, workers: usize) -> Result<MatchCacheAllocation> {
        let workers = workers.max(1);
        let bytes = if self.0 == 0 {
            if workers <= 2 {
                512 * MIB
            } else if workers <= 8 || usize::BITS < 64 {
                1024 * MIB
            } else {
                (workers as u64)
                    .checked_next_power_of_two()
                    .and_then(|workers| workers.checked_mul(128 * MIB))
                    .ok_or_else(|| {
                        Error::InvalidFormat(
                            "automatic match-cache budget overflows for this worker count".into(),
                        )
                    })?
            }
        } else {
            self.0
        };
        let requested_bytes = addressable_bytes(bytes)?;
        let limit = requested_bytes / MATCH_CACHE_ENTRY_BYTES / workers;
        if limit == 0 {
            return Err(Error::InvalidFormat(format!(
                "match-cache budget must allow at least {MATCH_CACHE_ENTRY_BYTES} bytes per mapping worker ({workers} workers)"
            )));
        }
        let entries_per_worker = 1usize << (usize::BITS - 1 - limit.leading_zeros());
        // The slot hash supplies 32 bits; a larger table would contain entries
        // that cannot be addressed. Reject it instead of silently wasting RAM.
        if entries_per_worker as u64 > 1u64 << 32 {
            return Err(Error::InvalidFormat(
                "match-cache table exceeds the supported 128 GiB per worker".into(),
            ));
        }
        Ok(MatchCacheAllocation {
            policy: self,
            requested_bytes,
            workers,
            entries_per_worker,
        })
    }
}

fn addressable_bytes(bytes: u64) -> Result<usize> {
    usize::try_from(bytes)
        .ok()
        .filter(|&value| value <= isize::MAX as usize)
        .ok_or_else(|| {
            Error::InvalidFormat(format!(
                "match-cache budget {bytes} bytes exceeds this target's addressable allocation size"
            ))
        })
}

impl MatchCacheAllocation {
    pub fn policy(self) -> MatchCacheSize {
        self.policy
    }

    pub fn requested_bytes(self) -> usize {
        self.requested_bytes
    }

    pub fn workers(self) -> usize {
        self.workers
    }

    pub fn bytes_per_worker(self) -> usize {
        self.entries_per_worker * MATCH_CACHE_ENTRY_BYTES
    }

    pub fn total_bytes(self) -> usize {
        self.bytes_per_worker() * self.workers
    }

    pub(crate) fn entries_per_worker(self) -> usize {
        self.entries_per_worker
    }
}

impl fmt::Display for MatchCacheSize {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.0 == 0 {
            f.write_str("auto")
        } else {
            write!(f, "{}B", self.0)
        }
    }
}

impl FromStr for MatchCacheSize {
    type Err = String;

    fn from_str(value: &str) -> std::result::Result<Self, Self::Err> {
        let value = value.trim();
        if value.eq_ignore_ascii_case("auto") {
            return Ok(Self::AUTO);
        }
        let split = value
            .find(|c: char| !c.is_ascii_digit())
            .unwrap_or(value.len());
        let (number, unit) = value.split_at(split);
        let number: u64 = number.parse().map_err(|_| {
            "expected auto or a positive integer size such as 256MiB or 4GiB".to_owned()
        })?;
        let multiplier = match unit.trim().to_ascii_lowercase().as_str() {
            "" | "b" => 1,
            "kib" => 1 << 10,
            "mib" => 1 << 20,
            "gib" => 1 << 30,
            "tib" => 1 << 40,
            "kb" => 1_000,
            "mb" => 1_000_000,
            "gb" => 1_000_000_000,
            "tb" => 1_000_000_000_000,
            _ => {
                return Err(
                    "unknown size unit; use B, KiB, MiB, GiB, TiB, KB, MB, GB or TB".into(),
                );
            }
        };
        let bytes = number
            .checked_mul(multiplier)
            .filter(|&bytes| bytes >= MATCH_CACHE_ENTRY_BYTES as u64)
            .ok_or_else(|| {
                "size overflows or is smaller than one 32-byte cache entry".to_owned()
            })?;
        Self::from_bytes(bytes).map_err(|err| err.to_string())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn budgets_are_aggregate_rounded_and_do_not_allocate() {
        assert_eq!(std::mem::size_of::<MatchCacheSize>(), 8);
        assert!(MatchCacheSize::from_bytes(0).is_err());
        assert!(MatchCacheSize::from_bytes(31).is_err());
        for (workers, requested_mib, per_worker_mib, total_mib) in [
            (0, 512, 512, 512),
            (1, 512, 512, 512),
            (2, 512, 256, 512),
            (3, 1024, 256, 768),
            (4, 1024, 256, 1024),
            (8, 1024, 128, 1024),
        ] {
            let allocation = MatchCacheSize::AUTO.allocation(workers).unwrap();
            assert_eq!(allocation.requested_bytes() as u64, requested_mib * MIB);
            assert_eq!(allocation.bytes_per_worker() as u64, per_worker_mib * MIB);
            assert_eq!(allocation.total_bytes() as u64, total_mib * MIB);
            assert!(allocation.entries_per_worker().is_power_of_two());
            assert!(allocation.total_bytes() <= allocation.requested_bytes());
        }
        if usize::BITS >= 64 {
            for (workers, requested_mib, per_worker_mib, total_mib) in [
                (9, 2048, 128, 1152),
                (16, 2048, 128, 2048),
                (17, 4096, 128, 2176),
                (32, 4096, 128, 4096),
                (33, 8192, 128, 4224),
                (64, 8192, 128, 8192),
                (65, 16384, 128, 8320),
                (128, 16384, 128, 16384),
                (129, 32768, 128, 16512),
            ] {
                let allocation = MatchCacheSize::AUTO.allocation(workers).unwrap();
                assert_eq!(allocation.requested_bytes() as u64, requested_mib * MIB);
                assert_eq!(allocation.bytes_per_worker() as u64, per_worker_mib * MIB);
                assert_eq!(allocation.total_bytes() as u64, total_mib * MIB);
            }
            let override_size = MatchCacheSize::from_bytes(4 << 30).unwrap();
            assert_eq!(
                override_size.allocation(1).unwrap().total_bytes() as u64,
                4u64 << 30
            );
            assert_eq!(
                override_size.allocation(8).unwrap().bytes_per_worker() as u64,
                512u64 << 20
            );
            // The aggregate can exceed the per-worker hash-width bound when
            // enough workers divide it; parsing must not assume one worker.
            let large: MatchCacheSize = "1TiB".parse().unwrap();
            assert!(large.allocation(1).is_err());
            assert_eq!(
                large.allocation(8).unwrap().bytes_per_worker() as u64,
                128u64 << 30
            );
        } else {
            assert_eq!(
                MatchCacheSize::AUTO.allocation(16).unwrap().total_bytes(),
                1usize << 30
            );
            assert!(MatchCacheSize::from_bytes(4 << 30).is_err());
        }
        assert!(
            MatchCacheSize::from_bytes(32)
                .unwrap()
                .allocation(2)
                .is_err()
        );
        assert!(MatchCacheSize::from_bytes(u64::MAX).is_err());
        assert!(MatchCacheSize::AUTO.allocation(usize::MAX).is_err());
    }

    #[test]
    fn sizes_parse_without_overflow_or_ambiguous_multipliers() {
        for (text, expected) in [
            ("auto", MatchCacheSize::AUTO),
            (" AUTO ", MatchCacheSize::AUTO),
            ("256MiB", MatchCacheSize::from_bytes(256 << 20).unwrap()),
            ("1gib", MatchCacheSize::from_bytes(1 << 30).unwrap()),
            ("64 MiB", MatchCacheSize::from_bytes(64 << 20).unwrap()),
            ("256MB", MatchCacheSize::from_bytes(256_000_000).unwrap()),
            ("4096", MatchCacheSize::from_bytes(4096).unwrap()),
            ("32B", MatchCacheSize::from_bytes(32).unwrap()),
        ] {
            assert_eq!(text.parse::<MatchCacheSize>().unwrap(), expected);
            assert_eq!(
                expected.to_string().parse::<MatchCacheSize>().unwrap(),
                expected
            );
        }
        for text in [
            "",
            "0",
            "31B",
            "-1GiB",
            "1.5GiB",
            "1watts",
            "18446744073709551615GiB",
        ] {
            assert!(text.parse::<MatchCacheSize>().is_err(), "{text}");
        }
        if usize::BITS >= 64 {
            assert_eq!(
                "4GiB".parse::<MatchCacheSize>().unwrap(),
                MatchCacheSize::from_bytes(4 << 30).unwrap()
            );
        } else {
            assert!("4GiB".parse::<MatchCacheSize>().is_err());
        }
    }
}
