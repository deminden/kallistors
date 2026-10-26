//! Exact inverse-CDF resampling with a bounded shared search directory.

use super::SplitMix64;
use strength_reduce::StrengthReducedU64;

pub(super) struct SamplingPlan {
    total: u64,
    cumulative: Vec<u64>,
    reduced: Option<StrengthReducedU64>,
    zone: u64,
    shift: u32,
    ranges: Vec<(usize, usize)>,
    inverse: Vec<u32>,
}

impl SamplingPlan {
    pub(super) fn new(counts: &[u32]) -> Self {
        // Retain the original total and saturating cumulative equations.
        let total = counts.iter().map(|&value| u64::from(value)).sum::<u64>();
        let mut running = 0u64;
        let cumulative = counts
            .iter()
            .map(|&value| {
                running = running.saturating_add(u64::from(value));
                running
            })
            .collect();
        Self::from_cumulative(total, cumulative)
    }

    fn from_cumulative(total: u64, cumulative: Vec<u64>) -> Self {
        let shift = (64 - total.saturating_sub(1).leading_zeros()).saturating_sub(12);
        let ranges = if total > 0 && cumulative.len() >= 128 {
            let buckets = ((total - 1) >> shift) + 1;
            (0..buckets)
                .map(|bucket| {
                    let low = bucket << shift;
                    let high = low.saturating_add((1u64 << shift) - 1).min(total - 1);
                    let start = cumulative.partition_point(|&bound| bound <= low);
                    let end = cumulative.partition_point(|&bound| bound <= high);
                    (start, end)
                })
                .collect()
        } else {
            Vec::new()
        };
        Self {
            inverse: Vec::new(),
            total,
            cumulative,
            reduced: (total > 0).then(|| StrengthReducedU64::new(total)),
            zone: if total > 0 {
                u64::MAX - u64::MAX % total
            } else {
                0
            },
            shift,
            ranges,
        }
    }

    // The serial caller can amortize one bounded exact table across samples.
    // Requested multiworker pools keep the existing compact search directory.
    pub(super) fn with_exact_inverse(mut self) -> Self {
        let total = self.total;
        let cumulative = &self.cumulative;
        // Exact inverse CDF, shared by all bootstrap workers. Keep larger
        // libraries on the bounded directory; at most 64 MiB is allocated.
        let inverse = if total <= (64 * 1024 * 1024 / 4)
            && total > 0
            && cumulative.len() >= 128
            && cumulative.len() <= u32::MAX as usize
        {
            let mut table = Vec::with_capacity(total as usize);
            for (index, &bound) in cumulative.iter().enumerate() {
                table.resize(bound.min(total) as usize, index as u32);
            }
            // Synthetic/incomplete cumulative tables can leave draws beyond
            // the last bound. The original search returns len in that case.
            table.resize(total as usize, cumulative.len() as u32);
            table
        } else {
            Vec::new()
        };
        self.inverse = inverse;
        self
    }

    fn draw(&self, rng: &mut SplitMix64) -> u64 {
        let reduced = self.reduced.expect("positive sampling total");
        loop {
            let value = rng.next_u64();
            if value < self.zone {
                return value % reduced;
            }
        }
    }

    fn locate(&self, draw: u64) -> usize {
        if !self.inverse.is_empty() {
            self.inverse[draw as usize] as usize
        } else if self.ranges.is_empty() {
            self.cumulative.partition_point(|&bound| bound <= draw)
        } else {
            let (start, end) = self.ranges[(draw >> self.shift) as usize];
            // All preceding bounds are <= bucket low; the excluded end bound
            // (if present) is > bucket high. Thus returning end is exact too.
            start + self.cumulative[start..end].partition_point(|&bound| bound <= draw)
        }
    }

    // Commit independently addressed SplitMix chunks only if every candidate
    // passes the original rejection zone. A rejection replays the serial
    // sampler from the untouched RNG, retaining its exact draws and state.
    #[inline(never)]
    pub(super) fn sample_with_pool(
        &self,
        rng: &mut SplitMix64,
        pool: &rayon::ThreadPool,
    ) -> Vec<u32> {
        use rayon::prelude::*;
        let workers = pool.current_num_threads();
        let bytes = self
            .cumulative
            .len()
            .saturating_mul(4)
            .saturating_mul(workers);
        if workers <= 1 || self.total < 262_144 || bytes > 64 * 1024 * 1024 {
            return self.sample(rng);
        }
        let reduced = self.reduced.expect("positive sampling total");
        let initial = rng.state;
        let histograms: Option<Vec<_>> = pool.install(|| {
            (0..workers)
                .into_par_iter()
                .map(|worker| {
                    let base = self.total / workers as u64;
                    let remainder = self.total % workers as u64;
                    let start = base * worker as u64 + (worker as u64).min(remainder);
                    let count = base + u64::from((worker as u64) < remainder);
                    let mut local = SplitMix64::new(initial);
                    local.advance(start);
                    let mut sampled = vec![0u32; self.cumulative.len()];
                    for _ in 0..count {
                        let value = local.next_u64();
                        if value >= self.zone {
                            return None;
                        }
                        let index = self.locate(value % reduced);
                        if let Some(slot) = sampled.get_mut(index) {
                            *slot = slot.saturating_add(1);
                        }
                    }
                    Some(sampled)
                })
                .collect()
        });
        let Some(mut histograms) = histograms else {
            return self.sample(rng);
        };
        let mut sampled = histograms.pop().expect("positive worker count");
        for histogram in histograms {
            for (slot, count) in sampled.iter_mut().zip(histogram) {
                *slot = slot.saturating_add(count);
            }
        }
        rng.advance(self.total);
        sampled
    }

    pub(super) fn sample(&self, rng: &mut SplitMix64) -> Vec<u32> {
        let mut sampled = vec![0u32; self.cumulative.len()];
        // RNG rejection order remains serial and identical. Only the
        // independent integer lookups are exposed together to the compiler.
        // Histogram increments remain in draw order, including saturation.
        for _ in 0..self.total / 4 {
            let draws: [u64; 4] = std::array::from_fn(|_| self.draw(rng));
            let indices = draws.map(|draw| self.locate(draw));
            for index in indices {
                if let Some(slot) = sampled.get_mut(index) {
                    *slot = slot.saturating_add(1);
                }
            }
        }
        for _ in 0..self.total % 4 {
            let index = self.locate(self.draw(rng));
            if let Some(slot) = sampled.get_mut(index) {
                *slot = slot.saturating_add(1);
            }
        }
        sampled
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn original_draw(rng: &mut SplitMix64, upper: u64) -> u64 {
        let zone = u64::MAX - u64::MAX % upper;
        loop {
            let value = rng.next_u64();
            if value < zone {
                return value % upper;
            }
        }
    }

    #[test]
    fn directory_preserves_all_small_draws_and_u64_tail_boundaries() {
        let mut exercised = [false; 2];
        for length in [0, 1, 127, 128, 1024] {
            for total in [0, 1, 2, 4095, 4096, 4097, 1 << 31, (1 << 63) + 1, u64::MAX] {
                let cumulative: Vec<_> = (0..length)
                    .map(|i| ((i as u128 * total as u128) / length as u128) as u64)
                    .collect();
                for dense in [false, true] {
                    let plan = SamplingPlan::from_cumulative(total, cumulative.clone());
                    let plan = if dense {
                        plan.with_exact_inverse()
                    } else {
                        plan
                    };
                    assert!(plan.inverse.len() <= 64 * 1024 * 1024 / 4);
                    assert!(plan.ranges.len() <= 4096);
                    if total == 0 {
                        continue;
                    }
                    exercised[usize::from(!plan.ranges.is_empty())] = true;
                    let mut draws = vec![0, total / 2, total - 1];
                    if total <= 8192 {
                        draws.extend(0..total);
                    }
                    for &bound in &plan.cumulative {
                        for value in [bound.saturating_sub(1), bound, bound.saturating_add(1)] {
                            if value < total {
                                draws.push(value);
                            }
                        }
                    }
                    for (bucket, _) in plan.ranges.iter().enumerate() {
                        let low = (bucket as u64) << plan.shift;
                        for value in [low.saturating_sub(1), low, low.saturating_add(1)] {
                            if value < total {
                                draws.push(value);
                            }
                        }
                    }
                    for draw in draws {
                        assert_eq!(
                            plan.locate(draw),
                            plan.cumulative.partition_point(|&bound| bound <= draw)
                        );
                    }
                }
            }
        }
        assert_eq!(exercised, [true, true]);
    }

    #[test]
    fn reduced_draw_and_sample_preserve_rejections_rng_state_and_zero_counts() {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(3)
            .build()
            .unwrap();
        for counts in [vec![65_536, 0, 131_072, 65_536], vec![1_024; 257]] {
            let plan = SamplingPlan::new(&counts);
            for seed in [0, 42, u64::MAX, 0x31628af67b2131ab] {
                let mut expected_rng = SplitMix64::new(seed);
                let mut actual_rng = SplitMix64::new(seed);
                assert_eq!(
                    plan.sample_with_pool(&mut actual_rng, &pool),
                    plan.sample(&mut expected_rng)
                );
                assert_eq!(actual_rng.state, expected_rng.state);
                assert_eq!(actual_rng.next_u64(), expected_rng.next_u64());
            }
        }
        // This fixed state emits MAX first, which the original zone rejects
        // even when upper==1 or MAX. Keep the exact boundary observable.
        assert_eq!(SplitMix64::new(0x31628af67b2131ab).next_u64(), u64::MAX);
        for upper in [1, 2, 3, 17, 4097, (1 << 32) + 3, (1 << 63) + 1, u64::MAX] {
            let plan = SamplingPlan::from_cumulative(upper, vec![upper]);
            for seed in [0, 42, u64::MAX, 0x31628af67b2131ab] {
                let mut a = SplitMix64::new(seed);
                let mut b = SplitMix64::new(seed);
                for _ in 0..1000 {
                    assert_eq!(plan.draw(&mut a), original_draw(&mut b, upper));
                    assert_eq!(a.state, b.state);
                }
            }
        }
        for length in [0, 1, 127, 128, 1024] {
            let counts: Vec<_> = (0..length)
                .map(|i| if i % 3 == 0 { 0 } else { (i % 7) as u32 })
                .collect();
            let plan = SamplingPlan::new(&counts).with_exact_inverse();
            for seed in [0, 42, u64::MAX] {
                let mut a = SplitMix64::new(seed);
                let mut b = SplitMix64::new(seed);
                let actual = plan.sample(&mut a);
                let mut expected = vec![0u32; length];
                // Independent original cumulative/rejection/search equations.
                let total: u64 = counts.iter().map(|&x| u64::from(x)).sum();
                let mut running = 0u64;
                let mut cumulative = Vec::new();
                for &count in &counts {
                    running = running.saturating_add(u64::from(count));
                    cumulative.push(running);
                }
                for _ in 0..total {
                    let draw = original_draw(&mut b, total);
                    let index = cumulative.partition_point(|&bound| bound <= draw);
                    if let Some(slot) = expected.get_mut(index) {
                        *slot = slot.saturating_add(1);
                    }
                }
                assert_eq!(actual, expected);
                assert_eq!(a.state, b.state);
            }
        }
    }
}
