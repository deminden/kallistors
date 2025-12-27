use std::collections::HashMap;

use crate::bias::BiasCounts;

use super::EcCounts;

pub(crate) fn merge_ec_counts(
    target: &mut EcCounts,
    target_map: &mut HashMap<Vec<u32>, usize>,
    other: EcCounts,
) {
    target.reads_processed = target.reads_processed.saturating_add(other.reads_processed);
    target.reads_aligned = target.reads_aligned.saturating_add(other.reads_aligned);
    for (idx, ec) in other.ec_list.into_iter().enumerate() {
        let count = other.counts.get(idx).copied().unwrap_or(0);
        if let Some(&id) = target_map.get(&ec) {
            if let Some(slot) = target.counts.get_mut(id) {
                *slot = slot.saturating_add(count);
            }
        } else {
            let id = target.ec_list.len();
            target_map.insert(ec.clone(), id);
            target.ec_list.push(ec);
            target.counts.push(count);
        }
    }
    match (&mut target.bias, other.bias) {
        (Some(target_bias), Some(other_bias)) => merge_bias_counts(target_bias, &other_bias),
        (None, Some(other_bias)) => {
            let mut bias = BiasCounts::new();
            merge_bias_counts(&mut bias, &other_bias);
            target.bias = Some(bias);
        }
        _ => {}
    }
    match (
        &mut target.fragment_length_stats,
        other.fragment_length_stats,
    ) {
        (Some(target_stats), Some(other_stats)) => target_stats.merge(&other_stats),
        (None, Some(other_stats)) => target.fragment_length_stats = Some(other_stats),
        _ => {}
    }
    match (&mut target.fragment_length_hist, other.fragment_length_hist) {
        (Some(target_hist), Some(other_hist)) => merge_fragment_hist(target_hist, &other_hist),
        (None, Some(other_hist)) => target.fragment_length_hist = Some(other_hist),
        _ => {}
    }
}

fn merge_bias_counts(target: &mut BiasCounts, other: &BiasCounts) {
    let len = target.counts.len().min(other.counts.len());
    for i in 0..len {
        target.counts[i] = target.counts[i].saturating_add(other.counts[i]);
    }
    target.total = target.total.saturating_add(other.total);
}

fn merge_fragment_hist(target: &mut [u32], other: &[u32]) {
    let len = target.len().min(other.len());
    for i in 0..len {
        target[i] = target[i].saturating_add(other[i]);
    }
}
