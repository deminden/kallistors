//! Prepared EC work, incoming edges and ordered weighted sums.

use super::*;

pub(super) struct EmMultiEc {
    pub(super) transcript_range: Range<usize>,
    pub(super) count: f64,
}

pub(super) struct EmWork {
    pub(super) singleton_alpha: Vec<f64>,
    pub(super) active_transcripts: Vec<usize>,
    pub(super) multi_ecs: Vec<EmMultiEc>,
    pub(super) multi_transcripts: Vec<usize>,
    pub(super) multi_weights: Vec<f64>,
}

impl EmWork {
    pub(super) fn compact_transcript_ids(&mut self) {
        let mut dense_ids = vec![0usize; self.singleton_alpha.len()];
        let mut singleton_alpha = Vec::with_capacity(self.active_transcripts.len());
        for (dense, &original) in self.active_transcripts.iter().enumerate() {
            dense_ids[original] = dense;
            singleton_alpha.push(self.singleton_alpha[original]);
        }
        for tr in &mut self.multi_transcripts {
            *tr = dense_ids[*tr];
        }
        self.singleton_alpha = singleton_alpha;
    }
}

pub(super) enum IncomingEcIds {
    Compact(Vec<u32>),
    Wide(Vec<usize>),
}

impl IncomingEcIds {
    pub(super) fn new(edges: usize, ecs: usize) -> Self {
        if ecs <= u32::MAX as usize {
            Self::Compact(vec![0; edges])
        } else {
            Self::Wide(vec![0; edges])
        }
    }

    pub(super) fn set(&mut self, slot: usize, ec: usize) {
        match self {
            Self::Compact(ids) => ids[slot] = ec as u32,
            Self::Wide(ids) => ids[slot] = ec,
        }
    }
}

pub(super) trait IncomingEcId: Copy {
    fn as_usize(self) -> usize;
}

impl IncomingEcId for u32 {
    #[inline]
    fn as_usize(self) -> usize {
        self as usize
    }
}

impl IncomingEcId for usize {
    #[inline]
    fn as_usize(self) -> usize {
        self
    }
}

#[inline]
pub(super) fn weighted_alpha_sum<I: IncomingEcId>(
    ids: &[I],
    weights: &[f64],
    alpha: &[f64],
) -> f64 {
    let mut sum = 0.0;
    for (&tr, &weight) in ids.iter().zip(weights) {
        sum += alpha[tr.as_usize()] * weight;
    }
    sum
}

pub(super) struct EcScales {
    pub(super) values: Vec<f64>,
    pub(super) present: Vec<u8>,
}

impl EcScales {
    pub(super) fn new(len: usize) -> Self {
        Self {
            values: vec![0.0; len],
            present: vec![0; len],
        }
    }
}

#[inline]
pub(super) fn sum_scaled_edges<I: IncomingEcId>(
    ids: &[I],
    weights: &[f64],
    scales: &EcScales,
    singleton: f64,
    alpha: f64,
) -> f64 {
    let mut sum = singleton;
    for (&ec, &weight) in ids.iter().zip(weights) {
        let ec = ec.as_usize();
        if scales.present[ec] != 0 {
            sum += weight * alpha * scales.values[ec];
        }
    }
    sum
}

// Each transcript's incoming edges retain their original EC/edge order.
// Independent transcript sums can then run concurrently without reordering f64 additions.
pub(super) struct EmIncoming {
    pub(super) offsets: Vec<usize>,
    pub(super) edge_indices: Vec<usize>,
    pub(super) ec_indices: IncomingEcIds,
    pub(super) weights: Vec<f64>,
}

impl EmIncoming {
    pub(super) fn new(work: &EmWork, track_edges: bool) -> Self {
        let mut offsets = vec![0; work.singleton_alpha.len() + 1];
        for &tr in &work.multi_transcripts {
            offsets[tr + 1] += 1;
        }
        for tr in 0..work.singleton_alpha.len() {
            offsets[tr + 1] += offsets[tr];
        }
        let mut cursors = offsets[..work.singleton_alpha.len()].to_vec();
        let edges = work.multi_transcripts.len();
        let mut incoming = Self {
            offsets,
            edge_indices: if track_edges {
                vec![0; edges]
            } else {
                Vec::new()
            },
            ec_indices: IncomingEcIds::new(edges, work.multi_ecs.len()),
            weights: vec![0.0; edges],
        };
        for (ec_id, ec) in work.multi_ecs.iter().enumerate() {
            for edge in ec.transcript_range {
                let tr = work.multi_transcripts[edge];
                let slot = cursors[tr];
                cursors[tr] += 1;
                if track_edges {
                    incoming.edge_indices[slot] = edge;
                }
                incoming.ec_indices.set(slot, ec_id);
                incoming.weights[slot] = work.multi_weights[edge];
            }
        }
        incoming
    }

    #[inline]
    pub(super) fn transcript_sum(
        &self,
        tr: usize,
        singleton: f64,
        alpha: f64,
        scales: &EcScales,
    ) -> f64 {
        let range = self.offsets[tr]..self.offsets[tr + 1];
        let weights = &self.weights[range.clone()];
        match &self.ec_indices {
            IncomingEcIds::Compact(ids) => {
                sum_scaled_edges(&ids[range], weights, scales, singleton, alpha)
            }
            IncomingEcIds::Wide(ids) => {
                sum_scaled_edges(&ids[range], weights, scales, singleton, alpha)
            }
        }
    }

    pub(super) fn refresh_weights(&mut self, work: &EmWork) {
        for (weight, &edge) in self.weights.iter_mut().zip(&self.edge_indices) {
            *weight = work.multi_weights[edge];
        }
    }
}

pub(super) fn prepare_em_work(ec_classes: &[Vec<u32>], counts: &[u32], eff_lens: &[f64]) -> EmWork {
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
                transcript_range: Range {
                    start,
                    end: multi_transcripts.len(),
                },
                count,
            });
        }
    }
    let mut active = vec![false; eff_lens.len()];
    for (tr, &count) in singleton_alpha.iter().enumerate() {
        active[tr] = count != 0.0;
    }
    for &tr in &multi_transcripts {
        active[tr] = true;
    }
    let active_transcripts = active
        .into_iter()
        .enumerate()
        .filter_map(|(tr, observed)| observed.then_some(tr))
        .collect();
    EmWork {
        singleton_alpha,
        active_transcripts,
        multi_ecs,
        multi_transcripts,
        multi_weights,
    }
}

pub(super) fn refresh_multi_weights(work: &mut EmWork, eff_lens: &[f64]) {
    for ec in &work.multi_ecs {
        for i in ec.transcript_range {
            let tr = work.multi_transcripts[i];
            let len = eff_lens.get(tr).copied().unwrap_or(1.0);
            work.multi_weights[i] = ec.count / len;
        }
    }
}

impl EmWork {
    #[allow(clippy::too_many_arguments)]
    pub(super) fn scalar_step(
        &self,
        ecs: Range<usize>,
        active: &[usize],
        alpha: &mut [f64],
        next_alpha: &mut [f64],
        scalar_products: &mut [f64],
    ) -> (usize, bool) {
        self.scalar_step_mode::<2>(ecs, active, alpha, next_alpha, scalar_products)
    }

    // MODE 2 counts every changed target, MODE 1 only proves any change,
    // MODE 0 omits a count that the independent witness proves unobservable.
    // All modes retain the exact same updates and fixed-point predicate.
    pub(super) fn scalar_step_mode<const MODE: u8>(
        &self,
        ecs: Range<usize>,
        active: &[usize],
        alpha: &mut [f64],
        next_alpha: &mut [f64],
        scalar_products: &mut [f64],
    ) -> (usize, bool) {
        self.scalar_targets::<MODE>(
            ecs,
            active.iter().copied(),
            alpha,
            next_alpha,
            scalar_products,
        )
    }

    pub(super) fn scalar_dense_step_mode<const MODE: u8>(
        &self,
        ecs: Range<usize>,
        active: Range<usize>,
        alpha: &mut [f64],
        next_alpha: &mut [f64],
        scalar_products: &mut [f64],
    ) -> (usize, bool) {
        self.scalar_targets::<MODE>(
            ecs,
            active.start..active.end,
            alpha,
            next_alpha,
            scalar_products,
        )
    }

    // Both monomorphizations visit exactly the same active target sequence.
    // Dense component callers pass their proven compact coordinates directly,
    // avoiding loads through the identity-ID vector. Row arithmetic is shared.
    fn scalar_targets<const MODE: u8>(
        &self,
        ecs: Range<usize>,
        active: impl Iterator<Item = usize> + Clone,
        alpha: &mut [f64],
        next_alpha: &mut [f64],
        scalar_products: &mut [f64],
    ) -> (usize, bool) {
        for tr in active.clone() {
            next_alpha[tr] = self.singleton_alpha[tr];
        }
        let ids = &self.multi_transcripts;

        let rows = &self.multi_ecs[ecs];
        let (pairs, tail) = rows.as_chunks::<2>();
        for [left, right] in pairs {
            let left_ids = &ids[left.transcript_range];
            let right_ids = &ids[right.transcript_range];
            let (left_products, rest) = scalar_products.split_at_mut(left_ids.len());
            let right_products = &mut rest[..right_ids.len()];
            let (left_denom, right_denom) = scalar_pair_products(
                left_ids,
                &self.multi_weights[left.transcript_range],
                left_products,
                right_ids,
                &self.multi_weights[right.transcript_range],
                right_products,
                alpha,
            );
            // Scatter order is exactly the original EC order, including when
            // these two rows contribute to the same transcript.
            scalar_scatter(left, left_ids, left_products, left_denom, next_alpha);
            scalar_scatter(right, right_ids, right_products, right_denom, next_alpha);
        }
        for ec in tail {
            let transcripts = &ids[ec.transcript_range];
            let weights = &self.multi_weights[ec.transcript_range];
            let products = &mut scalar_products[..transcripts.len()];
            let mut denom = 0.0;
            for ((&tr, &weight), product) in
                transcripts.iter().zip(weights).zip(products.iter_mut())
            {
                *product = alpha[tr.as_usize()] * weight;
                denom += *product;
            }
            scalar_scatter(ec, transcripts, products, denom, next_alpha);
        }
        let mut changes = 0;
        let mut identical = true;
        for tr in active {
            let next = next_alpha[tr];
            let old = alpha[tr];
            if MODE == 2 || (MODE == 1 && changes == 0) {
                changes += usize::from(
                    next > ALPHA_CHANGE_LIMIT && (next - old).abs() / next > ALPHA_CHANGE,
                );
            }
            identical &= next.to_bits() == old.to_bits();
            alpha[tr] = next;
        }
        (changes, identical)
    }
}

// Alpha is immutable until this component's complete scatter has finished.
// Interleave two independent denominator chains, preserving the original
// multiplication and addition sequence within each row.
#[inline]
#[allow(clippy::too_many_arguments)]
fn scalar_pair_products<I: IncomingEcId>(
    left_ids: &[I],
    left_weights: &[f64],
    left_products: &mut [f64],
    right_ids: &[I],
    right_weights: &[f64],
    right_products: &mut [f64],
    alpha: &[f64],
) -> (f64, f64) {
    let shared = left_ids.len().min(right_ids.len());
    let mut left_denom = 0.0;
    let mut right_denom = 0.0;
    for (((&left, &lw), lp), ((&right, &rw), rp)) in left_ids
        .iter()
        .zip(left_weights)
        .zip(&mut left_products[..shared])
        .zip(
            right_ids
                .iter()
                .zip(right_weights)
                .zip(&mut right_products[..shared]),
        )
    {
        *lp = alpha[left.as_usize()] * lw;
        left_denom += *lp;
        *rp = alpha[right.as_usize()] * rw;
        right_denom += *rp;
    }
    for ((&tr, &weight), product) in left_ids[shared..]
        .iter()
        .zip(&left_weights[shared..])
        .zip(&mut left_products[shared..])
    {
        *product = alpha[tr.as_usize()] * weight;
        left_denom += *product;
    }
    for ((&tr, &weight), product) in right_ids[shared..]
        .iter()
        .zip(&right_weights[shared..])
        .zip(&mut right_products[shared..])
    {
        *product = alpha[tr.as_usize()] * weight;
        right_denom += *product;
    }
    (left_denom, right_denom)
}

#[inline]
fn scalar_scatter<I: IncomingEcId>(
    ec: &EmMultiEc,
    ids: &[I],
    products: &[f64],
    denom: f64,
    next: &mut [f64],
) {
    if denom < TOLERANCE {
        return;
    }
    let scale = ec.count / denom;
    for (&tr, &product) in ids.iter().zip(products) {
        next[tr.as_usize()] += product * scale;
    }
}
