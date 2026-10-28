//! Exact grouped transcript positions without an allocation per transcript.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub(super) struct BlockPositions {
    ends: Vec<usize>,
    values: Vec<u32>,
}
impl BlockPositions {
    pub(super) fn with_capacity(groups: usize, values: usize) -> Self {
        Self {
            ends: Vec::with_capacity(groups),
            values: Vec::with_capacity(values),
        }
    }
    pub(super) fn len(&self) -> usize {
        self.ends.len()
    }
    pub(super) fn iter(&self) -> impl Iterator<Item = &[u32]> {
        let mut start = 0;
        self.ends.iter().map(move |&end| {
            let values = &self.values[start..end];
            start = end;
            values
        })
    }
    pub(super) fn push(&mut self, values: impl IntoIterator<Item = u32>) {
        self.start_group();
        for value in values {
            self.push_last_unique(value);
        }
    }
    pub(super) fn start_group(&mut self) {
        self.ends.push(self.values.len());
    }
    pub(super) fn push_last_unique(&mut self, value: u32) {
        // Input for each group is sorted, so adjacent comparison is sufficient.
        // Compare within this group only: neighboring transcripts may share a position.
        let group = self.ends.len() - 1;
        let start = if group == 0 { 0 } else { self.ends[group - 1] };
        if self.values.len() == start || self.values.last() != Some(&value) {
            self.values.push(value);
            self.ends[group] = self.values.len();
        }
    }
    pub(super) fn clear(&mut self) {
        self.ends.clear();
        self.values.clear();
    }
}
