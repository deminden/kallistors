use crate::{Error, Result};

pub(super) struct IndexWorkers {
    owned: Option<rayon::ThreadPool>,
    threads: usize,
}

impl IndexWorkers {
    pub(super) fn new(threads: usize) -> Result<Self> {
        let reuse = threads <= 1
            || (rayon::current_thread_index().is_some() && rayon::current_num_threads() == threads);
        let owned = if reuse {
            None
        } else {
            Some(
                rayon::ThreadPoolBuilder::new()
                    .num_threads(threads)
                    .build()
                    .map_err(|err| {
                        Error::InvalidFormat(format!("failed to create index worker pool: {err}"))
                    })?,
            )
        };
        Ok(Self { owned, threads })
    }

    pub(super) fn collect_unused<T>(&self, allocation: &[T]) {
        // SAFETY: this predicate accepts any live pointer and transfers no
        // ownership. System-allocator hosts must not initialize worker heaps.
        if allocation.is_empty()
            || !unsafe { libmimalloc_sys::mi_is_in_heap_region(allocation.as_ptr().cast()) }
        {
            return;
        }

        let collect = || {
            // SAFETY: mimalloc's thread-local collection reclaims only freed
            // blocks. No pointers, live values or ownership are transferred.
            unsafe {
                libmimalloc_sys::mi_collect(true);
            }
        };
        if self.threads > 1 {
            self.install(|| {
                rayon::broadcast(|_| collect());
            });
        } else {
            collect();
        }
    }

    pub(super) fn install<R: Send>(&self, operation: impl FnOnce() -> R + Send) -> R {
        match &self.owned {
            Some(pool) => pool.install(operation),
            None => operation(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn nested_stages_reuse_only_the_requested_parallel_budget() {
        IndexWorkers::new(3).unwrap().install(|| {
            assert_eq!(rayon::current_num_threads(), 3);
            let reused = IndexWorkers::new(3).unwrap();
            assert!(reused.owned.is_none());
            reused.install(|| assert_eq!(rayon::current_num_threads(), 3));
            IndexWorkers::new(2)
                .unwrap()
                .install(|| assert_eq!(rayon::current_num_threads(), 2));
            assert_eq!(rayon::current_num_threads(), 3);
        });
        assert!(IndexWorkers::new(1).unwrap().owned.is_none());
    }
}
