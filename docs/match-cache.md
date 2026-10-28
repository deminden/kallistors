# Match-cache sizing

`quant` and `pseudoalign` reuse exact k-mer match answers in a private table for
each mapping worker. Complete sequence keys and matching flags are checked on
every cache hit. Changing table capacity changes retention and memory use, not
mapping rules or equivalence-class membership.

## Automatic budgets

The default, `--match-cache-size auto`, uses these **aggregate** allowances:

| Mapping workers | Total allowance |
|---:|---:|
| 1–2 | 512 MiB |
| 3–8 | 1 GiB |
| 9–16 | 2 GiB |
| 17–32 | 4 GiB |
| 33–64 | 8 GiB |
| 65–128 | 16 GiB |
| Higher counts | Double at each power-of-two worker boundary |

Each worker receives an equal share. Its entry count is rounded down to a power
of two, so the actual total can be smaller than the allowance. For example, eight
workers use 128 MiB each; sixteen workers also use 128 MiB each. Three workers
use 256 MiB each, totaling 768 MiB within the 1 GiB allowance. A zero requested
thread count is treated as one worker, matching the mapping API. Above eight
workers, the nominal allowance is `128 MiB × next_power_of_two(workers)`.
Together with table rounding, this keeps each worker's effective table at
128 MiB from five workers upward on 64-bit targets. Arithmetic and addressability
are checked before allocating or starting workers.

These defaults are an allocation policy. They do not imply measured speed gains
at every thread count or on every dataset. The k-mer tables are not shared among
workers; the immutable index is shared. BUS processing retains its existing
cache policy and does not accept this option.

## Explicit budgets

Set a smaller or larger **total** allowance independently of worker count:

```sh
kallistors quant -i transcripts.idx -o counts -t 8 \
  --match-cache-size 4GiB reads_1.fastq.gz reads_2.fastq.gz

kallistors quant -i transcripts.idx -o counts -t 1 \
  --match-cache-size 256MiB reads_1.fastq.gz reads_2.fastq.gz
```

512 MiB is an automatic one-core default, not a required minimum. An explicit
4 GiB allowance gives eight workers 512 MiB each, or one worker 4 GiB. Explicit
allowances are not capped at the automatic tier sizes.

Sizes accept positive integers with binary units (`KiB`, `MiB`, `GiB`, `TiB`),
decimal units (`KB`, `MB`, `GB`, `TB`), or bytes (`B`, or no suffix). Units are
case-insensitive. For example, `1GiB` means 1,073,741,824 bytes; `1GB` means
1,000,000,000 bytes and can round to a smaller table. Fractions, overflow,
unrepresentable allocations and budgets smaller than one 32-byte entry per
worker are rejected.

Each worker's effective table is limited to 128 GiB by the cache's 32-bit slot index.
Larger resolved tables are rejected even on a machine with enough addressable memory.

On 32-bit targets, automatic allowances are capped at 1 GiB. Explicit sizes that
exceed the target's addressable allocation size are rejected rather than silently
replaced with a smaller allowance.

## Memory accounting

This allowance covers the exact k-mer match tables. **It is not a process RSS
limit:** the reference index, candidate selectors, input buffers, EC counts and
quantification arrays need additional memory. Larger tables can increase peak
RSS even when their allocation is released before later quantification stages.

The command reports the policy and effective per-worker and total allocation.
`quant` also records a `match_cache` object in `run_info.json`, with `policy`,
`requested_bytes`, `allocated_bytes`, `bytes_per_worker` and `workers`. Library
callers can set `PseudoalignOptions::match_cache_size` and inspect the resolved
layout with `MatchCacheSize::allocation` without allocating a table.
