# physmerge C executable: performance and validation report

Machine: Apple M4 Pro, 24 GB RAM, macOS 26.6.2, Apple clang 17 (`-O2`),
R 4.5.3 with data.table. Timings are the median of three runs; memory is
`maximum resident set size` from `/usr/bin/time -l`.

"R pipeline" = `read_sumstat()` → `physical_merge()` → `annotate_blocks()` →
`export_snp_list()`. "C" = one `physmerge` invocation doing all four.

## 1. Speed and memory

| Dataset | Rows | File | R wall | R CPU | R peak RSS | C wall | C peak RSS | Speed-up | Memory ratio |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| synthetic, 22 chr | 1 M | 74 MB | 4.11 s | 4.1 s | 996 MB | **0.14 s** | **2.4 MB** | 29× | 415× |
| synthetic, 22 chr | 5 M | 371 MB | 21.04 s | 21.1 s | 3 839 MB | **0.73 s** | **2.4 MB** | 29× | 1 600× |
| synthetic, 22 chr | 10 M | 743 MB | 44.50 s | 44.6 s | 7 271 MB | **1.49 s** | **2.4 MB** | 30× | 3 030× |
| real: chr22 HbA1c PLINK2 (16 columns, TEST filter) | 18.8 M | 1.85 GB | 10.89 s | 36.5 s | 4 502 MB | **2.81 s** | **2.4 MB** | 3.9× wall, 12.9× CPU | 1 876× |

Throughput of the C tool is flat at ~6.7 million rows per second in all four
cases; it is bound by line parsing, not by the merge itself.

Two observations explain the shape of the table:

- **R memory grows linearly with the file; C memory does not.** The R pipeline
  materializes the whole table plus the `position`/`value` copies appended by
  `read_sumstat()`, so 10 M SNPs cost 7.3 GB. The C tool holds one line and one
  pending block, so its footprint is 2.4 MB (dominated by the 1 MB read buffer)
  whether the input is 74 MB or 1.85 GB.
- **The R bottleneck is the scan loop, not I/O.** On the real file the TEST
  filter reduces 18.8 M rows to 1.25 M before merging, so R spends most of its
  wall time inside multithreaded `fread` and finishes in 10.9 s wall but 36.5 s
  CPU. On the synthetic files all rows survive to the merge, and the pure-R
  `for` loop in `.physical_merge_single()` dominates; that is where the ~30×
  gap comes from.

## 2. Correctness

| Suite | Cases | Result |
|---|---:|---|
| `validate_vs_R.R`, randomized inputs against `physical_merge()` | 212 | 212 match |
| `validate_pipeline_R.R`, full `read_sumstat` to `export_snp_list` chain | 7 | 7 match |
| real chr22 HbA1c file, window 25 kb, `reset_on="best"` | 542 blocks | 542/542 identical, including `rps_ID` |

The randomized suite covers every combination of `reward` ∈ {min, max},
`reset_on` ∈ {best, any}, window ∈ {1 kb, 50 kb, 500 kb}, one and four
chromosomes, n ∈ {1, 5, 200, 5 000}, duplicated positions, and both sorted and
shuffled input (the latter through `--sort`). `start`, `end` and `rps_BP` are
compared with zero tolerance and `rps_value` as its 17-digit text form, so the
match is bit-exact rather than approximate.

## 3. Known differences from the R package

1. **Chromosome output order.** `read_sumstat()` sorts the whole table by
   `position` alone, so `physical_merge()` receives chromosomes interleaved and
   emits them in order of their smallest coordinate (e.g. 22, 2, 1 for a
   1/2/22 file). The C tool keeps input order. Per-chromosome block content is
   identical; only the `serial` numbering differs.
2. **`CHROM` column.** R omits it when the input holds a single chromosome; the
   C tool always emits it when a chromosome column is present, so the output
   schema does not depend on the data.
3. **Lead SNP at multi-allelic sites, found by this port and since fixed in R.**
   On the real chr22 file 2 of 542 blocks originally got a different `rps_ID`.
   Both are positions carrying two variants: 22:21092217 holds
   `22:21092217_CT_C` (P = 0.79) and `rs546487622` (P = 1.5 × 10⁻²⁵), and
   22:30939062 holds `rs142544112` (P = 0.15) and `rs545363690`
   (P = 2.3 × 10⁻⁹). `annotate_blocks()` deduplicated the input on position and
   kept the first row, so it labeled the block with the non-significant
   variant. `physical_merge()` now records the input row index of every
   representative SNP in `attr(blocks, "rps_row")` and `annotate_blocks()` uses
   it, so R and C agree on all 542 blocks including `rps_ID`. Blocks that do
   not carry the attribute still go through the old position join, so existing
   code is unaffected.

## 4. Reproducing

```
cd cli && make
Rscript validate_vs_R.R
Rscript validate_pipeline_R.R
cd bench && ./bench.sh 1000000 5000000 10000000
```
