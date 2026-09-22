# Edge-case probes

Exploratory scripts, not a regression suite: they print what each function does
rather than asserting. `FINDINGS.txt` records what the 2026-09-22 run turned up.

| script | covers |
|---|---|
| `probe_read_sumstat.R` | `read_sumstat()`: formats, TEST and chromosome filters, NA, gzip, odd positions, p = 0/1/-1, duplicate and colliding column names, empty files |
| `probe_physical_merge.R` | `physical_merge()`: 30 cases -- empty, single row, all/none significant, NA, duplicate positions, unsorted input, window 1 and 3e8, negative positions, `reward=max`, both `reset_on` rules, multi-chromosome including factors, strings, and NA |
| `probe_export.R` | `annotate_blocks()` and `export_snp_list()`: multi-allelic leads, the fallback join, missing ID columns, `keep_*` toggles, zip export, bad paths |
| `probe_cli.sh` | the C tool: quoted CSV, BOM, blank lines, underflow, fractional and negative positions, `--no-chrom`, `--chrom`, `--annotate-full`, name collisions, unwritable outputs |
| `diff_c_vs_r.R` | 400 randomised inputs through both the C tool and R, comparing block tables |
| `gen_sas_vectors.R` | regenerates `sas/testdata/` from the R implementation |

Run from the package root, after `make -C cli`:

```bash
Rscript tests/qa/probe_read_sumstat.R
Rscript tests/qa/probe_physical_merge.R
Rscript tests/qa/probe_export.R
Rscript tests/qa/diff_c_vs_r.R
bash    tests/qa/probe_cli.sh
```
