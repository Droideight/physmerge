# physmerge — C command-line executable

A single-file C99 port of the physmerge R package, for running physical locus
merging on full-size summary statistics without loading them into R.

## Build

The source is one portable C99 file. It needs a C compiler and, only if you
want to read `.gz` input directly, zlib.

**macOS / Linux**

```
make                       # cc -O2 -std=c99 physmerge.c -lz
make nozlib                # no gzip support, zero dependencies
make install PREFIX=~/.local
```

macOS ships zlib and a compiler with the Xcode command line tools
(`xcode-select --install`). On Debian/Ubuntu: `sudo apt install build-essential zlib1g-dev`.

**Windows**

```
build.bat
```

`build.bat` uses MSVC (`cl`) if you are in a *x64 Native Tools Command Prompt
for VS*, otherwise MinGW-w64 `gcc`. Either way it produces `physmerge.exe` and
builds without zlib, so `.gz` input is not supported there; decompress first, or
add zlib and drop `/DPHYSMERGE_NO_ZLIB`. Manually:

```
cl /O2 /DPHYSMERGE_NO_ZLIB /Fephysmerge.exe physmerge.c        REM MSVC
gcc -O2 -std=c99 -DPHYSMERGE_NO_ZLIB -o physmerge.exe physmerge.c   REM MinGW-w64
```

There is one binary per platform, built from the same source; nothing else
differs. Windows line endings (CRLF) in the input are handled. Under WSL or
Git Bash, use the macOS/Linux instructions instead.

## Example

`example/demo.sh` runs the whole tool on `example/demo.glm.linear`, a 59-row
PLINK2 file that contains two chromosomes, covariate `TEST` rows, a missing
p-value and a multi-allelic site. On Windows run the commands inside it by
hand (`..\physmerge.exe --input demo.glm.linear --format plink2`).

## Usage

```
physmerge --input gwas.glm.linear --format plink2 \
          --sig-th 5e-8 --window 500000 --reset-on any \
          --out blocks.tsv --snp-list lead_snps.txt
```

`--help` lists every flag. The flags map one-to-one onto the R arguments:

| R | C |
|---|---|
| `read_sumstat(path, format=)` | `--input`, `--format plink2\|gpcm\|custom` |
| `chrom_col`, `pos_col`, `id_col`, `value_col` | `--chrom-col`, `--pos-col`, `--id-col`, `--value-col` |
| `test_filter`, `test_col`, `test_val` | `--test-filter` / `--no-test-filter`, `--test-col`, `--test-val` |
| `chrom = c(1,2)` | `--chrom 1,2` |
| `physical_merge(sig_th=, window=)` | `--sig-th`, `--window` |
| `reward = "min"/"max"` | `--reward min\|max` |
| `reset_on = "best"/"any"` | `--reset-on best\|any` |
| `annotate_blocks()` | on by default (`rps_ID`); `--annotate-full` for every original column |
| `export_snp_list()` | `--snp-list FILE`, `--snp-list-dir DIR` |

Input may be plain text, gzip, or `-` for stdin. The separator is auto-detected
from the header (tab, comma, whitespace) and can be forced with `--sep`.

## Memory model

The default path is a single streaming pass: only the current line and one
pending block are held, so resident memory is constant (~2.4 MB) regardless of
file size. This requires the input to be position-sorted within each
chromosome and each chromosome to occupy one contiguous block of lines, which
is true of PLINK2, REGENIE, BOLT-LMM and GCTA output. If that does not hold,
the tool stops with a message and `--sort` buffers and stable-sorts the records
exactly as R's `order()` would, at the cost of ~32 bytes per SNP.

## Fidelity to the R implementation

`validate_vs_R.R` runs 212 randomised comparisons against
`physmerge::physical_merge()` (all combinations of `reward`, `reset_on`,
window size, single/multi chromosome, duplicate positions, sorted and shuffled
input) and requires bit-identical `start`, `end`, `rps_BP` and `rps_value`.
`validate_pipeline_R.R` checks the whole `read_sumstat` → `physical_merge` →
`annotate_blocks` → `export_snp_list` chain on PLINK2-shaped input, including
the TEST filter, NA p-values, chromosome subsetting and gzip.

Both suites pass with zero mismatches. Two intentional differences:

1. **Chromosome output order.** `read_sumstat()` sorts the whole table by
   `position` alone, so `physical_merge()` emits chromosomes in order of their
   smallest coordinate (e.g. 22, 2, 1). The C tool keeps input order
   (1, 2, 22). Block contents are identical per chromosome; only `serial`
   numbering differs.
2. **`CHROM` column.** R drops it when the input holds a single chromosome;
   the C tool always writes it when a chromosome column is present.
3. **Lead SNP at multi-allelic sites.** This port exposed a bug in
   `annotate_blocks()`, which deduplicated the input on position and could
   therefore label a block with the non-significant variant at a multi-allelic
   site. Fixed in physmerge 0.3.0 via `attr(blocks, "rps_row")`; R and C now
   agree on every block. See `PERFORMANCE.md` §3.

## Benchmark

See `bench/bench.sh`. Results are in `PERFORMANCE.md`.
