# physmerge

Linkage disequilibrium (LD) makes neighbouring SNPs in a GWAS significant
together, so one association signal shows up as many correlated hits. This is
commonly addressed by clumping against an LD reference panel. When only summary
statistics are available, distance-based locus definition is the panel-free
alternative, but it is usually written as ad-hoc code for each study.

physmerge collapses significant SNPs into contiguous, strictly non-overlapping
locus blocks, using summary statistics alone. It applies a forward
sliding-window rule: an open block is extended whenever the next significant SNP
lies within one window of the previous one. This repository carries three
implementations that take the same options and return the same blocks: an R
package, a standalone executable, and a Base SAS port. The executable reads the
file in one streaming pass, so its memory use does not grow with the input; a
1.85 GB file ran in 2.4 MB of memory.

## How it works

physmerge makes a single forward pass over position-sorted summary statistics.
A block opens at the first significant SNP, with its start placed one window
upstream (`max(0, position - window)`). The block carries a window-sized
budget. The distance travelled spends the budget, and every significant SNP
refills it to the full window. The block stays open until the budget runs out,
which happens when the next significant SNP lies one window or more beyond the
previous one. The block then closes one window downstream of its last
significant SNP, mirroring its start, and the next significant SNP opens a new
block. The result is the union of the ±window intervals around all significant
SNPs.

The representative of each block is its most significant SNP. Each block is
padded by one window on both sides, so two adjacent blocks overlap whenever the
gap between the last significant SNP of one and the first significant SNP of the
next is less than two windows. A final trim step moves the end of the earlier
block back to the start of the next one. The result is contiguous, strictly
non-overlapping blocks. A block is two windows wide only when it is isolated; a
trimmed block is narrower. With a chromosome column the merge runs per
chromosome.

---

## 1. Install

### R package

```r
install.packages("devtools")
devtools::install_github("Droideight/physmerge")
```

The R package needs `data.table` and `R.utils`, which `install_github()` pulls
in; everything else is base R and `utils`.

### C executable, macOS

```bash
git clone https://github.com/Droideight/physmerge.git
cd physmerge/cli
make
make install PREFIX=~/.local
```

This needs the Xcode command line tools (`xcode-select --install`). If the shell
cannot find `physmerge` afterwards, put `~/.local/bin` on the PATH:

```bash
echo 'export PATH="$HOME/.local/bin:$PATH"' >> ~/.zshrc && source ~/.zshrc
```

### C executable, Linux

Same as macOS. On Debian/Ubuntu install the toolchain first:

```bash
sudo apt install build-essential zlib1g-dev
```

### C executable, Windows

```
git clone https://github.com/Droideight/physmerge.git
cd physmerge\cli
build.bat
```

`build.bat` uses MSVC (`cl`) inside an *x64 Native Tools Command Prompt for VS*,
and MinGW-w64 `gcc` otherwise. Either produces `physmerge.exe`. The Windows build
carries no zlib, so `.gz` input has to be decompressed first; plain text, which is
what PLINK2 writes, is unaffected. Under WSL or Git Bash, follow the Linux
instructions instead.

Check it works:

```bash
physmerge --version
```

### SAS

Point `PMDIR` at the `sas` directory of a clone and include the file:

```sas
%let PMDIR = /path/to/physmerge/sas;
%include "&PMDIR/physmerge.sas";
%pm_version;
```

Two SAS-specific traps: a missing value compares below every number, and
`BEST32.` does not read p-values far below 1e-300 reliably. If your p-values
reach that range, merge on the `NEG_LOG10_P` column instead. Reading `.gz`
files needs SAS 9.4M5 or later.

---

## 2. Quick start

The repository ships a small example file. From the `cli` directory:

```bash
cd example
physmerge --input demo.glm.linear
```

From the top of the repository it is `cd cli/example` instead.

```
physmerge: TEST filter: kept 52 of 59 rows where TEST = 'ADD'.
physmerge: 1 row(s) dropped (NA in position or value).
physmerge: 51 SNPs -> 3 blocks (window=500000, sig_th=5e-08, reward=min).
serial  CHROM  start    end      rps_BP   rps_ID    rps_P
1       1      700000   1855000  1265000  rs100007  2.2e-14
2       1      2300000  3300000  2800000  rs100029  1.4e-15
3       2      0        1020000  500000   rs100042  9.9e-20
```

`./demo.sh` in that directory runs six variations on the same file.

The same in R:

```r
library(physmerge)
d <- read_sumstat("demo.glm.linear")
b <- physical_merge(d)
b <- annotate_blocks(b, d)
export_snp_list(b, file = "lead_snps.txt")
```

And in SAS:

```sas
%pm_read(path=demo.glm.linear, out=ss);
%physmerge(data=ss, out=blocks);
%pm_annotate(blocks=blocks, data=ss, out=blocks_annot);
%pm_export(blocks=blocks, file=lead_snps.txt);
```

All three print the same progress messages and return the same three blocks.

---

## 3. Output

The executable prints the block table to the terminal (stdout) and writes
nothing to disk unless asked. Progress messages go to stderr, so a redirect or a
pipe carries the table only. R returns a data frame and SAS a data set.

| What | R | Command line | SAS |
|---|---|---|---|
| block table | the return value | `--out FILE`, or stdout | `out=` |
| other original columns of each lead SNP | `annotate_blocks(b, d)` | `--annotate-full` | `%pm_annotate` |
| lead SNP ids, one per line | `export_snp_list(b, file = )` | `--snp-list FILE` | `%pm_export(file=)` |
| one id file per chromosome | `export_snp_list(b, dir = )` | `--snp-list-dir DIR` | `%pm_export(dir=)` |

A per-chromosome directory holds `snp_ch1.txt`, `snp_ch2.txt`, … for the
chromosomes that have blocks, and is created if it does not exist. The id is
`rps_ID` when the input has an id column, otherwise `rps_BP`.

### Block table columns

| Column | Meaning |
|---|---|
| `serial` | block number, running across chromosomes |
| `CHROM` | chromosome (present when the input has a chromosome column) |
| `start` | block start in bp, one window upstream of the first significant SNP |
| `end` | block end in bp; never reaches into the next block |
| `rps_BP` | position of the representative (most significant) SNP |
| `rps_ID` | its id, when the input has an id column |
| `rps_<VALUE>` | its value, named after the input column, e.g. `rps_P` or `rps_NEG_LOG10_P` |

Blocks do not overlap; within a chromosome, `end[i] <= start[i+1]`. Adding the
original columns appends every other column of the lead SNP's row, in file
order; the chromosome, position, id and value columns are not repeated.

Block ends are positions in base pairs, not clipped to any assembly, so the last
block of a chromosome can end past the chromosome's length by up to one window.

---

## 4. Recipes

For standard PLINK2 `.glm.*` output, nothing needs to be set. The `plink2`
format reads `#CHROM`, `POS`, `ID` and `P`, keeps only `TEST = ADD` rows and
drops rows with a missing p-value:

```bash
physmerge --input gwas.glm.linear --out blocks.tsv --snp-list lead_snps.txt
```

If PLINK2 was run with `--glm log10`, the p-values are in `NEG_LOG10_P`. Point
at that column; the threshold and direction follow from its name:

```bash
physmerge --input gwas.glm.logistic.hybrid --value-col NEG_LOG10_P
```

For a dominant or recessive model, PLINK2 writes `DOM` or `REC` in the `TEST`
column instead of `ADD`:

```bash
physmerge --input gwas.glm.linear --test-val DOM
```

For any other table, use `--format custom` and name the columns:

```bash
physmerge --input sumstats.txt --format custom \
  --chrom-col CHR --pos-col BP --id-col SNP --value-col P
```

To merge on a test statistic or any column whose name physmerge does not
recognise, set the direction and the threshold yourself:

```bash
physmerge --input sumstats.txt --format custom \
  --chrom-col CHR --pos-col BP --value-col CHISQ --reward max --sig-th 29.7
```

Genotype data split by chromosome, as in UK Biobank, takes one id list per
chromosome:

```bash
physmerge --input gwas.glm.linear --snp-list-dir lead_by_chrom
plink2 --bgen ukb_imp_chr1_v3.bgen ref-first --sample ukb.sample \
  --extract lead_by_chrom/snp_ch1.txt --make-pgen --out lead_chr1
```

A larger window merges more. Where significant SNPs are dense and never more
than one window apart, the default 500 kb window chains the whole region into a
single block; shrink the window for finer loci. In a chr22 HbA1c scan (1.25
million SNPs, 2,550 of them genome-wide significant), 500 kb returned 1 block
and 25 kb returned 542.

---

## 5. Options

| R | Command line | SAS | Default |
|---|---|---|---|
| `read_sumstat(path)` | `--input FILE` | `%pm_read(path=)` | |
| `format` | `--format` | `format=` | `plink2`; or `custom` |
| `chrom_col` | `--chrom-col` | `chrom_col=` | plink2: `#CHROM`; `NA` ignores it |
| `pos_col` | `--pos-col` | `pos_col=` | plink2: `POS` |
| `id_col` | `--id-col` | `id_col=` | plink2: `ID`; `NA` for none |
| `value_col` | `--value-col` | `value_col=` | plink2: `P` |
| `test_val` | `--test-val` | `test_val=` | `ADD` (plink2 only) |
| `chrom = c(1, 2)` | `--chrom 1,2` | `chrom=1 2` | all chromosomes |
| `sep` | `--sep` | `sep=` | read from the header |
| `physical_merge(sig_th = )` | `--sig-th` | `%physmerge(sig_th=)` | from the value column |
| `window` | `--window` | `window=` | `500000` |
| `reward` | `--reward min\|max` | `reward=` | from the value column |
| `sort` | `--sort` | `sort=1` | off |
| `quiet` | `--quiet` | `quiet=1` | off |

Threshold and direction. When they are not given, the name of the value column
decides. `P`, `PVAL`, `PVALUE`, `P_VALUE`, `P-VALUE` and any name starting with
`P_` give `reward = min` and `sig_th = 5e-8`. `NEG_LOG10_P` and `LOG10_P` give
`reward = max` and `sig_th = 7.30103`, which is -log10(5e-8). Any other column
stops with a message asking for both. A T or Z statistic is not guessed, because
its sign would hide strong negative effects.

Separator. `sep` takes `tab`, `space` or a single character. When it is not
given, the header decides: tab if it has a tab, otherwise comma, otherwise
semicolon, otherwise whitespace.

Chromosome. With `chrom_col` set to `NA`, the chromosome column is ignored and
the whole file is merged as one sequence.

Sorting. Input must be sorted by position within each chromosome, with each
chromosome in one stretch, as PLINK2 writes it. Otherwise all three stop with a
message such as `input is not position-sorted (chromosome 1: 2988324 after
2997970)`. With `sort` on, the input is sorted first; the executable then holds
the whole file in memory instead of streaming it.

TEST filter. If no row has the requested `TEST` value, all three stop and list
the values the file does have, for example `TEST values in this file: DOM, sex`.

Input may be plain text or gzip (`.gz`); the executable also reads `-` as stdin.
`quiet` hides the progress lines; warnings and errors are always shown.

---

## 6. More

- `physmerge --help`: every flag

MIT licensed.
