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
implementations that return the same blocks: an R package, a standalone
executable, and a Base SAS port. The executable reads the file in one streaming
pass, so its memory use does not grow with the input; a 1.85 GB file ran in
2.4 MB of memory.

## How it works

`physical_merge()` makes a single forward pass over position-sorted summary
statistics. A block opens at the first significant SNP, with its start placed one
window upstream (`max(0, position - window)`). The block carries a window-sized
budget. The distance travelled spends the budget, and every significant SNP
refills it to the full window. The block stays open until the budget runs out,
which happens when the next significant SNP lies one window or more beyond the
previous one. The block then closes one window downstream of its last
significant SNP, mirroring its start, and the next significant SNP opens a new
block. The pass ends at the last position.

This refill-at-every-significant-SNP rule is `reset_on = "any"`, the default.
With `reset_on = "best"`, only a SNP more significant than the block's current
representative refills the budget, so a block can close while significant SNPs
continue.

The representative of each block is its most significant SNP. Each block is
padded by one window on both sides, so two adjacent blocks overlap whenever the
gap between the last significant SNP of one and the first significant SNP of the
next is less than two windows. A final trim step moves the end of the earlier
block back to the start of the next one. The result is contiguous, strictly
non-overlapping blocks. A block is two windows wide only when it is isolated; a
trimmed block is narrower. With a chromosome column the algorithm runs per
chromosome.

---

## 1. Install

### R package

```r
install.packages("devtools")
devtools::install_github("Droideight/physmerge")
```

The R package needs `data.table`, which `install_github()` pulls in; everything
else is base R and `utils`.

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
reach that range, merge on `LOG10_P` with `reward=max`.

---

## 2. Quick start

The repository ships a small example file. From the `cli` directory:

```bash
cd example
physmerge --input demo.glm.linear --format plink2
```

From the top of the repository it is `cd cli/example` instead.

```
physmerge: TEST filter: kept 52 of 59 rows where TEST = 'ADD'.
physmerge: 1 row(s) dropped (NA in position or value).
physmerge: 51 SNPs -> 3 blocks (window=500000, sig_th=5e-08, reward=min, reset_on=any).
serial  CHROM  start    end      rps_BP   rps_ID    rps_P
1       1      700000   1855000  1265000  rs100007  2.2e-14
2       1      2300000  3300000  2800000  rs100029  1.4e-15
3       2      0        1020000  500000   rs100042  9.9e-20
```

`./demo.sh` in that directory runs six variations on the same file.

The equivalent in R:

```r
library(physmerge)
d <- read_sumstat("demo.glm.linear", format = "plink2")
b <- physical_merge(d$data, sig_th = 5e-8, window = 500000,
                    reward = d$reward, chrom_col = "CHROM")
b <- annotate_blocks(b, d$data)
export_snp_list(b, "lead_snps.txt")
```

And in SAS:

```sas
%pm_read(path=demo.glm.linear, out=ss, format=plink2);
%physmerge(data=ss, out=blocks, sig_th=5e-8, window=500000,
           reward=&pm_reward,
           chrom=&pm_chrom_var, pos=POSITION, value=VALUE, id=&pm_id_var);
%pm_export(blocks=blocks, file=lead_snps.txt);
```

---

## 3. Output destination

By default the block table is printed to the terminal (stdout) and nothing is
written to disk. Progress messages go to stderr, so a redirect or a pipe carries
the table only.

To write files, name them yourself:

| Flag | What it writes |
|---|---|
| `--out FILE` | the block table (tab-separated, with a header) |
| `--snp-list FILE` | one representative SNP id per line |
| `--snp-list-dir DIR` | `snp_ch1.txt`, `snp_ch2.txt`, … one per chromosome |

The files are written where the flag points, so give a full path to keep them out
of the current directory:

```bash
mkdir -p ~/physmerge_out
physmerge --input gwas.glm.linear --format plink2 \
  --out ~/physmerge_out/blocks.tsv \
  --snp-list ~/physmerge_out/lead_snps.txt
```

To pipe instead of writing a file:

```bash
physmerge --input gwas.glm.linear --format plink2 --quiet | head
```

### Output columns

In R, `read_sumstat()` reads and filters the file. `physical_merge()` returns
one row per block: serial, chromosome, start, end, representative position and
value. `annotate_blocks()` joins the original fields back to each
representative, and `export_snp_list()` writes the representative SNP IDs as a
flat file or a per-chromosome ZIP. The executable does all four steps in one
call.

| Column | Meaning |
|---|---|
| `serial` | block number, running across chromosomes |
| `CHROM` | chromosome (present when the input has a chromosome column) |
| `start` | block start in bp, one window upstream of the first significant SNP |
| `end` | block end in bp; guaranteed not to reach into the next block |
| `rps_BP` | position of the representative (most significant) SNP |
| `rps_ID` | its id, when the input has an id column |
| `rps_<VALUE>` | its p-value or statistic, named after the input column |

Blocks do not overlap; within a chromosome, `end[i] <= start[i+1]`.

---

## 4. Recipes

For standard PLINK2 `.glm.*` output, use the `plink2` format. It keeps only
`TEST=ADD` rows and drops rows with a missing p-value, so no pre-filtering is
needed:

```bash
physmerge --input gwas.glm.linear --format plink2 \
  --sig-th 5e-8 --window 500000 \
  --out blocks.tsv --snp-list lead_snps.txt
```

If the value column is `-log10(P)` rather than `P` (PLINK2 writes
`NEG_LOG10_P` for some runs), point at the column, flip the direction, and
convert the threshold (`-log10(5e-8) = 7.30103`):

```bash
physmerge --input gwas.glm.logistic.hybrid --format plink2 \
  --value-col NEG_LOG10_P --reward max --sig-th 7.30103 \
  --out blocks.tsv
```

Leave the TEST filter on. `--no-test-filter` admits the DOMDEV and RECESSIVE
rows as well, which can produce a block whose representative is not an additive
test; use it only for a file that has no TEST column.

For any other table, use `--format custom` and name the columns. The file can
be space-, tab- or comma-separated; the separator is read from the header:

```bash
physmerge --input sumstats.txt --format custom \
  --chrom-col CHR --pos-col POS --id-col SNP --value-col P \
  --sig-th 5e-8 --window 500000
```

Feed the lead SNPs straight into PLINK:

```bash
plink2 --pfile your_data --extract lead_snps.txt --make-pgen --out lead_only
```

A larger window merges more. Where significant SNPs are dense and never more
than one window apart, `--window 500000` chains the whole region into a single
block; shrink the window for finer loci. In a chr22 HbA1c scan (1.25 million
SNPs, 2,550 of them genome-wide significant), 500 kb returned 1 block and 25 kb
returned 542.

---

## 5. Options

| R argument | Command-line flag |
|---|---|
| `read_sumstat(path, format=)` | `--input`, `--format plink2\|gpcm\|custom` |
| `chrom_col`, `pos_col`, `id_col`, `value_col` | `--chrom-col`, `--pos-col`, `--id-col`, `--value-col` |
| `test_filter`, `test_col`, `test_val` | `--test-filter` / `--no-test-filter`, `--test-col`, `--test-val` |
| `chrom = c(1, 2)` | `--chrom 1,2` |
| `sig_th` | `--sig-th 5e-8` |
| `window` | `--window 500000` |
| `reward = "min"` / `"max"` | `--reward min\|max` |
| `reset_on = "any"` (default) / `"best"` | `--reset-on any\|best` |
| `annotate_blocks()` | the executable always adds the lead SNP id; `--annotate-full` appends every original column |
| `export_snp_list()` | `--snp-list FILE`, `--snp-list-dir DIR` |

Other flags: `--sep` to force a separator, `--sort` for input that is not
position-sorted, `--no-header`, `--quiet`, `--help`.

`--reset-on` controls how a block stays open. With `any`, the default, every
significant SNP refills the window, which gives the union of the ±window
intervals around all significant SNPs. With `best`, only a more significant SNP
refills it. A long run of comparably significant SNPs is then cut roughly every
window into contiguous pieces. Use `any` when such a region should come out as
one block.

To merge a file as one sequence, pass `--no-chrom` on the command line, or
`chrom_col = NA` to `read_sumstat()`. Both work whether or not the file has a
chromosome column, and both ignore it when it is there.

Input may be plain text, gzip (`.gz`), or `-` for stdin. A file that is not
position-sorted within a chromosome is rejected by the executable with a message
instead of being merged wrongly; add `--sort` in that case. `read_sumstat()`
sorts on read, so the R path needs nothing.

Block ends are positions in base pairs, not clipped to any assembly, so the last
block of a chromosome can end past the chromosome's length by up to one window.

---

## 6. More

- `physmerge --help`: every flag

MIT licensed.
