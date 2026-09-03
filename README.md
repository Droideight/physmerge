# physmerge

Linkage disequilibrium (LD) within a GWAS study can produce spurious hits; this is
commonly addressed using an LD reference panel to perform clumping. In situations
where only summary statistics are available, however, distance-based locus
definition is the panel-free alternative but is typically applied as ad-hoc
per-study code.

physmerge collapses significant SNPs into non-overlapping locus blocks from
summary statistics alone using a forward sliding-window rule (a key parameter
`reset_on = "any"`: an open block is extended whenever the next significant SNP
lies within the window of the current one), yielding contiguous, strictly
non-overlapping locus blocks directly from summary statistics. This program is
available both as an R package and an executable, in this repository. The two
return the same blocks; the executable reads the file in one streaming pass, so
its memory use stays at about 2.4 MB (1.85 GB input, 2.4 MB resident).

## How it works

`physical_merge()` makes a single forward pass over position-sorted summary
statistics. A block opens at the first significant SNP, with its start placed one
window upstream (`max(0, position - window)`). The block carries a window-sized
budget that is spent by the distance traveled and refilled to the full window at
every significant SNP (`reset_on = "any"`); it stays open until the budget runs
out, equivalently, until the next significant SNP lies one window or more beyond
the previous one, at which point it closes one window downstream of the last
significant SNP, mirroring its start. A new block opens upon the next significant
SNP until the last position is visited.

The representative of each block is its most significant SNP. By construction the
representatives of successive blocks are at least one window apart, but because
each block is padded by one window on both sides, adjacent blocks still overlap
whenever that gap is less than two windows; a final trim step therefore shortens
any block whose downstream-extended end runs past the next block's
upstream-extended start, giving contiguous, strictly non-overlapping blocks. With
a chromosome column the algorithm runs per chromosome.

---

## 1. Install

### R package

```r
install.packages("devtools")           # if you do not have it
devtools::install_github("Droideight/physmerge")
```

The R package depends only on base R and `utils` (`data.table` optional for fast
reading).

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

---

## 2. Quick start

The repository ships a small example file, so the commands below run without any
data of your own. From the `cli` directory you were left in by the build:

```bash
cd example
physmerge --input demo.glm.linear --format plink2
```

From the top of the repository it is `cd cli/example` instead.

```
physmerge: TEST filter: kept 52 of 59 rows where TEST = 'ADD'.
physmerge: 1 row(s) dropped (NA in position or value).
physmerge: 51 SNPs -> 3 blocks (window=500000, sig_th=5e-08, reward=min, reset_on=best)
serial  CHROM  start    end      rps_BP   rps_ID    rps_P
1       1      700000   1765000  1265000  rs100007  2.2e-14
2       1      2300000  3300000  2800000  rs100029  1.4e-15
3       2      0        1000000  500000   rs100042  9.9e-20
```

`./demo.sh` in that directory walks through six variations of the same file.

The equivalent in R:

```r
library(physmerge)
d <- read_sumstat("demo.glm.linear", format = "plink2")
b <- physical_merge(d$data, sig_th = 5e-8, window = 500000,
                    reward = d$reward, reset_on = "any", chrom_col = "CHROM")
b <- annotate_blocks(b, d$data)
export_snp_list(b, "lead_snps.txt")
```

---

## 3. Where the output goes

**By default the block table is printed to the terminal (stdout)** and nothing is
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

`physical_merge()` returns one row per block (serial, chromosome, start, end,
representative position and value); `annotate_blocks()` joins original fields back
to each representative, and `export_snp_list()` writes representative SNP IDs as a
flat file or per-chromosome ZIP. The executable does all four steps in one call.


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

**Standard PLINK2 `.glm.*` output.** The `plink2` format keeps only `TEST=ADD`
rows and drops rows with a missing p-value, so no pre-filtering is needed:

```bash
physmerge --input gwas.glm.linear --format plink2 \
  --sig-th 5e-8 --window 500000 --reset-on any \
  --out blocks.tsv --snp-list lead_snps.txt
```

**The value column is `-log10(P)` rather than `P`** (PLINK2 writes
`NEG_LOG10_P` for some runs). Point at the column, flip the direction, and
convert the threshold (`-log10(5e-8) = 7.30103`):

```bash
physmerge --input gwas.glm.logistic.hybrid --format plink2 \
  --value-col NEG_LOG10_P --reward max --sig-th 7.30103 --no-test-filter \
  --out blocks.tsv
```

**Any other table**, space-, tab- or comma-separated; the separator is read from
the header. Name the columns:

```bash
physmerge --input sumstats.txt --format custom \
  --chrom-col CHR --pos-col POS --id-col SNP --value-col P \
  --sig-th 5e-8 --window 500000
```

**Feed the lead SNPs straight into PLINK:**

```bash
plink2 --pfile your_data --extract lead_snps.txt --make-pgen --out lead_only
```

**Choosing a window.** A larger window merges more. Where significant SNPs are
dense and never more than one window apart, `--window 500000 --reset-on any`
chains the whole region into a single block; shrink the window for finer loci. In
a chr22 HbA1c scan (1.25 million SNPs, 2,550 of them genome-wide significant),
500 kb returned 1 block; in comparison, 25 kb returned 542.

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
| `reset_on = "best"` / `"any"` | `--reset-on best\|any` |
| `annotate_blocks()` | on by default; `--annotate-full` appends every original column |
| `export_snp_list()` | `--snp-list FILE`, `--snp-list-dir DIR` |

Other flags: `--sep` to force a separator, `--sort` for input that is not
position-sorted, `--no-header`, `--quiet`, `--help`.

`--reset-on` controls how a block stays open. `any` extends an open block whenever
the next significant SNP lies within the window of the current one, which is the
union of the ±window intervals around all significant SNPs; `best` refills the
window only when a more significant SNP appears.

Input may be plain text, gzip (`.gz`), or `-` for stdin. A file that is not
position-sorted within a chromosome is rejected with a message instead of being
merged wrongly; add `--sort` in that case.

---

## 6. More

- `cli/README.md`: build details, memory model, differences from the R package
- `cli/PERFORMANCE.md`: benchmarks and the validation suites
- `TECHNICAL_SPEC.md`: algorithm, data contract, edge cases (in Chinese)
- `physmerge --help`: every flag

MIT licensed.
