#!/bin/bash
# Walk-through of the physmerge executable on example/demo.glm.linear
# (59 rows, PLINK2 shape: 2 chromosomes, covariate TEST rows, one NA p-value,
#  and one multi-allelic site at 1:1265000).
set -e
cd "$(dirname "$0")"
PM=../physmerge

echo "### 1. defaults (reset_on=best, window=500kb, PLINK2 TEST filter)"
$PM --input demo.glm.linear --format plink2

echo
echo "### 2. reset_on=any -- the two chr1 peaks 155 kb apart now merge"
$PM --input demo.glm.linear --format plink2 --reset-on any --quiet

echo
echo "### 3. a 50 kb window splits them again"
$PM --input demo.glm.linear --format plink2 --reset-on any --window 50000 --quiet

echo
echo "### 4. write the block table and the lead-SNP list to files"
$PM --input demo.glm.linear --format plink2 --reset-on any \
    --out blocks.tsv --snp-list lead_snps.txt --quiet
cat lead_snps.txt

echo
echo "### 5. gzip input, chromosome 2 only, full original columns kept"
$PM --input demo.glm.linear.gz --format plink2 --chrom 2 --annotate-full --quiet

echo
echo "### 6. same thing in R, for comparison"
Rscript - <<'RS'
suppressMessages(library(physmerge))
d <- read_sumstat("demo.glm.linear", format = "plink2")
b <- physical_merge(d$data, sig_th = 5e-8, window = 500000,
                    reward = d$reward, reset_on = "any", chrom_col = "CHROM")
print(annotate_blocks(b, d$data)[, c("serial","CHROM","start","end","rps_BP","rps_ID","rps_value")])
RS
