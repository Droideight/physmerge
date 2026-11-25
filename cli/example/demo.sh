#!/bin/bash
set -e
cd "$(dirname "$0")"
PM=../physmerge

echo "### 1. defaults (window=500kb, threshold and direction from the P column, PLINK2 TEST filter)"
$PM --input demo.glm.linear

echo
echo "### 2. a stricter threshold keeps only the strongest peaks"
$PM --input demo.glm.linear --sig-th 1e-15 --quiet

echo
echo "### 3. a 50 kb window splits the chr1 peaks apart"
$PM --input demo.glm.linear --window 50000 --quiet

echo
echo "### 4. write the block table and the lead-SNP list to files"
$PM --input demo.glm.linear --out blocks.tsv --snp-list lead_snps.txt --quiet
cat lead_snps.txt

echo
echo "### 5. gzip input, chromosome 2 only, other original columns appended"
$PM --input demo.glm.linear.gz --chrom 2 --annotate-full --quiet

echo
echo "### 6. same thing in R, for comparison"
Rscript - <<'RS'
suppressMessages(library(physmerge))
d <- read_sumstat("demo.glm.linear")
print(annotate_blocks(physical_merge(d), d))
RS
