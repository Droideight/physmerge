#!/usr/bin/env Rscript
## R reference pipeline: read + merge + annotate + export ids
suppressMessages(library(physmerge))
a <- commandArgs(TRUE); f <- a[1]; ids <- a[2]
rs <- suppressMessages(read_sumstat(f, format = "plink2"))
b  <- physical_merge(rs$data, sig_th = 5e-8, window = 500000,
                     reward = rs$reward, reset_on = "any", chrom_col = "CHROM")
b  <- annotate_blocks(b, rs$data)
suppressMessages(export_snp_list(b, ids))
cat(sprintf("R blocks: %d\n", nrow(b)))
