#!/usr/bin/env Rscript
## Generate PLINK2-shaped benchmark files: n rows spread over 22 chromosomes,
## ~0.4% of SNPs genome-wide significant, sorted by (#CHROM, POS).
suppressMessages(library(data.table))
args <- commandArgs(TRUE); n <- as.numeric(args[1]); out <- args[2]
set.seed(1)
chr <- sort(sample.int(22L, n, replace = TRUE))
pos <- integer(n)
for (c in 1:22) { i <- which(chr == c); pos[i] <- sort(sample.int(2.5e8, length(i))) }
p <- runif(n)
hit <- sample.int(n, max(1L, round(n * 0.004)))
p[hit] <- 10^(-runif(length(hit), 8, 30))
dt <- data.table(`#CHROM` = chr, POS = pos, ID = paste0("rs", seq_len(n)),
                 REF = "A", ALT = "G", A1 = "G", TEST = "ADD",
                 OBS_CT = 350000L, BETA = round(rnorm(n), 5), SE = 0.01,
                 T_STAT = round(rnorm(n), 4), P = p)
fwrite(dt, out, sep = "\t", quote = FALSE)
cat(sprintf("wrote %s (%.0f rows)\n", out, n))
