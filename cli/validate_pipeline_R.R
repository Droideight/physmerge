#!/usr/bin/env Rscript
## End-to-end check: read_sumstat() + physical_merge() + annotate_blocks() in R
## versus the single C invocation, on PLINK2-shaped input (TEST filter, NA
## p-values, multiple chromosomes, gzip, chromosome subsetting).
suppressMessages({library(physmerge); library(data.table)})
BIN <- normalizePath("./physmerge", mustWork = TRUE)
TMP <- tempfile("pmpipe_"); dir.create(TMP)
set.seed(7)

n <- 40000
d <- data.table(`#CHROM` = sample(c(1, 2, 22), n, TRUE),
                POS = sample(1:5e6, n, TRUE),
                ID  = paste0("rs", seq_len(n)),
                TEST = "ADD",
                P = as.numeric(sprintf("%.17g", runif(n))))
hit <- sample(n, 1500); d$P[hit] <- as.numeric(sprintf("%.17g", 10^(-runif(1500, 8, 25))))
d$P[sample(n, 300)] <- NA                                  # NA p-values
extra <- copy(d[sample(n, 5000)]); extra$TEST <- "COV1"     # rows the TEST filter must drop
extra$P <- as.numeric(sprintf("%.17g", runif(5000)))
d <- rbind(d, extra)
setorderv(d, c("#CHROM", "POS"))
f <- file.path(TMP, "gwas.tsv")
dw <- copy(d); dw$P <- ifelse(is.na(dw$P), "NA", sprintf("%.17g", dw$P))
fwrite(dw, f, sep = "\t", quote = FALSE)
system(paste("gzip -kf", shQuote(f)))

check <- function(tag, chrom = NULL, reset = "any", win = 5e5, gz = FALSE) {
  rs <- suppressMessages(read_sumstat(if (gz) paste0(f, ".gz") else f,
                                      format = "plink2", chrom = chrom))
  rb <- physical_merge(rs$data, sig_th = 5e-8, window = win,
                       reward = rs$reward, reset_on = reset, chrom_col = "CHROM")
  rb <- annotate_blocks(rb, rs$data)
  args <- c("--input", shQuote(if (gz) paste0(f, ".gz") else f), "--format", "plink2",
            "--sig-th", "5e-8", "--window", format(win, scientific = FALSE),
            "--reset-on", reset, "--quiet")
  if (!is.null(chrom)) args <- c(args, "--chrom", paste(chrom, collapse = ","))
  out <- system(paste(shQuote(BIN), paste(args, collapse = " ")), intern = TRUE)
  cb <- read.delim(text = paste(out, collapse = "\n"), colClasses = "character")
  ok <- TRUE
  msg <- function(...) { cat(sprintf("  FAIL %s: %s\n", tag, paste0(...))); ok <<- FALSE }
  if (nrow(rb) != nrow(cb)) msg(sprintf("block count R=%d C=%d", nrow(rb), nrow(cb)))
  else if (nrow(rb) > 0) {
    ## read_sumstat() sorts globally by position, so R emits chromosomes in
    ## order of first appearance after that sort; the C tool keeps file order.
    ## Compare block content per chromosome instead of by serial.
    rb <- rb[order(as.character(rb$CHROM), rb$start), ]
    cb <- cb[order(cb$CHROM, as.numeric(cb$start)), ]
    if (!identical(as.character(rb$CHROM), cb$CHROM)) msg("CHROM")
    for (v in c("start", "end", "rps_BP"))
      if (!identical(as.numeric(rb[[v]]), as.numeric(cb[[v]]))) msg(v)
    if (!identical(as.character(rb$rps_ID), cb$rps_ID)) msg("rps_ID")
    ## R's sprintf/strtod are not always exact at extreme exponents; compare
    ## numerically at 1e-15 relative tolerance
    if (!isTRUE(all.equal(rb$rps_value, as.numeric(cb$rps_P), tolerance = 1e-15)))
      msg("rps_value")
  }
  cat(sprintf("%-38s R=%3d C=%3d %s\n", tag, nrow(rb), nrow(cb), if (ok) "OK" else "MISMATCH"))
  ok
}
res <- c(
  check("plain / reset=any / 500kb"),
  check("plain / reset=best / 500kb", reset = "best"),
  check("plain / reset=any / 50kb", win = 5e4),
  check("chrom subset 2,22", chrom = c(2, 22)),
  check("gzip input", gz = TRUE),
  check("gzip + chrom 1", chrom = 1, gz = TRUE)
)
## SNP-list export
rs <- suppressMessages(read_sumstat(f, format = "plink2"))
rb <- annotate_blocks(physical_merge(rs$data, 5e-8, 5e5, reward = "min",
                                     reset_on = "any", chrom_col = "CHROM"), rs$data)
p1 <- file.path(TMP, "r_ids.txt"); suppressMessages(export_snp_list(rb, p1))
p2 <- file.path(TMP, "c_ids.txt")
system(paste(shQuote(BIN), "--input", shQuote(f), "--format plink2 --sig-th 5e-8",
             "--window 500000 --reset-on any --quiet --out /dev/null --snp-list", shQuote(p2)))
same <- identical(sort(readLines(p1)), sort(readLines(p2)))
cat(sprintf("%-38s %s\n", "export_snp_list vs --snp-list", if (same) "OK" else "MISMATCH"))
res <- c(res, same)
cat(sprintf("\n%d/%d checks passed\n", sum(res), length(res)))
if (!all(res)) quit(status = 1)
