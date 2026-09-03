#!/usr/bin/env Rscript
## Fidelity harness: compare the C executable against the R reference
## implementation (physmerge::physical_merge) over randomised inputs.
suppressMessages(library(physmerge))
BIN <- normalizePath("./physmerge", mustWork = TRUE)
TMP <- tempfile("pmval_"); dir.create(TMP)
set.seed(20260903)

run_c <- function(df, sig, win, reward, reset, sorted = TRUE, multi = TRUE) {
  f <- file.path(TMP, "in.tsv")
  d2 <- df; d2$value <- sprintf("%.17g", df$value)   # full-precision round-trip
  d2$position <- sprintf("%.17g", df$position)
  write.table(d2, f, sep = "\t", quote = FALSE, row.names = FALSE)
  args <- c("--input", shQuote(f), "--format", "custom",
            "--pos-col", "position", "--value-col", "value", "--id-col", "NA",
            "--sig-th", format(sig, scientific = TRUE, digits = 17),
            "--window", format(win, digits = 17),
            "--reward", reward, "--reset-on", reset, "--quiet")
  args <- c(args, if (multi) c("--chrom-col", "CHROM") else "--no-chrom")
  if (!sorted) args <- c(args, "--sort")
  out <- system(paste(shQuote(BIN), paste(args, collapse = " ")), intern = TRUE)
  if (length(out) <= 1) return(data.frame())
  read.delim(text = paste(out, collapse = "\n"), colClasses = "character")
}

cmp <- function(rb, cb, tag) {
  if (nrow(rb) != nrow(cb))
    return(sprintf("%s: block count R=%d C=%d", tag, nrow(rb), nrow(cb)))
  if (nrow(rb) == 0L) return(NA_character_)
  rb <- rb[order(rb$serial), ]; cb <- cb[order(as.numeric(cb$serial)), ]
  for (v in c("start", "end", "rps_BP")) {
    if (!isTRUE(all.equal(as.numeric(rb[[v]]), as.numeric(cb[[v]]), tolerance = 0)))
      return(sprintf("%s: %s differs", tag, v))
  }
  ## compare the value as text: R's own parser is not always correctly
  ## rounded for extreme exponents, so a string check is the exact test
  if (!identical(cb$rps_value, sprintf("%.17g", rb$rps_value)))
    return(sprintf("%s: rps_value differs", tag))
  if ("CHROM" %in% names(rb) && "CHROM" %in% names(cb))
    if (!identical(as.character(rb$CHROM), as.character(cb$CHROM)))
      return(sprintf("%s: CHROM differs", tag))
  NA_character_
}

fails <- character(0); ntest <- 0L
gen <- function(n, nchr, tie_frac, sig_frac, span) {
  chrs <- sample(seq_len(nchr), n, replace = TRUE)
  pos  <- sample(seq_len(span), n, replace = TRUE)
  if (tie_frac > 0) {                       # inject duplicate positions
    k <- max(1L, floor(n * tie_frac)); ii <- sample(n, k)
    pos[ii] <- pos[sample(n, k)]
  }
  val <- as.numeric(sprintf("%.17g", runif(n)))
  hit <- sample(n, max(1L, floor(n * sig_frac)))
  val[hit] <- as.numeric(sprintf("%.17g", 10^(-runif(length(hit), 8, 30))))
  d <- data.frame(CHROM = chrs, position = pos, value = val)
  d[order(d$CHROM, d$position), ]
}

for (n in c(0L, 1L, 5L, 200L, 5000L)) {
 for (nchr in c(1L, 4L)) {
  for (reset in c("best", "any")) {
   for (reward in c("min", "max")) {
    for (win in c(1000, 50000, 500000)) {
      d <- if (n == 0L) data.frame(CHROM = integer(0), position = numeric(0), value = numeric(0))
           else gen(n, nchr, 0.1, 0.05, 2e6)
      if (n == 0L) next
      sig <- if (reward == "min") 5e-8 else 0.9
      if (reward == "max") d$value <- as.numeric(sprintf("%.17g", runif(nrow(d), 0, 1.2)))
      rb <- physical_merge(d, sig_th = sig, window = win, reward = reward,
                           reset_on = reset, chrom_col = "CHROM")
      cb <- run_c(d, sig, win, reward, reset, sorted = TRUE, multi = TRUE)
      tag <- sprintf("n=%d chr=%d %s/%s w=%g", n, nchr, reward, reset, win)
      r <- cmp(rb, cb, tag); ntest <- ntest + 1L
      if (!is.na(r)) fails <- c(fails, r)
      ## unsorted input through --sort
      ds <- d[sample(nrow(d)), ]
      rb2 <- physical_merge(ds, sig_th = sig, window = win, reward = reward,
                            reset_on = reset, chrom_col = "CHROM")
      cb2 <- run_c(ds, sig, win, reward, reset, sorted = FALSE, multi = TRUE)
      r2 <- cmp(rb2, cb2, paste(tag, "[--sort]")); ntest <- ntest + 1L
      if (!is.na(r2)) fails <- c(fails, r2)
    }
   }
  }
 }
}
## single-chromosome (no CHROM column) path
for (rep in 1:20) {
  d <- gen(1000, 1L, 0.05, 0.08, 5e6); d$CHROM <- NULL
  win <- sample(c(1e3, 1e4, 1e5, 1e6), 1); reset <- sample(c("best","any"),1)
  rb <- physical_merge(d, sig_th = 5e-8, window = win, reward = "min", reset_on = reset)
  cb <- run_c(d, 5e-8, win, "min", reset, sorted = TRUE, multi = FALSE)
  r <- cmp(rb, cb, sprintf("nochrom rep=%d w=%g %s", rep, win, reset)); ntest <- ntest + 1L
  if (!is.na(r)) fails <- c(fails, r)
}
cat(sprintf("\n%d comparisons, %d mismatches\n", ntest, length(fails)))
if (length(fails)) { cat(paste0("  ", fails, collapse = "\n"), "\n"); quit(status = 1) }
cat("ALL MATCH\n")
