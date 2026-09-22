# Regression tests for the findings in tests/qa/FINDINGS.txt (2026-09-22).
# Each test names the finding it locks down.

# ---------------------------------------------------------------- read_sumstat

write_tsv <- function(lines) {
  f <- tempfile(fileext = ".tsv")
  writeLines(lines, f)
  f
}

test_that("R1: format='custom' works without an id_col", {
  f <- write_tsv(c("CHR\tBP\tPV", "1\t100\t1e-9", "1\t200\t1e-9"))
  out <- expect_silent(
    suppressMessages(read_sumstat(f, "custom", chrom_col = "CHR",
                                  pos_col = "BP", value_col = "PV")))
  expect_equal(nrow(out$data), 2L)
  expect_equal(out$reward, "min")
})

test_that("R2: rows come back grouped by chromosome, in first-appearance order", {
  f <- write_tsv(c("#CHROM\tPOS\tID\tP",
                   "2\t500\ta\t1e-9", "1\t100\tb\t1e-9",
                   "2\t900\tc\t1e-9", "1\t700\td\t1e-9"))
  d <- suppressMessages(read_sumstat(f, "gpcm", value_col = "P"))$data
  expect_equal(as.character(d$CHROM), c("2", "2", "1", "1"))
  expect_equal(d$position, c(500, 900, 100, 700))
})

# -------------------------------------------------------------- physical_merge

test_that("R3: an NA value drops its row instead of aborting the run", {
  d <- data.frame(position = c(1, 500000, 900000),
                  value    = c(1e-9, NA, 1e-9))
  expect_message(physical_merge(d, 5e-8, 5e5), "1 row\\(s\\) dropped")
  expect_equal(nrow(suppressMessages(physical_merge(d, 5e-8, 5e5))), 2L)
})

test_that("R3: an NA position drops its row, and sig_th = NA is rejected", {
  d <- data.frame(position = c(1, NA, 900000), value = rep(1e-9, 3))
  expect_equal(nrow(suppressMessages(physical_merge(d, 5e-8, 5e5))), 2L)
  expect_error(physical_merge(data.frame(position = 1, value = 1e-9), NA, 500),
               "single numeric value")
})

test_that("R4: an NA chromosome is dropped with a message, not silently", {
  d <- data.frame(CHROM = c("1", NA), position = c(1000, 2000),
                  value = c(1e-9, 1e-9))
  expect_message(physical_merge(d, 5e-8, 500), "1 row\\(s\\) dropped")
})

test_that("R3/R4: rps_row still points at the caller's rows after a drop", {
  d <- data.frame(CHROM = "1", ID = c("drop", "lead", "tail"),
                  position = c(100, 300, 500),
                  value    = c(NA, 1e-12, 1e-9))
  b <- suppressMessages(physical_merge(d, 5e-8, 150))
  expect_equal(attr(b, "rps_row"), c(2L, 3L))
  expect_equal(annotate_blocks(b, d)$rps_ID, c("lead", "tail"))
})

test_that("R8: an empty data frame merges without warnings", {
  expect_silent(physical_merge(data.frame(position = numeric(0),
                                          value = numeric(0)), 5e-8, 5e5))
})

test_that("R11: a block's end is never below its start", {
  d <- data.frame(position = c(-1000, -200), value = c(1e-9, 1e-9))
  expect_warning(physical_merge(d, 5e-8, 500), "Negative position")
  b <- suppressWarnings(physical_merge(d, 5e-8, 500))
  expect_true(all(b$end >= b$start))
})

# ------------------------------------------------------------- annotate_blocks

test_that("R9: turning every keep_* off does not error", {
  d <- data.frame(CHROM = "1", ID = c("a", "b"), POS = c(1000, 2000),
                  position = c(1000, 2000), value = c(1e-9, 1e-9))
  b <- physical_merge(d, 5e-8, 500)
  expect_error(
    annotate_blocks(b, d, keep_serial = FALSE, keep_start = FALSE,
                    keep_end = FALSE, keep_rps_BP = FALSE,
                    keep_rps_value = FALSE, keep_rps_ID = FALSE),
    NA)
})

test_that("R10: annotating with the wrong data frame warns", {
  d  <- data.frame(CHROM = "1", ID = "a", position = 1000, value = 1e-9)
  d2 <- data.frame(CHROM = "1", ID = "z", position = 9e6,  value = 1e-9)
  b  <- physical_merge(d, 5e-8, 500)
  attr(b, "rps_row") <- NULL
  expect_warning(annotate_blocks(b, d2), "No block matched")
})

# ------------------------------------------------------------- export_snp_list

test_that("R5: a numeric id is not written in scientific notation", {
  b <- data.frame(serial = 1:2, CHROM = "1", start = 0, end = 1,
                  rps_BP = c(1e6, 9e5), rps_value = 1e-9)
  f <- tempfile()
  suppressMessages(export_snp_list(b, f))
  expect_equal(readLines(f), c("1000000", "900000"))
})

test_that("R6: by_chrom with a relative path writes where the caller is", {
  b <- data.frame(serial = 1:2, CHROM = c("1", "2"), start = 0, end = 1,
                  rps_BP = c(1e6, 2e6), rps_value = 1e-9)
  dir <- tempfile(); dir.create(dir)
  old <- setwd(dir); on.exit(setwd(old), add = TRUE)
  suppressMessages(export_snp_list(b, "rel.zip", by_chrom = TRUE))
  expect_true(file.exists(file.path(dir, "rel.zip")))
  expect_setequal(utils::unzip(file.path(dir, "rel.zip"), list = TRUE)$Name,
                  c("snp_ch1.txt", "snp_ch2.txt"))
  expect_equal(normalizePath(getwd()), normalizePath(dir))
})

test_that("R7: a chromosome name cannot escape the archive", {
  b <- data.frame(serial = 1:2, CHROM = c("1/alt", "2"), start = 0, end = 1,
                  rps_BP = c(1e6, 2e6), rps_value = 1e-9)
  f <- tempfile(fileext = ".zip")
  suppressMessages(export_snp_list(b, f, by_chrom = TRUE))
  expect_setequal(utils::unzip(f, list = TRUE)$Name,
                  c("snp_ch1_alt.txt", "snp_ch2.txt"))

  bad <- b; bad$CHROM <- c("a.b", "a/b")
  expect_error(suppressMessages(export_snp_list(bad, tempfile(), by_chrom = TRUE)),
               "collide")
})
