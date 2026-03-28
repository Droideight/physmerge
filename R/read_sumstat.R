#' Read and prepare GWAS summary statistics for physical merging
#'
#' A unified reader that handles PLINK2 summary statistics, REML-GPCM per-SNP
#' output, and arbitrary tabular formats.  Uses \code{data.table::fread} for
#' fast, separator-agnostic reading.
#'
#' @param path      Path to the file (.csv, .tsv, .txt, or compressed).
#' @param format    One of \code{"plink2"}, \code{"gpcm"}, or \code{"custom"}.
#'   \describe{
#'     \item{\code{"plink2"}}{PLINK2 \code{.glm.*} output.  Default columns:
#'       \code{#CHROM}, \code{POS}, \code{ID}, \code{P}.}
#'     \item{\code{"gpcm"}}{REML-GPCM per-SNP output.  Default columns:
#'       \code{#CHROM}, \code{POS}, \code{ID}, \code{P_HPI}.}
#'     \item{\code{"custom"}}{Specify all column names manually via
#'       \code{chrom_col}, \code{pos_col}, \code{id_col}, \code{value_col}.}
#'   }
#' @param value_col Name of the column to use as the merging value.  Overrides
#'   the format default when supplied.  Common choices: \code{"P"},
#'   \code{"LOG10_P"}, \code{"P_HPI"}, \code{"P_Direct"}, \code{"T_STAT"}.
#' @param chrom_col Chromosome column name.  Overrides format default.
#' @param pos_col   Position column name.  Overrides format default.
#' @param id_col    SNP ID column name.  Overrides format default.
#'   Set to \code{NA} if no ID column exists.
#' @param test_filter Logical.  If \code{TRUE} (default for
#'   \code{format = "plink2"}), filter rows by \code{test_col == test_val}
#'   before returning.  Useful for multi-covariate PLINK2 output where each
#'   SNP appears once per TEST value.
#' @param test_col  Name of the TEST column.  Default \code{"TEST"}.
#' @param test_val  Value to retain in \code{test_col}.  Default \code{"ADD"}.
#' @param chrom     Optional character/integer vector of chromosomes to retain.
#'   \code{NULL} (default) keeps all.
#' @param ...       Additional arguments passed to \code{data.table::fread}
#'   (e.g. \code{nThread}, \code{skip}).
#'
#' @return A data frame containing all original columns plus:
#' \describe{
#'   \item{\code{position}}{Numeric copy of the position column.}
#'   \item{\code{value}}{Numeric copy of the value column.}
#' }
#' Rows where \code{position} or \code{value} is \code{NA} are dropped.
#' The data frame is sorted by \code{position}.
#'
#' An attribute \code{"suggested_reward"} is attached: \code{"max"} if
#' \code{value_col} is \code{"LOG10_P"} or a test statistic, \code{"min"}
#' otherwise.  Pass this to \code{\link{physical_merge}} as \code{reward}.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # PLINK2 sumstats (TEST filter auto-enabled)
#' df <- read_sumstat("my_gwas.glm.linear", format = "plink2")
#' blocks <- physical_merge(df, sig_th = 5e-8, window = 500000,
#'                          reward = attr(df, "suggested_reward"))
#'
#' # REML-GPCM output
#' df <- read_sumstat("stage1_ch1_P_HPI.csv", format = "gpcm",
#'                    value_col = "P_HPI", chrom = 1)
#'
#' # Custom format
#' df <- read_sumstat("results.txt", format = "custom",
#'                    chrom_col = "CHR", pos_col = "BP",
#'                    id_col = "SNP", value_col = "P")
#' }
read_sumstat <- function(path,
                         format      = c("plink2", "gpcm", "custom"),
                         value_col   = NULL,
                         chrom_col   = NULL,
                         pos_col     = NULL,
                         id_col      = NULL,
                         test_filter = NULL,
                         test_col    = "TEST",
                         test_val    = "ADD",
                         chrom       = NULL,
                         ...) {

  format <- match.arg(format)

  # ── Format-specific defaults ──────────────────────────────────────────────────
  defaults <- list(
    plink2 = list(chrom = "#CHROM", pos = "POS", id = "ID",
                  value = "P",     test_filter = TRUE),
    gpcm   = list(chrom = "#CHROM", pos = "POS", id = "ID",
                  value = "P_HPI", test_filter = FALSE),
    custom = list(chrom = NULL,     pos = NULL,  id = NULL,
                  value = NULL,     test_filter = FALSE)
  )[[format]]

  chrom_col   <- chrom_col   %||% defaults$chrom
  pos_col     <- pos_col     %||% defaults$pos
  id_col      <- id_col      %||% defaults$id
  value_col   <- value_col   %||% defaults$value
  test_filter <- test_filter %||% defaults$test_filter

  if (is.null(chrom_col) || is.null(pos_col) || is.null(value_col))
    stop("For format = 'custom', you must supply chrom_col, pos_col, and value_col.")

  # ── Read (fread auto-detects separator) ───────────────────────────────────────
  df <- tryCatch(
    as.data.frame(data.table::fread(path, header = TRUE,
                                    stringsAsFactors = FALSE,
                                    data.table = FALSE, ...)),
    error = function(e) stop("Failed to read file: ", conditionMessage(e))
  )

  # Normalise #CHROM → CHROM everywhere
  names(df)[names(df) == "#CHROM"] <- "CHROM"
  if (chrom_col == "#CHROM") chrom_col <- "CHROM"

  # ── Validate columns ──────────────────────────────────────────────────────────
  needed <- c(chrom_col, pos_col, value_col)
  if (!is.na(id_col)) needed <- c(needed, id_col)
  missing <- setdiff(needed, names(df))
  if (length(missing) > 0L)
    stop("Column(s) not found: ", paste(missing, collapse = ", "))

  # ── TEST column filter ────────────────────────────────────────────────────────
  if (isTRUE(test_filter)) {
    if (!test_col %in% names(df)) {
      warning("test_col '", test_col, "' not found; TEST filter skipped.")
    } else {
      n_before <- nrow(df)
      df       <- df[as.character(df[[test_col]]) == test_val, ]
      message(sprintf("TEST filter: kept %d of %d rows where %s = '%s'.",
                      nrow(df), n_before, test_col, test_val))
      if (nrow(df) == 0L)
        stop("No rows remain after TEST filter.")
    }
  }

  # ── Chromosome filter ─────────────────────────────────────────────────────────
  if (!is.null(chrom)) {
    df <- df[as.character(df[[chrom_col]]) %in% as.character(chrom), ]
    if (nrow(df) == 0L) warning("No rows remain after chromosome filter.")
  }

  # ── Append interface columns ──────────────────────────────────────────────────
  df$position <- suppressWarnings(as.numeric(df[[pos_col]]))
  df$value    <- suppressWarnings(as.numeric(df[[value_col]]))

  ok <- !is.na(df$position) & !is.na(df$value)
  if (any(!ok))
    message(sum(!ok), " row(s) dropped (NA in position or value).")
  df <- df[ok, ]

  if (nrow(df) == 0L)
    stop("No usable rows after filtering.")

  df <- df[order(df$position), ]

  # ── Suggest reward direction ──────────────────────────────────────────────────
  stat_cols <- c("LOG10_P", "T_STAT", "Z_STAT", "CHISQ", "F_STAT",
                 "T_STAT_Direct", "T_STAT_TE", "HPI")
  suggested <- if (value_col %in% stat_cols) "max" else "min"
  if (value_col == "LOG10_P")
    message("LOG10_P detected: consider reward = 'max' for physical_merge().")
  list(
    data   = df,
    reward = suggested
  )
}

# Null-coalescing operator (internal)
`%||%` <- function(x, y) if (is.null(x)) y else x
