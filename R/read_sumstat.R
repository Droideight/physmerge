#' Read GWAS summary statistics for physical merging
#'
#' Reads a summary-statistics file and adds the role columns that
#' \code{\link{physical_merge}} works on: \code{position} and \code{value},
#' plus \code{chrom} and \code{id} when the file has them. All original
#' columns are kept unchanged. Rows keep the order of the file.
#'
#' @param path Path to the file. Plain text or gzip (\code{.gz}).
#' @param format \code{"plink2"} (default) or \code{"custom"}.
#' \describe{
#' \item{\code{"plink2"}}{PLINK2 \code{.glm.*} output. Columns \code{#CHROM},
#' \code{POS}, \code{ID} and \code{P}; only rows whose \code{TEST} is
#' \code{test_val} are kept.}
#' \item{\code{"custom"}}{Any other table. \code{pos_col}, \code{value_col}
#' and \code{chrom_col} must be named.}
#' }
#' @param chrom_col Chromosome column. \code{NA} ignores the chromosome
#' column, so the whole file is merged as one sequence.
#' @param pos_col Position column.
#' @param id_col SNP id column. \code{NA} for none.
#' @param value_col Column holding the p-value or statistic.
#' @param test_val For \code{format = "plink2"}, the \code{TEST} value to keep.
#' Default \code{"ADD"}; \code{"DOM"} or \code{"REC"} for dominant or
#' recessive models.
#' @param chrom Chromosomes to keep, e.g. \code{c(1, 2)}. \code{NULL} keeps all.
#' @param sep Field separator: \code{"tab"}, \code{"space"} or a single
#' character. \code{NULL} reads it from the header: tab if the header has a
#' tab, otherwise comma, otherwise semicolon, otherwise whitespace.
#' @param quiet \code{TRUE} hides progress messages. Warnings and errors are
#' still shown.
#'
#' @return A data frame: the original columns, then \code{position},
#' \code{value}, and \code{chrom} and \code{id} when present.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' d <- read_sumstat("my_gwas.glm.linear")
#' b <- physical_merge(d)
#'
#' d <- read_sumstat("results.txt", format = "custom",
#' chrom_col = "CHR", pos_col = "BP", id_col = "SNP", value_col = "P")
#' }
read_sumstat <- function(path,
                         format = c("plink2", "custom"),
                         chrom_col = NULL,
                         pos_col = NULL,
                         id_col = NULL,
                         value_col = NULL,
                         test_val = "ADD",
                         chrom = NULL,
                         sep = NULL,
                         quiet = FALSE) {

  format <- match.arg(format)
  say <- function(...) if (!quiet) message("physmerge: ", ...)

  if (format == "plink2") {
    chrom_col <- chrom_col %||% "#CHROM"
    pos_col <- pos_col %||% "POS"
    id_col <- id_col %||% "ID"
    value_col <- value_col %||% "P"
  }
  if (is.null(chrom_col) || is.null(pos_col) || is.null(value_col))
    stop("For format = 'custom', you must supply chrom_col, pos_col and ",
         "value_col; set chrom_col = NA to ignore the chromosome column.",
         call. = FALSE)
  if (is.null(id_col)) id_col <- NA_character_

  if (is.null(sep)) sep <- .pm_sep(path)
  else if (identical(sep, "tab")) sep <- "\t"
  else if (identical(sep, "space")) sep <- " "
  df <- tryCatch(
    as.data.frame(data.table::fread(path, header = TRUE, sep = sep,
                                    stringsAsFactors = FALSE,
                                    data.table = FALSE, check.names = FALSE)),
    error = function(e) stop("Failed to read file: ", conditionMessage(e),
                             call. = FALSE)
  )

  needed <- c(if (!is.na(chrom_col)) chrom_col, pos_col, value_col,
              if (!is.na(id_col)) id_col)
  missing <- setdiff(needed, names(df))
  if (length(missing) > 0L)
    stop("Column(s) not found: ", paste(missing, collapse = ", "), call. = FALSE)

  if (format == "plink2") {
    if (!"TEST" %in% names(df)) {
      warning("physmerge: TEST column not found; TEST filter skipped.",
              call. = FALSE)
    } else {
      tv <- as.character(df$TEST)
      n0 <- nrow(df)
      hit <- !is.na(tv) & tv == test_val
      if (!any(hit))
        stop(sprintf("no rows have TEST = '%s'; TEST values in this file: %s. Set test_val to one of them.",
                     test_val, paste(head(unique(tv[!is.na(tv)]), 10), collapse = ", ")),
             call. = FALSE)
      df <- df[hit, , drop = FALSE]
      say(sprintf("TEST filter: kept %d of %d rows where TEST = '%s'.",
                  nrow(df), n0, test_val))
    }
  }

  if (!is.null(chrom) && !is.na(chrom_col)) {
    hit <- as.character(df[[chrom_col]]) %in% as.character(chrom)
    if (any(!hit)) say(sum(!hit), " row(s) dropped by the chromosome filter.")
    df <- df[hit, , drop = FALSE]
  }

  df$position <- suppressWarnings(as.numeric(df[[pos_col]]))
  df$value <- suppressWarnings(as.numeric(df[[value_col]]))
  if (!is.na(chrom_col)) df$chrom <- as.character(df[[chrom_col]])
  if (!is.na(id_col)) df$id <- as.character(df[[id_col]])

  ok <- !is.na(df$position) & !is.na(df$value)
  if (any(!ok)) say(sum(!ok), " row(s) dropped (NA in position or value).")
  df <- df[ok, , drop = FALSE]
  rownames(df) <- NULL

  if (nrow(df) == 0L)
    stop("No usable rows after filtering.", call. = FALSE)
  attr(df, "pm_cols") <- c(chrom = chrom_col, pos = pos_col, id = id_col,
                           value = value_col)
  df
}

.pm_sep <- function(path) {
  con <- if (grepl("\\.gz$", path)) gzfile(path, "r") else file(path, "r")
  on.exit(close(con))
  h <- readLines(con, n = 1L, warn = FALSE)
  if (!length(h)) stop("empty input file", call. = FALSE)
  if (grepl("\t", h, fixed = TRUE)) "\t"
  else if (grepl(",", h, fixed = TRUE)) ","
  else if (grepl(";", h, fixed = TRUE)) ";"
  else " "
}

`%||%` <- function(x, y) if (is.null(x)) y else x
