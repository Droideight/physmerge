.is_significant <- function(val, sig_th, reward) {
  if (reward == "min") val < sig_th else val > sig_th
}

.is_more_significant <- function(val, best, reward) {
  if (reward == "min") val < best else val > best
}

#' Physical locus merging
#'
#' Scans a position-sorted genomic data frame and collapses nearby significant
#' signals into non-overlapping locus blocks using a forward sliding-window
#' approach.  No LD reference panel is required.
#'
#' @details
#' The algorithm works in three passes:
#'
#' \enumerate{
#'   \item \strong{Forward scan}: opens a block when a significant SNP is
#'     encountered and keeps it alive as long as the window has not been
#'     exhausted.  The behaviour when a significant (but not necessarily more
#'     significant) SNP is found inside the block is controlled by
#'     \code{reset_on}:
#'     \describe{
#'       \item{\code{"best"} (default)}{Steps reset only when a \emph{more}
#'         significant SNP is found.  The representative SNP is always the
#'         local maximum.}
#'       \item{\code{"any"}}{Steps reset whenever \emph{any} significant SNP
#'         is found, regardless of its value.  This is equivalent to taking
#'         the union of \eqn{\pm}\code{window} intervals around every
#'         significant SNP (i.e. locusDefiner-style logic).}
#'     }
#'   \item \strong{Collapse pass}: merges adjacent blocks whose representative
#'     SNPs (\code{rps_BP}) are fewer than \code{window} bp apart, retaining
#'     the more significant representative.
#'   \item \strong{Trim pass}: if after collapsing any block's \code{end}
#'     still overlaps the next block's \code{start}, the \code{end} is trimmed
#'     to \code{start} of the next block, guaranteeing zero overlap.
#' }
#'
#' When the input contains multiple chromosomes (detected via \code{chrom_col}
#' or an existing \code{CHROM} column), the algorithm is run independently per
#' chromosome to prevent cross-boundary merges.
#'
#' @param data      A data frame with (at least) two numeric columns:
#'   \describe{
#'     \item{\code{position}}{Base-pair coordinate, resorted internally.}
#'     \item{\code{value}}{Test statistic or p-value.}
#'   }
#' @param sig_th    Significance threshold (length-1 numeric).
#' @param window    Window size in base-pairs (positive numeric).
#' @param reward    \code{"min"} (default) for p-values; \code{"max"} for
#'   test statistics.
#' @param reset_on  \code{"best"} (default): steps reset only when a more
#'   significant SNP is encountered inside the current block.
#'   \code{"any"}: steps reset whenever any significant SNP is encountered,
#'   equivalent to the union-of-intervals logic used by locusDefiner.
#' @param chrom_col Name of the chromosome column in \code{data}.  If
#'   \code{NULL} (default), the function auto-detects a column named
#'   \code{"CHROM"}.  When a chromosome column is found and contains more
#'   than one unique value, the algorithm runs per chromosome.
#'
#' @return A data frame with one row per merged locus block:
#' \describe{
#'   \item{\code{serial}}{Sequential block index (1, 2, 3, …).}
#'   \item{\code{CHROM}}{Chromosome (present when a chromosome column is
#'     detected).}
#'   \item{\code{start}}{Block start in bp.}
#'   \item{\code{end}}{Block end in bp.}
#'   \item{\code{rps_BP}}{Position of the most significant representative SNP.}
#'   \item{\code{rps_value}}{Value of the representative SNP.}
#' }
#'
#' The returned data frame additionally carries an attribute
#' \code{"rps_row"}: the row index in \code{data} of each representative SNP.
#' \code{\link{annotate_blocks}} uses it to recover the exact input row, which
#' is the only reliable way to label a block when several variants share one
#' base-pair position (multi-allelic sites).  The attribute is not a column and
#' does not change the visible output.
#'
#' @export
#'
#' @examples
#' df <- data.frame(
#'   position = c(100, 200, 350, 5000, 5100, 5200, 9000),
#'   value    = c(0.04, 0.001, 0.03, 0.5, 0.02, 0.008, 0.04)
#' )
#' # default: reset only on more significant SNP
#' physical_merge(df, sig_th = 0.05, window = 500, reward = "min")
#'
#' # locusDefiner-equivalent: reset on any significant SNP
#' physical_merge(df, sig_th = 0.05, window = 500, reward = "min",
#'                reset_on = "any")
physical_merge <- function(data, sig_th, window, reward = "min",
                           reset_on = "best", chrom_col = NULL) {

  if (!is.data.frame(data))
    stop("`data` must be a data frame.")
  if (!all(c("position", "value") %in% names(data)))
    stop("`data` must contain columns named 'position' and 'value'.")
  if (!is.numeric(data$position) || !is.numeric(data$value))
    stop("Both 'position' and 'value' columns must be numeric.")
  if (!reward %in% c("min", "max"))
    stop("`reward` must be either 'min' or 'max'.")
  if (length(sig_th) != 1L || !is.numeric(sig_th) || is.na(sig_th))
    stop("`sig_th` must be a single numeric value.")
  if (length(window) != 1L || !is.numeric(window) || window <= 0)
    stop("`window` must be a single positive numeric value.")
  if (!reset_on %in% c("best", "any"))
    stop("`reset_on` must be either 'best' or 'any'.")

  chcol <- if (!is.null(chrom_col)) {
    if (!chrom_col %in% names(data))
      stop("chrom_col '", chrom_col, "' not found in data.")
    chrom_col
  } else if ("CHROM" %in% names(data)) {
    "CHROM"
  } else {
    NULL
  }

  keep <- !is.na(data$position) & !is.na(data$value)
  if (!is.null(chcol)) keep <- keep & !is.na(data[[chcol]])
  if (any(!keep)) {
    message(sum(!keep), " row(s) dropped (NA in position, value",
            if (!is.null(chcol)) " or chromosome", ").")
    data <- data[keep, , drop = FALSE]
  }
  oidx <- which(keep)

  if (nrow(data) == 0L) {
    empty <- data.frame(serial = integer(0), start = numeric(0), end = numeric(0),
                        rps_BP = numeric(0), rps_value = numeric(0),
                        rps_row = integer(0))
    if (!is.null(chcol))
      empty <- cbind(empty[, "serial", drop = FALSE], CHROM = character(0),
                     empty[, setdiff(names(empty), "serial"), drop = FALSE])
    return(.stash_rps_row(empty))
  }

  if (any(data$position < 0))
    warning("Negative position(s) found.  Block boundaries are clamped at 0, ",
            "so a block that starts below 0 has no width.")

  if (!is.null(chcol)) {
    chroms <- unique(data[[chcol]])
    if (length(chroms) > 1L) {
      res <- lapply(chroms, function(ch) {
        sel <- which(data[[chcol]] == ch)
        sub <- data[sel, ]
        blk <- .physical_merge_single(sub, sig_th, window, reward, reset_on)
        if (nrow(blk) == 0L) return(blk)
        blk$rps_row <- sel[blk$rps_row]
        blk$CHROM <- ch
        blk
      })
      out <- do.call(rbind, res)
      if (is.null(out) || nrow(out) == 0L) {
        return(.stash_rps_row(data.frame(serial = integer(0), CHROM = character(0),
                                         start = numeric(0), end = numeric(0),
                                         rps_BP = numeric(0), rps_value = numeric(0),
                                         rps_row = integer(0))))
      }
      out$serial <- seq_len(nrow(out))
      rownames(out) <- NULL
      cord <- c("serial", "CHROM",
                         setdiff(names(out), c("serial", "CHROM")))
      return(.stash_rps_row(.remap_rps_row(out[, cord], oidx)))
    }
  } else {
    rng <- diff(range(data$position))
    if (rng > 2.5e8)
      warning("Position range > 250 Mb detected but no chromosome column found. ",
              "If data spans multiple chromosomes, SNPs near chromosome ",
              "boundaries may be incorrectly merged into the same block. ",
              "Add a CHROM column or filter to one chromosome at a time.")
  }

  .stash_rps_row(.remap_rps_row(
    .physical_merge_single(data, sig_th, window, reward, reset_on), oidx))
}

.remap_rps_row <- function(blk, idx) {
  if (!is.null(blk$rps_row) && length(blk$rps_row))
    blk$rps_row <- idx[blk$rps_row]
  blk
}

.stash_rps_row <- function(blk) {
  rr <- blk$rps_row
  blk$rps_row <- NULL
  attr(blk, "rps_row") <- if (is.null(rr)) integer(0) else as.integer(rr)
  blk
}

.physical_merge_single <- function(data, sig_th, window, reward, reset_on) {

  ord <- order(data$position)
  data <- data[ord, ]
  n <- nrow(data)

  empt <- data.frame(
    serial = integer(0), start = numeric(0), end = numeric(0),
    rps_BP = numeric(0), rps_value = numeric(0), rps_row = integer(0)
  )
  if (n == 0L) return(empt)

  oser <- integer(n); ostart <- numeric(n)
  oend <- numeric(n); obp <- numeric(n)
  oval <- numeric(n); orow <- integer(n)
  nblk <- 0L

  inblk <- FALSE
  steps <- window
  best <- sig_th
  lpos <- data$position[1L]

  openb <- function(pos, val, i) {
    nblk <<- nblk + 1L
    oser[nblk] <<- nblk
    ostart[nblk] <<- max(0, pos - window)
    oend[nblk] <<- NA_real_
    obp[nblk] <<- pos
    oval[nblk] <<- val
    orow[nblk] <<- ord[i]
    inblk <<- TRUE
    steps <<- window
    best <<- val
  }

  closeb <- function(lastp) {
    oend[nblk] <<- lastp + steps
    inblk <<- FALSE
    steps <<- window
    best <<- sig_th
  }

  for (i in seq_len(n)) {
    pos <- data$position[i]
    val <- data$value[i]

    if (!inblk) {
      if (.is_significant(val, sig_th, reward)) openb(pos, val, i)

    } else {
      rem <- steps - (pos - lpos)

      if (rem <= 0) {
        closeb(lpos)
        if (.is_significant(val, sig_th, reward)) openb(pos, val, i)

      } else {
        steps <- rem

        if (reset_on == "any" && .is_significant(val, sig_th, reward)) {
          steps <- window
          if (.is_more_significant(val, best, reward)) {
            best <- val
            obp[nblk] <- pos
            oval[nblk] <- val
            orow[nblk] <- ord[i]
          }

        } else if (.is_more_significant(val, best, reward)) {
          best <- val
          steps <- window
          obp[nblk] <- pos
          oval[nblk] <- val
          orow[nblk] <- ord[i]
        }
      }
    }
    lpos <- pos
  }
  if (inblk) closeb(lpos)
  if (nblk == 0L) return(empt)

  raw <- data.frame(
    serial = oser[seq_len(nblk)],
    start = ostart[seq_len(nblk)],
    end = oend[seq_len(nblk)],
    rps_BP = obp[seq_len(nblk)],
    rps_value = oval[seq_len(nblk)],
    rps_row = orow[seq_len(nblk)],
    stringsAsFactors = FALSE
  )

  blk <- .collapse_blocks(raw, window, reward)

  if (nrow(blk) > 1L) {
    for (i in seq_len(nrow(blk) - 1L)) {
      if (blk$end[i] > blk$start[i + 1L])
        blk$end[i] <- blk$start[i + 1L]
    }
  }

  blk$end <- pmax(blk$end, blk$start)
  blk
}

.collapse_blocks <- function(blk, w, reward) {
  if (nrow(blk) <= 1L) return(blk)
  out <- blk[1L, ]
  for (i in seq(2L, nrow(blk))) {
    cur <- blk[i, ]
    if ((cur$rps_BP - out$rps_BP[nrow(out)]) < w) {
      last <- nrow(out)
      out$end[last] <- max(out$end[last], cur$end)
      if (.is_more_significant(cur$rps_value, out$rps_value[last], reward)) {
        out$rps_BP[last] <- cur$rps_BP
        out$rps_value[last] <- cur$rps_value
        out$rps_row[last] <- cur$rps_row
      }
    } else {
      out <- rbind(out, cur)
    }
  }
  out$serial <- seq_len(nrow(out))
  rownames(out) <- NULL
  out
}
