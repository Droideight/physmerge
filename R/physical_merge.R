# Internal helpers — not exported, not documented publicly

.is_significant <- function(val, sig_th, reward) {
  if (reward == "min") val < sig_th else val > sig_th
}

.is_more_significant <- function(val, best, reward) {
  if (reward == "min") val < best else val > best
}


# ==============================================================================

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
#'     encountered and keeps it alive as long as the next SNP falls within
#'     \code{window} bp of the previous one.  Whenever a more significant SNP
#'     is found inside the block, it becomes the new representative and the
#'     window resets.
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
#'     \item{\code{position}}{Base-pair coordinate.}
#'     \item{\code{value}}{Test statistic or p-value.}
#'   }
#' @param sig_th    Significance threshold (length-1 numeric).
#' @param window    Window size in base-pairs (positive numeric).
#' @param reward    \code{"min"} (default) for p-values; \code{"max"} for
#'   test statistics.
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
#' @export
#'
#' @examples
#' df <- data.frame(
#'   position = c(100, 200, 350, 5000, 5100, 5200, 9000),
#'   value    = c(0.04, 0.001, 0.03, 0.5, 0.02, 0.008, 0.04)
#' )
#' physical_merge(df, sig_th = 0.05, window = 500, reward = "min")
physical_merge <- function(data, sig_th, window, reward = "min",
                           chrom_col = NULL) {
  
  # ── Input validation ─────────────────────────────────────────────────────────
  if (!is.data.frame(data))
    stop("`data` must be a data frame.")
  if (!all(c("position", "value") %in% names(data)))
    stop("`data` must contain columns named 'position' and 'value'.")
  if (!is.numeric(data$position) || !is.numeric(data$value))
    stop("Both 'position' and 'value' columns must be numeric.")
  if (!reward %in% c("min", "max"))
    stop("`reward` must be either 'min' or 'max'.")
  if (length(sig_th) != 1L || !is.numeric(sig_th))
    stop("`sig_th` must be a single numeric value.")
  if (length(window) != 1L || !is.numeric(window) || window <= 0)
    stop("`window` must be a single positive numeric value.")
  
  # ── Per-chromosome dispatch ───────────────────────────────────────────────────
  resolved_chrom <- if (!is.null(chrom_col)) {
    if (!chrom_col %in% names(data))
      stop("chrom_col '", chrom_col, "' not found in data.")
    chrom_col
  } else if ("CHROM" %in% names(data)) {
    "CHROM"
  } else {
    NULL
  }
  
  if (!is.null(resolved_chrom)) {
    chroms <- unique(data[[resolved_chrom]])
    if (length(chroms) > 1L) {
      results <- lapply(chroms, function(ch) {
        sub <- data[data[[resolved_chrom]] == ch, ]
        blk <- .physical_merge_single(sub, sig_th, window, reward)
        if (nrow(blk) == 0L) return(blk)
        blk$CHROM <- ch
        blk
      })
      out <- do.call(rbind, results)
      if (is.null(out) || nrow(out) == 0L) {
        return(data.frame(serial = integer(0), CHROM = character(0),
                          start = numeric(0), end = numeric(0),
                          rps_BP = numeric(0), rps_value = numeric(0)))
      }
      out$serial    <- seq_len(nrow(out))
      rownames(out) <- NULL
      col_order     <- c("serial", "CHROM",
                         setdiff(names(out), c("serial", "CHROM")))
      return(out[, col_order])
    }
  } else {
    # No chrom info: warn if range suggests multi-chromosomal data
    pos_range <- diff(range(data$position, na.rm = TRUE))
    if (pos_range > 2.5e8)
      warning("Position range > 250 Mb detected but no chromosome column found. ",
              "If data spans multiple chromosomes, SNPs near chromosome ",
              "boundaries may be incorrectly merged into the same block. ",
              "Add a CHROM column or filter to one chromosome at a time.")
  }
  
  .physical_merge_single(data, sig_th, window, reward)
}


# Internal: single-chromosome merging
.physical_merge_single <- function(data, sig_th, window, reward) {
  
  data <- data[order(data$position), ]
  n    <- nrow(data)
  
  empty_out <- data.frame(
    serial    = integer(0), start = numeric(0), end   = numeric(0),
    rps_BP    = numeric(0), rps_value = numeric(0)
  )
  if (n == 0L) return(empty_out)
  
  out_serial  <- integer(n);  out_start   <- numeric(n)
  out_end     <- numeric(n);  out_rps_bp  <- numeric(n)
  out_rps_val <- numeric(n);  block_count <- 0L
  
  in_block       <- FALSE
  steps          <- window
  sig_this_block <- sig_th
  last_pos       <- data$position[1L]
  
  open_block <- function(pos, val) {
    block_count <<- block_count + 1L
    out_serial[block_count]  <<- block_count
    out_start[block_count]   <<- max(0, pos - window)
    out_end[block_count]     <<- NA_real_
    out_rps_bp[block_count]  <<- pos
    out_rps_val[block_count] <<- val
    in_block       <<- TRUE
    steps          <<- window
    sig_this_block <<- val
  }
  
  close_block <- function(last_inblock_pos) {
    out_end[block_count] <<- last_inblock_pos + steps  # remaining steps, not full window
    in_block             <<- FALSE
    steps                <<- window
    sig_this_block       <<- sig_th
  }
  
  for (i in seq_len(n)) {
    pos <- data$position[i]
    val <- data$value[i]
    
    if (!in_block) {
      if (.is_significant(val, sig_th, reward)) open_block(pos, val)
      
    } else {
      remaining <- steps - (pos - last_pos)
      
      if (remaining <= 0) {
        close_block(last_pos)
        if (.is_significant(val, sig_th, reward)) open_block(pos, val)
        
      } else {
        steps <- remaining
        if (.is_more_significant(val, sig_this_block, reward)) {
          sig_this_block           <- val
          steps                    <- window
          out_rps_bp[block_count]  <- pos
          out_rps_val[block_count] <- val
        }
      }
    }
    last_pos <- pos
  }
  if (in_block) close_block(last_pos)
  if (block_count == 0L) return(empty_out)
  
  raw_blocks <- data.frame(
    serial    = out_serial[seq_len(block_count)],
    start     = out_start[seq_len(block_count)],
    end       = out_end[seq_len(block_count)],
    rps_BP    = out_rps_bp[seq_len(block_count)],
    rps_value = out_rps_val[seq_len(block_count)],
    stringsAsFactors = FALSE
  )
  
  blk <- .collapse_blocks(raw_blocks, window, reward)
  
  if (nrow(blk) > 1L) {
    for (i in seq_len(nrow(blk) - 1L)) {
      if (blk$end[i] > blk$start[i + 1L])
        blk$end[i] <- blk$start[i + 1L]
    }
  }
  
  blk
}


# Internal: collapse pass
.collapse_blocks <- function(blk, w, reward) {
  if (nrow(blk) <= 1L) return(blk)
  out <- blk[1L, ]
  for (i in seq(2L, nrow(blk))) {
    cur <- blk[i, ]
    if ((cur$rps_BP - out$rps_BP[nrow(out)]) < w) {
      last          <- nrow(out)
      out$end[last] <- max(out$end[last], cur$end)
      if (.is_more_significant(cur$rps_value, out$rps_value[last], reward)) {
        out$rps_BP[last]    <- cur$rps_BP
        out$rps_value[last] <- cur$rps_value
      }
    } else {
      out <- rbind(out, cur)
    }
  }
  out$serial    <- seq_len(nrow(out))
  rownames(out) <- NULL
  out
}