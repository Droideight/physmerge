#' Physical locus merging
#'
#' Scans position-sorted summary statistics and collapses nearby significant
#' SNPs into contiguous, strictly non-overlapping locus blocks. No LD reference
#' panel is required.
#'
#' @details
#' A block opens at the first significant SNP, one window upstream of it. The
#' block carries a window-sized budget that the distance travelled spends, and
#' every significant SNP refills it, which gives the union of the
#' \eqn{\pm}\code{window} intervals around all significant SNPs. When the
#' budget runs out, the
#' block closes one window downstream of its last significant SNP. The end of
#' each block is then trimmed back to the start of the next one. With a
#' \code{chrom} column the merge runs per chromosome, in order of first
#' appearance.
#'
#' The reward direction and the threshold follow the name of the column the
#' values came from, when \code{reward} or \code{sig_th} is not given:
#' \code{P}, \code{PVAL}, \code{PVALUE}, \code{P_VALUE}, \code{P-VALUE} and any
#' name starting with \code{P_} give \code{"min"} and \code{5e-8};
#' \code{NEG_LOG10_P} and \code{LOG10_P} give \code{"max"} and
#' \code{-log10(5e-8)}. Any other column needs both set by hand.
#'
#' @param data A data frame with numeric columns \code{position} and
#' \code{value}, and optionally \code{chrom} and \code{id}, as returned by
#' \code{\link{read_sumstat}}.
#' @param sig_th Significance threshold. \code{NULL} sets it from the value
#' column's name.
#' @param window Window size in base pairs. Default 500000.
#' @param reward \code{"min"} when smaller values are more significant
#' (p-values), \code{"max"} when larger ones are. \code{NULL} sets it from the
#' value column's name.
#' @param sort \code{FALSE} (default) stops on input that is not
#' position-sorted within each chromosome. \code{TRUE} sorts it first.
#' @param quiet \code{TRUE} hides progress messages. Warnings and errors are
#' still shown.
#'
#' @return A data frame with one row per block: \code{serial}, \code{CHROM}
#' (when \code{data} has a \code{chrom} column), \code{start}, \code{end},
#' \code{rps_BP}, \code{rps_ID} (when \code{data} has an \code{id} column) and
#' the representative's value, named \code{rps_} followed by the name of the
#' value column, e.g. \code{rps_P}. The attribute \code{"rps_row"} holds the
#' row of \code{data} each representative came from.
#'
#' @export
#'
#' @examples
#' df <- data.frame(
#' position = c(100, 200, 350, 5000, 5100, 5200, 9000),
#' value = c(0.04, 0.001, 0.03, 0.5, 0.02, 0.008, 0.04)
#' )
#' physical_merge(df, sig_th = 0.05, window = 500, reward = "min")
physical_merge <- function(data, sig_th = NULL, window = 500000, reward = NULL,
                           sort = FALSE, quiet = FALSE) {

  .pm_flag(quiet, "quiet")
  .pm_flag(sort, "sort")
  say <- function(...) if (!quiet) message("physmerge: ", ...)

  if (!is.data.frame(data))
    stop("`data` must be a data frame.", call. = FALSE)
  if (!all(c("position", "value") %in% names(data)))
    stop("`data` must contain columns named 'position' and 'value'.", call. = FALSE)
  if (!is.numeric(data$position) || !is.numeric(data$value))
    stop("Both 'position' and 'value' columns must be numeric.", call. = FALSE)
  if (length(window) != 1L || !is.numeric(window) || is.na(window) || window <= 0)
    stop("`window` must be a single positive numeric value.", call. = FALSE)

  vname <- .pm_vname(data)
  kind <- .pm_kind(vname)
  if ((is.null(reward) || is.null(sig_th)) && kind == "other")
    stop(sprintf("cannot tell whether a larger or a smaller %s is more significant; set reward and sig_th.",
                 vname), call. = FALSE)
  if (is.null(reward)) reward <- if (kind == "p") "min" else "max"
  if (is.null(sig_th)) sig_th <- if (kind == "p") 5e-8 else -log10(5e-8)
  if (!identical(reward, "min") && !identical(reward, "max"))
    stop("`reward` must be either 'min' or 'max'.", call. = FALSE)
  if (length(sig_th) != 1L || !is.numeric(sig_th) || is.na(sig_th))
    stop("`sig_th` must be a single numeric value.", call. = FALSE)

  hasch <- "chrom" %in% names(data)
  hasid <- "id" %in% names(data)
  vcol <- paste0("rps_", vname)

  keep <- !is.na(data$position) & !is.na(data$value)
  if (any(!keep)) say(sum(!keep), " row(s) dropped (NA in position or value).")
  rows <- which(keep)
  pos <- data$position[rows]
  val <- data$value[rows]
  ch <- if (hasch) as.character(data$chrom[rows]) else rep("", length(rows))
  ch[is.na(ch)] <- ""
  n <- length(rows)

  if (n && any(pos < 0))
    warning("physmerge: negative position(s) found; block boundaries are clamped at 0.",
            call. = FALSE)
  if (n && !hasch && diff(range(pos)) > 2.5e8)
    warning("physmerge: positions span more than 250 Mb but there is no chromosome ",
            "column; blocks may cross chromosome boundaries.", call. = FALSE)

  if (n > 1L) {
    if (sort) {
      o <- order(match(ch, unique(ch)), pos)
      rows <- rows[o]; pos <- pos[o]; val <- val[o]; ch <- ch[o]
    } else {
      same <- c(FALSE, ch[-1L] == ch[-n])
      back <- same & c(FALSE, diff(pos) < 0)
      again <- !same & duplicated(ch)
      bad <- which(back | again)
      if (length(bad)) {
        i <- bad[1L]
        if (again[i])
          stop(sprintf("chromosome %s appears in more than one stretch of the input; set sort to merge it anyway.",
                       ch[i]), call. = FALSE)
        if (hasch)
          stop(sprintf("input is not position-sorted (chromosome %s: %.0f after %.0f); set sort to merge it anyway.",
                       ch[i], pos[i], pos[i - 1L]), call. = FALSE)
        stop(sprintf("input is not position-sorted (%.0f after %.0f); set sort to merge it anyway.",
                     pos[i], pos[i - 1L]), call. = FALSE)
      }
    }
  }

  rmax <- reward == "max"
  grp <- split(seq_len(n), factor(ch, levels = unique(ch)))
  parts <- lapply(seq_along(grp), function(j) {
    sel <- grp[[j]]
    cc <- names(grp)[j]
    b <- .pm_merge1(pos[sel], val[sel], sig_th, window, rmax)
    b$k <- sel[b$k]
    b$chrom <- rep(cc, length(b$k))
    b
  })
  start <- unlist(lapply(parts, `[[`, "start"))
  end <- unlist(lapply(parts, `[[`, "end"))
  bp <- unlist(lapply(parts, `[[`, "bp"))
  bv <- unlist(lapply(parts, `[[`, "val"))
  k <- unlist(lapply(parts, `[[`, "k"))
  bch <- unlist(lapply(parts, `[[`, "chrom"))
  if (is.null(k)) {
    start <- end <- bp <- bv <- numeric(0); k <- integer(0); bch <- character(0)
  }

  out <- data.frame(serial = seq_along(k), stringsAsFactors = FALSE)
  if (hasch) out$CHROM <- bch
  out$start <- start
  out$end <- end
  out$rps_BP <- bp
  if (hasid) out$rps_ID <- as.character(data$id[rows[k]])
  out[[vcol]] <- bv
  attr(out, "rps_row") <- as.integer(rows[k])

  say(sprintf("%d SNPs -> %d blocks (window=%.15g, sig_th=%g, reward=%s).",
              n, nrow(out), window, sig_th, reward))
  out
}

.pm_merge1 <- function(pos, val, sig, w, rmax) {
  issig <- if (rmax) function(v) v > sig else function(v) v < sig
  better <- if (rmax) function(v, b) v > b else function(v, b) v < b
  n <- length(pos)
  s <- e <- p <- v <- numeric(n); k <- integer(n); nb <- 0L
  inblk <- FALSE; steps <- w; best <- sig; lpos <- if (n) pos[1L] else 0
  for (i in seq_len(n)) {
    x <- pos[i]; y <- val[i]
    if (inblk) {
      rem <- steps - (x - lpos)
      if (rem <= 0) {
        e[nb] <- lpos + steps
        inblk <- FALSE
      } else {
        steps <- rem
        if (issig(y)) {
          steps <- w
          if (better(y, best)) { best <- y; p[nb] <- x; v[nb] <- y; k[nb] <- i }
        }
      }
    }
    if (!inblk && issig(y)) {
      nb <- nb + 1L
      s[nb] <- max(0, x - w); p[nb] <- x; v[nb] <- y; k[nb] <- i
      inblk <- TRUE; steps <- w; best <- y
    }
    lpos <- x
  }
  if (inblk) e[nb] <- lpos + steps
  s <- s[seq_len(nb)]; e <- e[seq_len(nb)]; p <- p[seq_len(nb)]
  v <- v[seq_len(nb)]; k <- k[seq_len(nb)]

  if (nb > 1L) {
    keep <- 1L
    for (i in 2:nb) {
      j <- keep[length(keep)]
      if (p[i] - p[j] < w) {
        e[j] <- max(e[j], e[i])
        if (better(v[i], v[j])) { p[j] <- p[i]; v[j] <- v[i]; k[j] <- k[i] }
      } else {
        keep <- c(keep, i)
      }
    }
    s <- s[keep]; e <- e[keep]; p <- p[keep]; v <- v[keep]; k <- k[keep]
    m <- length(s)
    if (m > 1L) e[-m] <- pmin(e[-m], s[-1L])
  }
  e <- pmax(e, s)
  list(start = s, end = e, bp = p, val = v, k = k)
}

.pm_vname <- function(data) {
  pc <- attr(data, "pm_cols")
  if (!is.null(pc) && !is.na(pc[["value"]])) return(pc[["value"]])
  v <- as.numeric(data$value)
  for (nm in setdiff(names(data), c("position", "value", "chrom", "id"))) {
    x <- data[[nm]]
    if (is.numeric(x) && identical(as.numeric(x), v)) return(nm)
  }
  "value"
}

.pm_kind <- function(name) {
  u <- toupper(name)
  if (u %in% c("P", "PVAL", "PVALUE", "P_VALUE", "P-VALUE") || startsWith(u, "P_"))
    "p"
  else if (u %in% c("NEG_LOG10_P", "LOG10_P"))
    "logp"
  else
    "other"
}
