#' Add the original columns of each representative SNP
#'
#' Appends to each block the original columns of its representative SNP. The
#' columns the block table already shows (chromosome, position, id and value)
#' are not repeated.
#'
#' @param blocks Data frame returned by \code{\link{physical_merge}}.
#' @param data The data frame that was passed to \code{\link{physical_merge}}.
#'
#' @return \code{blocks} with the representative's remaining original columns
#' appended, in the order of the input file.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' d <- read_sumstat("my_gwas.glm.linear")
#' b <- annotate_blocks(physical_merge(d), d)
#' b[, c("CHROM", "start", "end", "rps_ID", "BETA", "SE")]
#' }
annotate_blocks <- function(blocks, data) {

  if (!"position" %in% names(data) || !"value" %in% names(data))
    stop("`data` must contain columns named 'position' and 'value'.", call. = FALSE)

  roles <- c("position", "value", "chrom", "id")
  pc <- attr(data, "pm_cols")
  same <- function(x, y) {
    if (is.numeric(x) && is.numeric(y)) identical(as.numeric(x), as.numeric(y))
    else identical(as.character(x), as.character(y))
  }
  nm <- names(data)
  extra <- if (!is.null(pc)) {
    which(!nm %in% c(roles, pc[!is.na(pc)]))
  } else {
    Filter(function(j) {
      x <- data[[j]]
      !(nm[j] %in% roles || same(x, data$position) || same(x, data$value) ||
          ("chrom" %in% nm && same(x, data$chrom)) ||
          ("id" %in% nm && same(x, data$id)))
    }, seq_along(nm))
  }
  extra <- extra[!nm[extra] %in% names(blocks)]

  rr <- attr(blocks, "rps_row")
  ok <- !is.null(rr) && length(rr) == nrow(blocks) && !anyNA(rr) &&
    all(rr >= 1L) && all(rr <= nrow(data)) &&
    isTRUE(all.equal(as.numeric(data$position[rr]),
                     as.numeric(blocks$rps_BP), tolerance = 0))

  if (!ok) {
    key <- function(ch, p) paste(ch, vapply(p, format, "", digits = 15, scientific = FALSE, trim = TRUE))
    dk <- key(if ("chrom" %in% names(data)) data$chrom else "", data$position)
    bk <- key(if ("CHROM" %in% names(blocks)) blocks$CHROM else "", blocks$rps_BP)
    rr <- match(bk, dk)
    if (nrow(blocks) && all(is.na(rr)))
      warning("physmerge: no block matched a row of `data`; is this the data ",
              "frame the blocks were built from?", call. = FALSE)
    else if (anyNA(rr))
      warning("physmerge: ", sum(is.na(rr)), " of ", nrow(blocks),
              " block(s) matched no row of `data`.", call. = FALSE)
  }

  add <- data[rr, extra, drop = FALSE]
  names(add) <- nm[extra]
  rownames(add) <- NULL
  out <- cbind(blocks, add)
  attr(out, "rps_row") <- NULL
  out
}

#' Write the representative SNP ids to a file
#'
#' Writes one id per line: \code{rps_ID} when the blocks have it, otherwise
#' \code{rps_BP}. Give either \code{file} for a single list, or \code{dir} for
#' one file per chromosome (\code{snp_ch1.txt}, \code{snp_ch2.txt}, ...). The
#' directory is created if it does not exist.
#'
#' @param blocks Data frame returned by \code{\link{physical_merge}} or
#' \code{\link{annotate_blocks}}.
#' @param file Path of a single output file.
#' @param dir Directory for one file per chromosome.
#' @param quiet \code{TRUE} hides progress messages. Warnings and errors are
#' still shown.
#'
#' @return Invisibly, the path written.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' export_snp_list(b, file = "lead_snps.txt")
#' export_snp_list(b, dir = "lead_by_chrom")
#' }
export_snp_list <- function(blocks, file = NULL, dir = NULL, quiet = FALSE) {

  .pm_flag(quiet, "quiet")
  say <- function(...) if (!quiet) message("physmerge: ", ...)
  if (is.null(file) == is.null(dir))
    stop("give either `file` or `dir`, not both.", call. = FALSE)
  path <- if (is.null(file)) dir else file
  if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path))
    stop("an output path is empty.", call. = FALSE)
  if (!is.null(dir) && !"CHROM" %in% names(blocks))
    stop("a per-chromosome list needs a chromosome column.", call. = FALSE)
  if (nrow(blocks) == 0L)
    warning("physmerge: no blocks to export.", call. = FALSE)

  ids <- if ("rps_ID" %in% names(blocks)) as.character(blocks$rps_ID)
         else vapply(blocks$rps_BP, format, "", digits = 15, scientific = FALSE, trim = TRUE)

  if (!is.null(file)) {
    writeLines(ids, file)
    say("wrote ", length(ids), " IDs to ", file, ".")
    return(invisible(file))
  }

  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  chroms <- unique(as.character(blocks$CHROM))
  safe <- gsub("[/\\\\.]", "_", chroms)
  if (anyDuplicated(safe))
    stop("chromosome names collide once made safe for a file name: ",
         paste(unique(chroms[safe %in% safe[duplicated(safe)]]), collapse = ", "),
         call. = FALSE)
  for (i in seq_along(chroms))
    writeLines(ids[as.character(blocks$CHROM) == chroms[i]],
               file.path(dir, paste0("snp_ch", safe[i], ".txt")))
  say("wrote ", length(chroms), " chromosome file(s) to ", dir, ".")
  invisible(dir)
}
