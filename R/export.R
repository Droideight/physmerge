#' Annotate merged blocks with full rows from the original input
#'
#' Returns one row per block containing all original columns for the
#' representative SNP.  When \code{blocks} still carries the \code{"rps_row"}
#' attribute set by \code{\link{physical_merge}}, the representative row is
#' taken directly by index, so blocks at multi-allelic sites are labelled with
#' the variant that actually led the block.  Otherwise the function falls back
#' to joining on \code{CHROM} and \code{rps_BP}, which can only keep the first
#' row at a duplicated position.  Block metadata columns (\code{serial},
#' \code{start}, \code{end}, \code{rps_BP}, \code{rps_value}) can be
#' individually included or dropped.
#'
#' @param blocks    Data frame returned by \code{\link{physical_merge}}.
#' @param data      The original input data frame passed to
#'   \code{\link{read_sumstat$data}}.  Must contain a \code{position} column.
#' @param chrom_col Name of the chromosome column in \code{data}.  If
#'   \code{NULL} (default), auto-detects \code{"CHROM"} then \code{"#CHROM"}.
#' @param id_col    Name of the SNP ID column in \code{data} used to populate
#'   \code{rps_ID}.  If \code{NULL} (default), auto-detects \code{"ID"} then
#'   \code{"SNP"}.  Set to \code{NA} to skip.
#' @param keep_serial    Logical. Include \code{serial} column. Default \code{TRUE}.
#' @param keep_start     Logical. Include \code{start} column. Default \code{TRUE}.
#' @param keep_end       Logical. Include \code{end} column. Default \code{TRUE}.
#' @param keep_rps_BP    Logical. Include \code{rps_BP} column. Default \code{TRUE}.
#' @param keep_rps_value Logical. Include \code{rps_value} column. Default \code{TRUE}.
#' @param keep_rps_ID    Logical. Include \code{rps_ID} column (when available).
#'   Default \code{TRUE}.
#'
#' @return The \code{blocks} data frame merged with full original columns for
#'   each representative SNP row, with block metadata columns selectively
#'   retained based on \code{keep_*} arguments.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' df     <- read_sumstat("my_gwas.glm.linear", format = "plink2")
#' blocks <- physical_merge(df$data, sig_th = 5e-8, window = 500000)
#'
#' # Keep all block metadata (default)
#' annotate_blocks(blocks, df$data)
#'
#' # Drop start/end/rps_value — only keep serial, rps_BP, rps_ID + original cols
#' annotate_blocks(blocks, df$data,
#'                 keep_start = FALSE, keep_end = FALSE, keep_rps_value = FALSE)
#' }
annotate_blocks <- function(blocks, data,
                            chrom_col       = NULL,
                            id_col          = NULL,
                            keep_serial     = TRUE,
                            keep_start      = TRUE,
                            keep_end        = TRUE,
                            keep_rps_BP     = TRUE,
                            keep_rps_value  = TRUE,
                            # note, although function to toggle false/true is
                            # given, if annotate_blocks is not use terminally,
                            # this should always be set to true as rps_BP is
                            # used when mergining external information.
                            keep_rps_ID     = TRUE) {
  
  if (nrow(blocks) == 0L) return(blocks)
  if (!"position" %in% names(data))
    stop("`data` must contain a 'position' column.")
  
  # Resolve chrom column
  if (is.null(chrom_col))
    chrom_col <- if ("CHROM" %in% names(data)) "CHROM" else
      if ("#CHROM" %in% names(data)) "#CHROM" else NULL
  # reading #CHROM is provided in case user did not use read_sumstat
  
  # Resolve ID column
  if (is.null(id_col))
    id_col <- if ("ID" %in% names(data)) "ID" else
      if ("SNP" %in% names(data)) "SNP" else NA
  
  # ── Exact representative rows, when physical_merge() recorded them ──────────
  # physical_merge() stores the input row index of every representative SNP in
  # attr(blocks, "rps_row").  Using it is the only way to get the right row when
  # several variants share one base-pair position (multi-allelic sites): the
  # position-based lookup below can only keep the first row at that position,
  # which may be a non-significant variant.
  rr <- attr(blocks, "rps_row")
  use_rr <- !is.null(rr) && length(rr) == nrow(blocks) && !anyNA(rr) &&
    all(rr >= 1L) && all(rr <= nrow(data)) &&
    isTRUE(all.equal(as.numeric(data$position[rr]),
                     as.numeric(blocks$rps_BP), tolerance = 0))

  if (use_rr) {
    repr   <- data[rr, , drop = FALSE]
    repr$rps_BP <- repr$position
    has_id <- !is.na(id_col) && id_col %in% names(repr)
    if (has_id) repr$rps_ID <- repr[[id_col]]
    lead  <- c(if (has_id) "rps_ID",
               if (!is.null(chrom_col) && chrom_col %in% names(repr)) chrom_col)
    extra <- setdiff(names(repr), c(lead, "rps_BP", "position", "value",
                                    if (!is.na(id_col)) id_col, chrom_col))
    repr  <- repr[, c(lead, extra), drop = FALSE]
    if (!is.null(chrom_col) && chrom_col != "CHROM" &&
        chrom_col %in% names(repr) && "CHROM" %in% names(blocks))
      names(repr)[names(repr) == chrom_col] <- "CHROM"
    add  <- repr[, setdiff(names(repr), names(blocks)), drop = FALSE]
    rownames(add) <- NULL
    out  <- cbind(blocks, add)
    attr(out, "rps_row") <- NULL
    return(.finish_annotation(out, chrom_col, keep_serial, keep_start, keep_end,
                              keep_rps_BP, keep_rps_value, keep_rps_ID))
  }

  # ── Fallback: position lookup (blocks not produced by this physical_merge) ──
  if (!is.null(chrom_col) && chrom_col %in% names(data)) {
    data_dedup <- data[!duplicated(data[, c(chrom_col, "position")]), ]
  } else {
    data_dedup <- data[!duplicated(data$position), ]
  }
  # here, when multiple SNPs share one bp, only the first instance is kept
  if ("CHROM" %in% names(blocks) && !is.null(chrom_col) && chrom_col %in% names(data_dedup)) {
    keys_block <- paste(blocks$CHROM, blocks$rps_BP, sep = ":")
    keys_data  <- paste(data_dedup[[chrom_col]], data_dedup$position, sep = ":")
    repr       <- data_dedup[keys_data %in% keys_block, ]
  } else {
    repr <- data_dedup[data_dedup$position %in% blocks$rps_BP, ]
  }
  # only keep rows that have representative SNPs
  repr$rps_BP     <- repr$position
  
  # Add rps_ID if available
  has_id <- !is.na(id_col) && id_col %in% names(repr)
  if (has_id) repr$rps_ID <- repr[[id_col]]
  
  # Columns to bring in from original data
  join_cols <- c("rps_BP",
                 if (!is.null(chrom_col) && chrom_col %in% names(repr)) chrom_col,
                 if (has_id) "rps_ID")
  extra     <- setdiff(names(repr), c(join_cols, "position", "value", id_col))
  repr      <- repr[, c(join_cols, extra), drop = FALSE]
  
  # Merge
  merge_keys <- "rps_BP"
  if ("CHROM" %in% names(blocks) && !is.null(chrom_col) && chrom_col %in% names(repr)) {
    if (chrom_col != "CHROM") names(repr)[names(repr) == chrom_col] <- "CHROM"
    merge_keys <- c("CHROM", "rps_BP")
  }
  out <- merge(blocks, repr, by = merge_keys, all.x = TRUE)
  .finish_annotation(out, chrom_col, keep_serial, keep_start, keep_end,
                     keep_rps_BP, keep_rps_value, keep_rps_ID)
}


# Internal: shared column selection / ordering for annotate_blocks()
.finish_annotation <- function(out, chrom_col, keep_serial, keep_start, keep_end,
                               keep_rps_BP, keep_rps_value, keep_rps_ID) {
  meta <- c(
    if (keep_serial)    "serial",
    if (!is.null(chrom_col) && chrom_col %in% names(out)) chrom_col,
    if ("CHROM" %in% names(out) && !identical(chrom_col, "CHROM")) "CHROM",
    if (keep_start)     "start",
    if (keep_end)       "end",
    if (keep_rps_BP)    "rps_BP",
    if (keep_rps_ID && "rps_ID" %in% names(out)) "rps_ID",
    if (keep_rps_value) "rps_value"
  )
  meta <- unique(meta)
  rest <- setdiff(names(out), c(meta, "serial", "start", "end",
                                "rps_BP", "rps_ID", "rps_value"))
  out  <- out[, c(meta, rest), drop = FALSE]
  rownames(out) <- NULL
  out[order(out$serial), ]
}


# ==============================================================================

#' Export a SNP ID list from merged blocks
#'
#' Writes the representative SNP IDs to one or more plain-text files (one ID
#' per line).  When \code{by_chrom = TRUE}, a separate file is written for
#' each chromosome and bundled into a ZIP archive.
#'
#' @param blocks   Data frame returned by \code{\link{annotate_blocks}}.
#'   Should contain \code{rps_ID} (added by \code{annotate_blocks}) or at
#'   least \code{rps_BP}.
#' @param path     Output file path.
#'   \itemize{
#'     \item \code{by_chrom = FALSE}: path to a \code{.txt} file.
#'     \item \code{by_chrom = TRUE}: path to a \code{.zip} archive.
#'   }
#' @param by_chrom Logical.  \code{FALSE} (default) writes a single merged
#'   file; \code{TRUE} writes one file per chromosome bundled in a ZIP.
#'   Requires a \code{CHROM} column in \code{blocks}.
#' @param id_col   Name of the ID column to write.  Defaults to \code{"rps_ID"}
#'   if present, otherwise \code{"rps_BP"}.
#'
#' @return Invisibly returns \code{path}.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' export_snp_list(blocks, "snp_ids.txt")
#' export_snp_list(blocks, "snp_ids_by_chr.zip", by_chrom = TRUE)
#' }
export_snp_list <- function(blocks, path, by_chrom = FALSE, id_col = NULL) {
  
  if (nrow(blocks) == 0L) {
    warning("No blocks to export.")
    return(invisible(path))
  }
  
  if (is.null(id_col))
    id_col <- if ("rps_ID" %in% names(blocks)) "rps_ID" else "rps_BP"
  if (!id_col %in% names(blocks))
    stop("Column '", id_col, "' not found in blocks.")
  
  ids <- as.character(blocks[[id_col]])
  
  if (!by_chrom) {
    writeLines(ids, path)
    message("Wrote ", length(ids), " IDs to ", path)
    
  } else {
    if (!"CHROM" %in% names(blocks))
      stop("by_chrom = TRUE requires a 'CHROM' column in blocks.")
    
    path    <- normalizePath(path, mustWork = FALSE)  # fix relative path before setwd
    tmp_dir <- tempfile(pattern = "physmerge_export_")
    dir.create(tmp_dir)
    old_wd <- getwd()
    on.exit({
      setwd(old_wd)                        # restore wd even if zip fails
      unlink(tmp_dir, recursive = TRUE)    # then clean up tmp dir
    }, add = FALSE)
    
    chroms <- sort(unique(as.character(blocks$CHROM)))
    for (ch in chroms) {
      ch_ids  <- ids[as.character(blocks$CHROM) == ch]
      writeLines(ch_ids, file.path(tmp_dir, paste0("snp_ch", ch, ".txt")))
    }
    
    setwd(tmp_dir)
    utils::zip(path, files = list.files(tmp_dir), flags = "-j")
    
    message("Wrote ", length(chroms), " chromosome file(s) to ", path)
  }
  
  invisible(path)
}