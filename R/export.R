#' Annotate merged blocks with full rows from the original input
#'
#' Joins the original input data back to the merged block table on
#' \code{rps_BP}, returning one row per block containing all original columns
#' for the representative SNP.  Block metadata columns (\code{serial},
#' \code{start}, \code{end}, \code{rps_BP}, \code{rps_value}) can be
#' individually included or dropped.
#'
#' @param blocks    Data frame returned by \code{\link{physical_merge}}.
#' @param data      The original input data frame passed to
#'   \code{\link{read_sumstat}}.  Must contain a \code{position} column.
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
                            keep_rps_ID     = TRUE) {
  
  if (nrow(blocks) == 0L) return(blocks)
  if (!"position" %in% names(data))
    stop("`data` must contain a 'position' column.")
  
  # Resolve chrom column
  if (is.null(chrom_col))
    chrom_col <- if ("CHROM" %in% names(data)) "CHROM" else
      if ("#CHROM" %in% names(data)) "#CHROM" else NULL
  
  # Resolve ID column
  if (is.null(id_col))
    id_col <- if ("ID" %in% names(data)) "ID" else
      if ("SNP" %in% names(data)) "SNP" else NA
  
  # Deduplicate data on position
  data_dedup      <- data[!duplicated(data$position), ]
  repr            <- data_dedup[data_dedup$position %in% blocks$rps_BP, ]
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
  out <- merge(blocks, repr, by = "rps_BP", all.x = TRUE)
  
  # Build ordered column list based on keep_* flags
  meta <- c(
    if (keep_serial)    "serial",
    if (!is.null(chrom_col) && chrom_col %in% names(out)) chrom_col,
    if (keep_start)     "start",
    if (keep_end)       "end",
    if (keep_rps_BP)    "rps_BP",
    if (keep_rps_ID && "rps_ID" %in% names(out)) "rps_ID",
    if (keep_rps_value) "rps_value"
  )
  rest <- setdiff(names(out), c(meta,
                                "serial", "start", "end",
                                "rps_BP", "rps_ID", "rps_value"))
  out  <- out[, c(meta, rest), drop = FALSE]
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
    
    tmp_dir <- tempfile(pattern = "physmerge_export_")
    dir.create(tmp_dir)
    on.exit(unlink(tmp_dir, recursive = TRUE))
    
    chroms <- sort(unique(as.character(blocks$CHROM)))
    for (ch in chroms) {
      ch_ids  <- ids[as.character(blocks$CHROM) == ch]
      writeLines(ch_ids, file.path(tmp_dir, paste0("snp_ch", ch, ".txt")))
    }
    
    old_wd <- getwd()
    setwd(tmp_dir)
    utils::zip(path, files = list.files(tmp_dir), flags = "-j")
    setwd(old_wd)
    
    message("Wrote ", length(chroms), " chromosome file(s) to ", path)
  }
  
  invisible(path)
}
