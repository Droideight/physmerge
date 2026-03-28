#' Annotate merged blocks with full rows from the original input
#'
#' Joins the original input data back to the merged block table on
#' \code{rps_BP}, returning one row per block containing all original columns
#' for the representative SNP.
#'
#' @param blocks  Data frame returned by \code{\link{physical_merge}} (or
#'   \code{\link{run_physical_merge_from_csv}}).
#' @param data    The original input data frame passed to (or returned by)
#'   \code{\link{read_sumstat}} / \code{\link{read_gpcm_csv}}.  Must contain a
#'   \code{position} column.
#' @param chrom_col Name of the chromosome column in \code{data}.  If
#'   \code{NULL} (default), the function looks for \code{"CHROM"} then
#'   \code{"#CHROM"}.
#' @param id_col  Name of the SNP ID column in \code{data} used to populate
#'   \code{rps_ID} in the output.  If \code{NULL} (default), the function
#'   looks for \code{"ID"} then \code{"SNP"}.  Set to \code{NA} to skip.
#'
#' @return The \code{blocks} data frame with additional columns from
#'   \code{data} for the representative SNP row.  A column \code{rps_ID} is
#'   added when an ID column is found.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' df     <- read_sumstat("my_gwas.glm.linear", format = "plink2")
#' blocks <- physical_merge(df, sig_th = 5e-8, window = 500000)
#' blocks_annotated <- annotate_blocks(blocks, df)
#' }
annotate_blocks <- function(blocks, data,
                             chrom_col = NULL,
                             id_col    = NULL) {

  if (nrow(blocks) == 0L) return(blocks)
  if (!"position" %in% names(data))
    stop("`data` must contain a 'position' column.")

  # Resolve chrom column
  if (is.null(chrom_col)) {
    chrom_col <- if ("CHROM" %in% names(data)) "CHROM" else
                 if ("#CHROM" %in% names(data)) "#CHROM" else NULL
  }

  # Resolve ID column
  if (is.null(id_col)) {
    id_col <- if ("ID" %in% names(data)) "ID" else
              if ("SNP" %in% names(data)) "SNP" else NA
  }

  # Deduplicate data on position (keep first occurrence per BP)
  data_dedup <- data[!duplicated(data$position), ]

  # Build join key: all original columns for representative SNPs
  repr <- data_dedup[data_dedup$position %in% blocks$rps_BP, ]

  # Add rps_BP as join key
  repr$rps_BP <- repr$position

  # Select columns to attach: CHROM + ID (if available) + everything else
  keep <- "rps_BP"
  if (!is.null(chrom_col) && chrom_col %in% names(repr))
    keep <- c(keep, chrom_col)
  if (!is.na(id_col) && id_col %in% names(repr)) {
    repr$rps_ID <- repr[[id_col]]
    keep <- c(keep, "rps_ID")
  }

  # Attach all remaining original columns (excluding interface cols)
  extra <- setdiff(names(repr), c(keep, "position", "value"))
  keep  <- c(keep, extra)
  repr  <- repr[, keep, drop = FALSE]

  # Merge onto blocks
  out <- merge(blocks, repr, by = "rps_BP", all.x = TRUE)

  # Restore sensible column order: serial, CHROM, start, end, rps_BP, rps_ID, rps_value, rest
  priority <- intersect(c("serial", chrom_col, "start", "end",
                           "rps_BP", "rps_ID", "rps_value"), names(out))
  rest      <- setdiff(names(out), priority)
  out       <- out[, c(priority, rest), drop = FALSE]
  out[order(out$serial), ]
}


# ==============================================================================

#' Export a SNP ID list from merged blocks
#'
#' Writes the representative SNP IDs to one or more plain-text files (one ID
#' per line).  When \code{by_chrom = TRUE}, a separate file is written for
#' each chromosome and the results are bundled into a ZIP archive.
#'
#' @param blocks   Data frame returned by \code{\link{physical_merge}} or
#'   \code{\link{annotate_blocks}}.  Should contain an \code{rps_ID} column
#'   (added by \code{\link{annotate_blocks}}) or at least \code{rps_BP}.
#' @param path     Output file path.
#'   \itemize{
#'     \item When \code{by_chrom = FALSE}: path to a \code{.txt} file.
#'     \item When \code{by_chrom = TRUE}: path to a \code{.zip} archive.
#'   }
#' @param by_chrom Logical.  If \code{FALSE} (default), write a single merged
#'   file.  If \code{TRUE}, write one file per chromosome bundled in a ZIP.
#'   Requires a \code{CHROM} column in \code{blocks}.
#' @param id_col   Name of the ID column to write.  Defaults to \code{"rps_ID"}
#'   if present, otherwise falls back to \code{"rps_BP"}.
#'
#' @return Invisibly returns \code{path}.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Single file
#' export_snp_list(blocks_annotated, "snp_ids.txt")
#'
#' # Per-chromosome zip
#' export_snp_list(blocks_annotated, "snp_ids_by_chr.zip", by_chrom = TRUE)
#' }
export_snp_list <- function(blocks, path, by_chrom = FALSE, id_col = NULL) {

  if (nrow(blocks) == 0L) {
    warning("No blocks to export.")
    return(invisible(path))
  }

  # Resolve ID column
  if (is.null(id_col))
    id_col <- if ("rps_ID" %in% names(blocks)) "rps_ID" else "rps_BP"
  if (!id_col %in% names(blocks))
    stop("Column '", id_col, "' not found in blocks.")

  ids <- as.character(blocks[[id_col]])

  if (!by_chrom) {
    # ── Single merged file ────────────────────────────────────────────────────
    writeLines(ids, path)
    message("Wrote ", length(ids), " IDs to ", path)

  } else {
    # ── Per-chromosome ZIP ────────────────────────────────────────────────────
    if (!"CHROM" %in% names(blocks))
      stop("by_chrom = TRUE requires a 'CHROM' column in blocks.")

    tmp_dir <- tempfile(pattern = "physmerge_export_")
    dir.create(tmp_dir)
    on.exit(unlink(tmp_dir, recursive = TRUE))

    chroms <- sort(unique(as.character(blocks$CHROM)))
    for (ch in chroms) {
      ch_ids  <- ids[as.character(blocks$CHROM) == ch]
      outfile <- file.path(tmp_dir, paste0("snp_ch", ch, ".txt"))
      writeLines(ch_ids, outfile)
    }

    # Zip without directory paths (-j = junk paths)
    old_wd <- getwd()
    setwd(tmp_dir)
    utils::zip(path, files = list.files(tmp_dir), flags = "-j")
    setwd(old_wd)

    message("Wrote ", length(chroms), " chromosome file(s) to ", path)
  }

  invisible(path)
}
