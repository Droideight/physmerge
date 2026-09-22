PKG <- Sys.getenv("PHYSMERGE_PKG", unset = normalizePath("."))
for (f in list.files(file.path(PKG,"R"), pattern="\\.R$", full.names=TRUE)) source(f)
t <- function(id, expr) {
  ws <- character()
  r <- withCallingHandlers(
      tryCatch(list(ok=TRUE, v=suppressMessages(eval(expr, parent.frame()))),
               error=function(e) list(ok=FALSE, v=conditionMessage(e))),
      warning=function(w){ ws <<- c(ws, conditionMessage(w)); invokeRestart("muffleWarning")})
  cat(sprintf("\n[%s] %s\n", id, if (r$ok) "" else paste("ERROR:", r$v)))
  for (w in ws) cat("   WARN:", w, "\n"); if (r$ok) print(r$v); invisible(r$v)
}
cat("================ M3 annotate_blocks ================\n")

# N1 normal single chrom
dat <- data.frame(CHROM=1, POS=c(1000,1200,900000), ID=c("rsA","rsB","rsC"),
                  A1=c("A","C","G"), position=c(1000,1200,900000), value=c(1e-9,1e-10,1e-12))
b <- physical_merge(dat, 5e-8, 5e5)
t("N1 normal", quote(annotate_blocks(b, dat)))

# E1 multi-allelic: two variants at the same bp, the 2nd is the lead
dat2 <- data.frame(CHROM=1, POS=c(1000,1000), ID=c("rs_ref","rs_lead"),
                   position=c(1000,1000), value=c(0.4,1e-12))
b2 <- physical_merge(dat2, 5e-8, 5e5)
t("E1 multi-allelic lead", quote(annotate_blocks(b2, dat2)))

# E2 fallback path: blocks built by hand (no rps_row attribute)
b3 <- b; attr(b3,"rps_row") <- NULL
t("E2 no rps_row attribute", quote(annotate_blocks(b3, dat)))

# E3 empty blocks
t("E3 empty blocks", quote(annotate_blocks(physical_merge(data.frame(position=1,value=0.5),5e-8,500), dat)))

# E4 no ID column
dat4 <- dat[, c("CHROM","POS","position","value")]
t("E4 no ID column", quote(annotate_blocks(physical_merge(dat4,5e-8,5e5), dat4)))

# E5 ID column called SNP
dat5 <- dat; names(dat5)[names(dat5)=="ID"] <- "SNP"
t("E5 ID column named SNP", quote(annotate_blocks(physical_merge(dat5,5e-8,5e5), dat5)))

# E6 multi-chrom
dat6 <- data.frame(CHROM=c(1,1,2,2), POS=c(1000,2000,1000,2000), ID=c("a","b","c","d"),
                   position=c(1000,2000,1000,2000), value=c(1e-9,0.5,1e-9,1e-12))
t("E6 multi-chrom", quote(annotate_blocks(physical_merge(dat6,5e-8,500), dat6)))

# E7 mismatched data (blocks from one data set, annotated with another)
t("E7 wrong data frame", quote(annotate_blocks(b, dat6)))

# E8 #CHROM not normalised
dat8 <- dat; names(dat8)[1] <- "#CHROM"
t("E8 '#CHROM' column name", quote(annotate_blocks(physical_merge(dat8,5e-8,5e5,chrom_col="#CHROM"), dat8)))

# E9 data has a column literally called rps_BP already
dat9 <- dat; dat9$rps_BP <- 999
t("E9 data already has rps_BP", quote(annotate_blocks(physical_merge(dat9,5e-8,5e5), dat9)))

# E10 keep_* toggles all FALSE
t("E10 all keeps FALSE", quote(annotate_blocks(b, dat, keep_serial=FALSE, keep_start=FALSE,
     keep_end=FALSE, keep_rps_BP=FALSE, keep_rps_value=FALSE, keep_rps_ID=FALSE)))

cat("\n================ M4 export_snp_list ================\n")
tmp <- tempfile("exp"); dir.create(tmp)
ab <- annotate_blocks(physical_merge(dat6,5e-8,500), dat6)
t("N1 plain list", quote({p<-file.path(tmp,"snp.txt"); export_snp_list(ab,p); readLines(p)}))
t("E1 by_chrom zip", quote({p<-file.path(tmp,"snp.zip"); export_snp_list(ab,p,by_chrom=TRUE); c(p, if(file.exists(p)) unzip(p,list=TRUE)$Name else "NO FILE")}))
t("E2 relative path by_chrom", quote({owd<-setwd(tmp); on.exit(setwd(owd)); export_snp_list(ab,"rel.zip",by_chrom=TRUE); file.exists(file.path(tmp,"rel.zip"))}))
t("E3 empty blocks", quote(export_snp_list(ab[0,], file.path(tmp,"e.txt"))))
t("E4 no CHROM but by_chrom", quote(export_snp_list(annotate_blocks(physical_merge(dat,5e-8,5e5),dat), file.path(tmp,"x.zip"), by_chrom=TRUE)))
t("E5 no rps_ID -> falls back to rps_BP", quote({b<-physical_merge(dat,5e-8,5e5); p<-file.path(tmp,"bp.txt"); export_snp_list(b,p); readLines(p)}))
t("E6 chrom name with a slash", quote({ab2<-ab; ab2$CHROM<-c("1/2","1/2"); p<-file.path(tmp,"sl.zip"); export_snp_list(ab2,p,by_chrom=TRUE); file.exists(p)}))
t("E7 unwritable path", quote(export_snp_list(ab, "/nope/dir/snp.txt")))
t("E8 IDs containing NA", quote({ab3<-ab; ab3$rps_ID<-c(NA,"x"); p<-file.path(tmp,"na.txt"); export_snp_list(ab3,p); readLines(p)}))
