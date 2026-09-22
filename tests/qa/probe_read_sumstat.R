# Probe suite: R modules, normal + extreme cases.
PKG <- Sys.getenv("PHYSMERGE_PKG", unset = normalizePath("."))
for (f in list.files(file.path(PKG,"R"), pattern="\\.R$", full.names=TRUE)) source(f)

tmp <- tempfile("qa"); dir.create(tmp)
res <- list()
t <- function(id, expr) {
  r <- tryCatch(list(ok=TRUE, v=suppressMessages(suppressWarnings(eval(expr)))),
                error=function(e) list(ok=FALSE, v=conditionMessage(e)))
  w <- tryCatch({withCallingHandlers(suppressMessages(eval(expr)), warning=function(w) {
        res[[id]] <<- c(res[[id]], paste("WARN:", conditionMessage(w))); invokeRestart("muffleWarning")}); NULL},
      error=function(e) NULL)
  cat(sprintf("[%-28s] %s\n", id, if (r$ok) "ok" else paste("ERROR:", r$v)))
  invisible(r)
}
show <- function(x) { print(x); cat("\n") }

cat("\n================ M1 read_sumstat ================\n")

# N1 plink2 normal (with TEST column)
p1 <- file.path(tmp,"n1.glm.linear")
writeLines(c("#CHROM\tPOS\tID\tTEST\tP",
             "1\t1000\trsA\tADD\t1e-9",
             "1\t1000\trsA\tSEX\t0.5",
             "1\t900000\trsB\tADD\t2e-10",
             "1\t900100\trsC\tADD\t0.3"), p1)
r <- t("N1 plink2 normal", quote(read_sumstat(p1, "plink2"))); show(r$v$data[,c("CHROM","POS","ID","P","position","value")]); cat("reward:", r$v$reward, "\n")

# E1 custom with no id_col supplied
p2 <- file.path(tmp,"e1.txt")
writeLines(c("CHR BP SNP PV","1 100 s1 1e-9","1 200 s2 1e-9"), p2)
t("E1 custom, id_col not given", quote(read_sumstat(p2, "custom", chrom_col="CHR", pos_col="BP", value_col="PV")))
t("E1b custom, id_col=NA", quote(read_sumstat(p2, "custom", chrom_col="CHR", pos_col="BP", value_col="PV", id_col=NA)))
t("E1c custom, id_col given", quote(read_sumstat(p2, "custom", chrom_col="CHR", pos_col="BP", value_col="PV", id_col="SNP")))

# E2 multi-chromosome: is the sort chromosome-aware?
p3 <- file.path(tmp,"e2.tsv")
writeLines(c("#CHROM\tPOS\tID\tP","1\t500\ta\t1e-9","2\t100\tb\t1e-9","1\t900\tc\t1e-9","2\t700\td\t1e-9"), p3)
r <- t("E2 multi-chrom row order", quote(read_sumstat(p3,"gpcm", value_col="P")))
show(r$v$data[,c("CHROM","POS","ID")])

# E3 NA / '.' / NA-string values
p4 <- file.path(tmp,"e3.tsv")
writeLines(c("#CHROM\tPOS\tID\tP","1\t100\ta\t1e-9","1\tNA\tb\t1e-9","1\t300\tc\tNA","1\t400\td\t.","1\t500\te\t1e-12"), p4)
r <- t("E3 NA handling", quote(read_sumstat(p4,"gpcm", value_col="P"))); show(r$v$data[,c("POS","ID","value")])

# E4 TEST filter removes everything
p5 <- file.path(tmp,"e4.tsv")
writeLines(c("#CHROM\tPOS\tID\tTEST\tP","1\t100\ta\tSEX\t1e-9"), p5)
t("E4 TEST filter empties", quote(read_sumstat(p5,"plink2")))

# E5 chrom filter empties
t("E5 chrom filter empties", quote(read_sumstat(p3,"gpcm", value_col="P", chrom=22)))

# E6 gz
p6 <- file.path(tmp,"e6.tsv.gz"); con <- gzfile(p6,"w")
writeLines(c("#CHROM\tPOS\tID\tP","1\t100\ta\t1e-9"), con); close(con)
t("E6 gzip input", quote(read_sumstat(p6,"gpcm", value_col="P")))

# E7 scientific / huge / negative positions
p7 <- file.path(tmp,"e7.tsv")
writeLines(c("#CHROM\tPOS\tID\tP","1\t-50\ta\t1e-9","1\t1e6\tb\t1e-9","1\t3000000000\tc\t1e-9"), p7)
r <- t("E7 odd positions", quote(read_sumstat(p7,"gpcm", value_col="P"))); show(r$v$data[,c("POS","position")])

# E8 p-values of 0 and 1 and negative
p8 <- file.path(tmp,"e8.tsv")
writeLines(c("#CHROM\tPOS\tID\tP","1\t100\ta\t0","1\t200\tb\t1","1\t300\tc\t-1"), p8)
r <- t("E8 p = 0 / 1 / -1", quote(read_sumstat(p8,"gpcm", value_col="P"))); show(r$v$data[,c("POS","value")])

# E9 duplicate column names
p9 <- file.path(tmp,"e9.tsv")
writeLines(c("#CHROM\tPOS\tID\tP\tP","1\t100\ta\t1e-9\t0.5"), p9)
t("E9 duplicate column names", quote(read_sumstat(p9,"gpcm", value_col="P")))

# E10 file where value column is already named 'value' / 'position'
p10 <- file.path(tmp,"e10.tsv")
writeLines(c("#CHROM\tposition\tID\tvalue","1\t100\ta\t1e-9","1\t200\tb\t1e-9"), p10)
r <- t("E10 cols named position/value", quote(read_sumstat(p10,"custom", chrom_col="#CHROM", pos_col="position", id_col="ID", value_col="value"))); show(r$v$data)

# E11 LOG10_P
p11 <- file.path(tmp,"e11.tsv")
writeLines(c("#CHROM\tPOS\tID\tLOG10_P","1\t100\ta\t9.3"), p11)
r <- t("E11 LOG10_P reward", quote(read_sumstat(p11,"custom", chrom_col="#CHROM", pos_col="POS", id_col="ID", value_col="LOG10_P"))); cat("reward:", r$v$reward, "\n")

# E12 empty file / header only
p12 <- file.path(tmp,"e12.tsv"); writeLines("#CHROM\tPOS\tID\tP", p12)
t("E12 header-only file", quote(read_sumstat(p12,"gpcm", value_col="P")))
p13 <- file.path(tmp,"e13.tsv"); file.create(p13)
t("E13 completely empty file", quote(read_sumstat(p13,"gpcm", value_col="P")))
t("E14 nonexistent file", quote(read_sumstat(file.path(tmp,"nope.tsv"),"gpcm", value_col="P")))
