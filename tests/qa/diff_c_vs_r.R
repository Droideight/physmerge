# Randomised differential test: C CLI vs R implementation.
PKG <- Sys.getenv("PHYSMERGE_PKG", unset = normalizePath("."))
for (f in list.files(file.path(PKG,"R"), pattern="\\.R$", full.names=TRUE)) source(f)
PM <- file.path(PKG,"cli","physmerge")
tmp <- tempfile("dt"); dir.create(tmp)
set.seed(20260922)

run_c <- function(file, sig, win, reward, reset, extra="") {
  cmd <- sprintf("%s -i %s -f custom --chrom-col '#CHROM' --pos-col POS --id-col ID --value-col P -s %g -w %g -r %s --reset-on %s --sort -q %s",
                 shQuote(PM), shQuote(file), sig, win, reward, reset, extra)
  out <- system(cmd, intern=TRUE, ignore.stderr=TRUE)
  if (length(out) <= 1) return(data.frame())
  read.delim(text=paste(out, collapse="\n"), colClasses="character")
}
run_r <- function(df, sig, win, reward, reset) {
  d <- df; d$position <- as.numeric(d$POS); d$value <- as.numeric(d$P); d$CHROM <- as.character(d$`X.CHROM`)
  b <- physical_merge(d, sig, win, reward=reward, reset_on=reset)
  if (nrow(b)==0L) return(data.frame())
  a <- annotate_blocks(b, d)
  a
}
mism <- 0; n <- 0
for (rep in 1:400) {
  nchr <- sample(1:3,1); nsnp <- sample(1:60,1)
  win  <- sample(c(1,10,100,1000,5e4,5e5),1)
  sig  <- sample(c(5e-8,1e-5,0.05,0.5),1)
  reward <- sample(c("min","max"),1); reset <- sample(c("best","any"),1)
  ch  <- sample(as.character(1:nchr), nsnp, TRUE)
  pos <- sample(seq_len(max(10, nsnp*200)), nsnp, TRUE)
  val <- if (reward=="min") 10^-runif(nsnp,0,12) else runif(nsnp,0,60)
  if (reward=="max") sig <- sample(c(5,30,55),1)
  df <- data.frame(`#CHROM`=ch, POS=pos, ID=paste0("rs",seq_len(nsnp)), P=format(val, scientific=TRUE, digits=15), check.names=FALSE)
  f <- file.path(tmp,"in.tsv"); write.table(df, f, sep="\t", quote=FALSE, row.names=FALSE)
  cres <- run_c(f, sig, win, reward, reset)
  names(df)[1] <- "X.CHROM"
  rres <- run_r(df, sig, win, reward, reset)
  n <- n+1
  ck <- if (nrow(cres)) paste(cres$CHROM, cres$start, cres$end, cres$rps_BP, cres$rps_ID) else character(0)
  rk <- if (nrow(rres)) paste(rres$CHROM, rres$start, rres$end, rres$rps_BP, rres$rps_ID) else character(0)
  if (!identical(sort(ck), sort(rk))) {
    mism <- mism + 1
    if (mism <= 3) {
      cat("\n=== MISMATCH", mism, " win=",win," sig=",sig," reward=",reward," reset=",reset,"\n")
      cat("--- C ---\n"); print(cres); cat("--- R ---\n"); print(rres[,c("serial","CHROM","start","end","rps_BP","rps_ID","rps_value")])
      cat("--- input ---\n"); print(df)
    }
  }
}
cat(sprintf("\nDifferential test: %d cases, %d mismatches\n", n, mism))
