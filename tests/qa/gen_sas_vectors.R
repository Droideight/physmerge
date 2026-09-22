# Build the reference test vectors for the SAS self-test, using the R package
# as the source of truth.
PKG <- Sys.getenv("PHYSMERGE_PKG", unset = normalizePath("."))
for (f in list.files(file.path(PKG,"R"), pattern="\\.R$", full.names=TRUE)) source(f)

set.seed(11)
cases <- list()
add <- function(name, chrom, pos, val, id, sig, win, reward="min", reset="best")
  cases[[length(cases)+1]] <<- list(name=name, chrom=as.character(chrom), pos=pos, val=val,
                                    id=id, sig=sig, win=win, reward=reward, reset=reset)

add("t01_two_peaks", rep(1,4), c(1000,1200,900000,901000), c(1e-9,1e-10,1e-12,0.4),
    paste0("rs",1:4), 5e-8, 5e5)
add("t02_none_sig",  rep(1,3), c(1,2,3), c(.9,.8,.7), paste0("rs",1:3), 5e-8, 5e5)
add("t03_single",    1, 1000, 1e-9, "rs1", 5e-8, 5e5)
add("t04_dup_pos",   rep(1,3), c(1000,1000,1000), c(1e-9,1e-12,1e-10), paste0("rs",1:3), 5e-8, 5e5)
add("t05_tiny_win",  rep(1,3), c(1000,1001,1002), rep(1e-9,3), paste0("rs",1:3), 5e-8, 1)
add("t06_chain_best",rep(1,11), seq(1000,5000,by=400), rep(1e-9,11), paste0("rs",1:11), 5e-8, 1000)
add("t07_chain_any", rep(1,11), seq(1000,5000,by=400), rep(1e-9,11), paste0("rs",1:11), 5e-8, 1000,
    reset="any")
add("t08_multichrom",c(1,1,2,2,"X"), c(1000,2000,1000,2000,5000), c(1e-9,0.5,1e-9,1e-12,1e-20),
    paste0("rs",1:5), 5e-8, 500)
add("t09_reward_max",rep(1,3), c(1000,1200,900000), c(30,45,60), paste0("rs",1:3), 5.45, 5e5, reward="max")
add("t10_collapse",  rep(1,4), c(1000,1400,1800,2200), c(1e-9,0.5,1e-9,1e-9), paste0("rs",1:4), 5e-8, 500)
add("t11_thresh_eq", rep(1,2), c(1000,2000), c(5e-8,4.9e-8), c("rs1","rs2"), 5e-8, 500)
add("t12_improve_ext",rep(1,4), c(1000,1400,1800,2200), c(1e-9,1e-10,1e-11,1e-12), paste0("rs",1:4), 5e-8, 500)
# a bigger random one
n <- 120
add("t13_random", sample(c("1","2"),n,TRUE), sample(1:40000,n,TRUE),
    10^-runif(n,0,12), paste0("rs",1:n), 1e-5, 2000)

inp <- do.call(rbind, lapply(cases, function(k)
  data.frame(case=k$name, CHROM=k$chrom, POS=sprintf("%.0f", k$pos), ID=k$id, P=k$val, stringsAsFactors=FALSE)))
par <- do.call(rbind, lapply(cases, function(k)
  data.frame(case=k$name, sig_th=format(k$sig, scientific=TRUE, digits=17, trim=TRUE), window=sprintf("%.0f", k$win), reward=k$reward, reset_on=k$reset,
             stringsAsFactors=FALSE)))

exp <- do.call(rbind, lapply(cases, function(k) {
  d <- data.frame(CHROM=k$chrom, POS=k$pos, ID=k$id, position=as.numeric(k$pos), value=k$val,
                  stringsAsFactors=FALSE)
  d <- d[order(match(d$CHROM, unique(d$CHROM)), d$position), ]
  b <- physical_merge(d, k$sig, k$win, reward=k$reward, reset_on=k$reset)
  if (nrow(b)==0L) return(NULL)
  a <- annotate_blocks(b, d)
  if (!"CHROM" %in% names(a)) a$CHROM <- d$CHROM[1]
  data.frame(case=k$name, serial=a$serial, CHROM=as.character(a$CHROM),
             start=a$start, end=a$end, rps_BP=a$rps_BP, rps_ID=a$rps_ID,
             rps_value=a$rps_value, stringsAsFactors=FALSE)
}))

out <- file.path(PKG, "sas", "testdata")
dir.create(out, recursive=TRUE, showWarnings=FALSE)
w <- function(x, f) write.table(x, file.path(out,f), sep="\t", quote=FALSE, row.names=FALSE, na="")
inp$P <- format(inp$P, scientific=TRUE, digits=17, trim=TRUE)
exp$rps_value <- format(exp$rps_value, scientific=TRUE, digits=17, trim=TRUE)
exp$start <- sprintf("%.0f", exp$start); exp$end <- sprintf("%.0f", exp$end)
exp$rps_BP <- sprintf("%.0f", exp$rps_BP)
w(inp, "vectors_input.tsv"); w(par, "vectors_params.tsv"); w(exp, "vectors_expected.tsv")
cat("cases:", nrow(par), " input rows:", nrow(inp), " expected blocks:", nrow(exp), "\n")
print(head(exp, 15))
