bed <-read.table("1174T.HPV_barcode_reads.bed",
                 header = FALSE, sep="\t", stringsAsFactors=FALSE, quote="")
colnames(bed) <- c('chrom', 'start', 'end', 'seq_id', 'NA1', 'strand', 'barcode')
bed$barcode <- factor(bed$barcode)
bed$strand <- factor(bed$strand)

sort(tab <- table(bed$chrom))
stdchrom <- c(paste("chr", c(1:22, "X", "Y", "M"), sep=''), "HPV16")

odd <- bed[!(bed$chrom %in% stdchrom),]
dim(odd)

bed <- bed[(bed$chrom %in% stdchrom),]
dim(bed)
bed$chrom <- factor(bed$chrom, levels=stdchrom)
summary(bed)

starts <- tapply(bed$start, list(bed$chrom), min)
ends <- tapply(bed$end, list(bed$chrom), max)

BC <- levels(bed$barcode)



showloc <- function(bc) {
  oneset <- bed[bed$barcode == bc,]
  opar <- par(mfrow=c(4,7))
  plot(c(0,1), c(0,1), type='n', xaxt='n', yaxt='n', xlab='', ylab='')
  text(0.5, 0.67, bc)
  text(0.5, 0.33, paste("N =", nrow(oneset)))
  for (ch in stdchrom) {
    x <- oneset[oneset$chrom == ch,]
    if(nrow(x) ==0) {
      plot(0,0, type='n', main=ch)
    } else {
      plot(x$start, main=ch, ylim=c(starts[ch], ends[ch]))
    }
  }
  par(opar)
  nrow(oneset)
}
showhist <- function(bc) {
  oneset <- bed[bed$barcode == bc,]
  N <- NULL
  opar <- par(mfrow=c(4,7))
  plot(c(0,1), c(0,1), type='n', xaxt='n', yaxt='n', xlab='', ylab='')
  text(0.5, 0.67, bc)
  text(0.5, 0.33, paste("N =", nrow(oneset)))
  for (ch in stdchrom) {
    x <- oneset[oneset$chrom == ch,]
    if(nrow(x) ==0) {
      plot(0,0, type='n', main=ch)
    } else {
      counts <- hist(x$start, breaks=splits[ch], main=ch, ylim=c(0,200))$counts
      N <-  c(N, counts)
    }
  }
  par(opar)
  sum(N > 20) # MAGIC
}

bc <- "GTTGCACAGCGATG" # good
showloc(bc)
showhist(bc)

bc <- "GTATGCAGGTTCAA" # bad
showloc(bc)
showhist(bc)

bc <- "CGATACAGAATGCT" # mixed
showloc(bc)
showhist(bc)

bc <- sample(BC,1) # random
nreads <- showloc(bc)
nmols <- showhist(bc)
c(nreads, nmols, nreads/nmols)

wc <- hist(wid, breaks=seq(17.5, 219.5))$counts
