idiogram <- read.table("cytoband.txt", sep="\t") # chr, start, end, band, gstain
idiogram$V1 <- factor(idiogram$V1, levels=paste("chr", c(1:22, "X", "Y"), sep=''))
chlens <- as.vector(tapply(idiogram$V3, list(idiogram$V1), max))

save(chlens, file = "../../Data/chlens.rda")
