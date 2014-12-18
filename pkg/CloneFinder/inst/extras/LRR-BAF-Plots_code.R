setwd("C:\\Users\\Mark\\Documents\\Lab\\Data\\BAF sets\\Data Sets")
dfbaf<-na.omit(read.table("BAF.txt"))
#BAF.txt is a file contatining the BAF at each locus in the SNP chip for all samples
dfseg<-na.omit(read.table("Segments.txt"))
#Segments.txt is basically the 'rewind' file but with an added column denoting segment number (i.e., the fifth segment of chrom 3 from CL001 is number '1.3.5')
dfbaf$chrom<-dfbaf$Chr

func<-function(cid){
df1<-dfbaf[,1:3]
df1$chrom<-df1$Chr
df1<-cbind(df1, dfbaf[[cid]])
df2<-dfseg[dfseg$SamID==cid,]
require(data.table)
DT1 <- data.table(df1, key = c('chrom','Position'))
DT2 <- data.table(df2, key = c('chrom','loc.start'))
# this will perform the join you want (but retain all the 
# columns with names names of DT2)
# DT2[DT1, roll=TRUE]
# which is why I have renamed and subset here) 
BAF<-dfbaf[[cid]]
segdata<-DT2[DT1, roll=TRUE][ ,list(chrom,Position = loc.start,seg)]
}

funccid<-function(cid){
dfbaf[[cid]]
}

library(reshape)

funcz<-function(cid){
d<-cbind(func(cid), funccid(cid))
rename(d, c("V2"="BAF"))
}

#Here I'm just exporting a file for each sample which contains the BAF at each locus as well as which segment each locus belongs to;
#This intermediate exporting of data only to be read back in was done because, though I have apparently functional code to do everything
#without this intermediate step, it seems to take forever to complete the iterations for the EM below; for some reason when I include two 
#intermediate steps it goes quicker.
funcw<-function(cid){
write.table(funcz(cid), paste("C:\\Users\\Mark\\Documents\\Lab\\Data\\BAF sets\\EM\\Seg_data\\seg.", cid, ".txt", sep=""))
}

targfunc <- function(param, q, z) {
    a <- param[1]
    b <- param[2]
    temp11 <- sum(z[,1]*(-lbeta(a, b) + (a-1) * log(q) + (b-1) * log(1-q)))
    temp12 <- sum(z[,2]*(-lbeta(b, a) + (b-1) * log(q) + (a-1) * log(1-q)))
    - temp11 - temp12
}

#Your EM algorithm code
EMalg <- function(datavec, forever=100, epsilon=0.001) {
    # initialize cluster assignments
    n <- length(datavec)
    z1 <- rep(0, n)
    z1[datavec < 0.5] <- 1
    z <- matrix(c(z1, 1-z1), ncol=2)
    # loop control
    lastlike <- -10^5
    currlike <- 0
    count <- 0
    while ((count < forever) & (abs(lastlike - currlike) > epsilon)) {
        count = count + 1
        ## M-step
        maxlike <- nlm(targfunc, c(1, 1), q=datavec, z=z, stepmax=10000, print.level=1)
        mle <- maxlike$estimate
        lastlike <- currlike
        currlike <- maxlike$minimum
        if (any(mle < 0)) stop("One of the mle's is negative!")
        ## E-step
        ppi <- apply(z, 2, mean)
        mix1 <- ppi[1]*dbeta(datavec, mle[1], mle[2])
        mix2 <- ppi[2]*dbeta(datavec, mle[2], mle[1])
        z[,1] <- mix1/(mix1+mix2)
        z[,2] <- mix2/(mix1+mix2)
    }
    list(ppi=ppi, mle=mle, z=z)
}

setwd("C:\\Users\\Mark\\Documents\\Lab\\Data\\BAF sets\\Data Sets")
segments<-read.table("Segments.txt")
#Notice this is the folder I exported the segmented BAF data above; now reading it back in to apply EM algorith to get 'a' and 'b' for each segment.
setwd("C:\\Users\\Mark\\Documents\\Lab\\Data\\BAF sets\\EM\\Seg_data")

cidfunc<-function(cid){
d<-read.table(file=paste("seg.", cid, ".txt", sep=""))
func<-function(seg){
r<-d[d$seg==seg,]
q<-r$BAF
q<-q[q<0.97 & q>0.03]
stuff <- EMalg(q)
ppi <- stuff$ppi
z <- stuff$z
mle <- stuff$mle
a <- mle[1]
b <- mle[2]
xx <- seq(0,1, length=length(q))
yy <- ppi[1]*dbeta(xx, a, b) + ppi[2]*dbeta(xx, b, a)
baf<-a/(a+b)
baf
}
seglist<-unique(d$seg)
v_baf<-c(sapply(seglist, func))
baf_df<-data.frame(seglist, v_baf)
#Here I create another folder with all the BAF 'a' and 'b' parameters for each segment with a seperate txt file for each sample
write.table(baf_df, file=paste("C:\\Users\\Mark\\Documents\\Lab\\Data\\BAF sets\\EM\\Seg_data_2\\Segments_BAF.", cid, ".txt", sep=""))
}

setwd("C:\\Users\\Mark\\Documents\\Lab\\Data\\BAF sets\\EM\\Seg_data_2")
rfunc<-function(cid){
data<-read.table(file=paste("Segments_BAF.", cid, ".txt", sep=""))
data
}

segment_BAFs<-rbind(lapply(cidlist, rfunc))
write.table(segment_BAFs, "C:\\Users\\Mark\\Documents\\Lab\\Data\\BAF sets\\EM\\CLL_BAFs.txt")

setwd("C:\\Users\\Mark\\Documents\\Lab\\Data\\BAF sets\\Data Sets")
load("colorset.rda")
#the plot code you sent me a while back
ourPlot <- function(cid, dataset) {
    temp <- dataset[dataset$SamID==cid,]
    if(nrow(temp) < 1) {
        warning("no data for", cid)
        return(invisible(NULL))
    }
    temp <- temp[temp$chrom != "X",]
    size <- 1 + round(sqrt(temp$num.mark)/15)/2
    filt <- size > log(128)/3
    colvec <- colorset[temp$Call]
    plot(temp$seg.median, temp$AvgBAF, type='n', main=cid,
         xlim=c(-.7, 0.7), ylim=c(0, 1),
         xlab="Log R Ratio", ylab="B Allele Frequency")
    abline(v=0)
    points(0, 0.5, cex=15, col=colorset["N02.BalHet"])
    text(0, 0.5, "2", col=colorset["N02.BalHet"])
    n <- 3:6
    points(log10(n/2), 1/n, cex=15, col=colorset["N03.UnbalHet"])
    text(log10(n/2), 1/n, n, col=colorset["N03.UnbalHet"])
    points(log10(1/2), 0, cex=15, col=colorset["N01.Homozygous"])
    text(log10(1/2), 0, 1, col=colorset["N01.Homozygous"])
    points(0, 0, cex=15, col=colorset["N02.Homozygous"])
    text(0, 0, "LOH", col=colorset["N02.Homozygous"])
    points(temp$seg.median[filt], temp$AvgBAF[filt],
           cex=size[filt], col=colvec[filt], pch=16)
    invisible(temp)
}

#Now code to produce and export plots
setwd("C:\\Users\\Mark\\Documents\\LabData\\BAF sets\\EM\\CLL_BAFs.txt")
df<-read.table("CLL_BAFs.txt")
AvgBAF<-df$AvgBAF

df<-df[df$num.mark>=200,]

pfunc<-function(cid){
ourPlot(cid, df)
}

pngfunc<-function(cid){
png(filename=paste("C:\\Users\\Mark\\Documents\\Lab\\Data\\BAF sets\\EM\\CLL_plots\\", cid, ".png"))
pngfunc(cid)
dev.off()
}

#Also, I'll include some hastily and clumbsily written code to create data sets containing the 'a' and 'b' values for each segment rather than just the mean
#a/(a+b)
setwd("C:\\Users\\Mark\\Documents\\Lab\\Data\\BAF sets\\Data Sets")
cidlist<-read.table("CID_list.txt")
setwd("C:\\Users\\Mark\\Documents\\Lab\\New Data Sets")
segments<-read.table("Data_CLL.txt")
setwd("C:\\Users\\Mark\\Documents\\Lab\\Data\\BAF sets\\EM\\Seg_data")

cidfunc<-function(cid){
d<-read.table(file=paste("seg.", cid, ".txt", sep=""))
afunc<-function(seg){
r<-d[d$seg==seg,]
q<-r$BAF
q<-q[q<0.97 & q>0.03]
stuff <- EMalg(q)
ppi <- stuff$ppi
z <- stuff$z
mle <- stuff$mle
a <- mle[1]
b <- mle[2]
xx <- seq(0,1, length=length(q))
yy <- ppi[1]*dbeta(xx, a, b) + ppi[2]*dbeta(xx, b, a)
a
}
bfunc<-function(seg){
r<-d[d$seg==seg,]
q<-r$BAF
q<-q[q<0.97 & q>0.03]
stuff <- EMalg(q)
ppi <- stuff$ppi
z <- stuff$z
mle <- stuff$mle
a <- mle[1]
b <- mle[2]
xx <- seq(0,1, length=length(q))
yy <- ppi[1]*dbeta(xx, a, b) + ppi[2]*dbeta(xx, b, a)
b
}
seglist<-unique(d$seg)
avec<-sapply(seglist, afunc)
bvec<-sapply(seglist, bfunc)
params_df<-data.frame(seglist, avec, bvec)
#Again, I make new files for each sample before aggregating everthing into one 'master file.'
write.table(params_df, file=paste("C:\\Users\\Mark\\Documents\\Lab\\Data\\BAF sets\\EM\\CLL_parameters\\", cid, ".txt", sep=""))
}

#cidlist is just a list of all the SamIDs.
setwd("C:\\Users\\Mark\\Documents\\Lab\\Data\\BAF sets\\Data Sets")
cidlist<-read.table("CID_list.txt")

setwd("C:\\Users\\Mark\\Documents\\Lab\\Data\\BAF sets\\EM\\CLL_parameters")
rfunc<-function(cid){
data<-read.table(file=paste(cid, ".txt", sep=""))
data$SamID<-cid
na.omit(data)
}

d<-rbind(lapply(cidlist, rfunc))
write.table(d, file"C:\\Users\\Mark\\Documents\\Lab\\New Data Sets\\params_CLL\\.txt")
}
#I will also upload my files 'params_CLL' and 'params_control' which was made doing this exact same method only using the corresponding Hapmap data sets.