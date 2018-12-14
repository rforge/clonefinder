library(survival)
#library(prodlim)
library(quantmod)
library(ClassComparison)
cll.path <- 'C:/Users/Mark/OneDrive - The Ohio State University/clonetools/dat/data-CLL'
pos <- get(load('C:/Users/Mark/OneDrive - The Ohio State University/clonetools/pos.rda'))
cll.files <- list.files(cll.path)
samples <- sapply(1:length(cll.files),function(i){strsplit(strsplit(cll.files[i],split='dat-')[[1]][2],split='.rda')[[1]][1]})
cf.results <- lapply(1:length(cll.files), function(j){get(load(paste(cll.path, '/', cll.files[j], sep='')))$res})
cll.ids <- sapply(1:length(cll.files),function(j){strsplit(strsplit(cll.files[j],split='-')[[1]][4],split='.rda')[[1]][1]})
pars <- cf.results[[1]]$pars
rewind <- get(load('C:/Users/Mark/OneDrive - The Ohio State University/clonetools/SNP_array_data/rewind_both.rda'))
total.markers <- sum(rewind[rewind$SamID==samples[1],]$num.mark)
clinical <- get(load('C:/Users/Mark/Documents/Lab/clinok.rda'))
#clinical <- get(load('D:/CLL/ClinicalOK.rda'))
eligible <- get(load('C:/Users/Mark/OneDrive - The Ohio State University/clinical/eligible.rda'))
clinical$ttp.eligible <- eligible
cids <- rownames(clinical)
common.cids <- intersect(c(cids),samples)
clinical <- clinical[sapply(1:length(common.cids),function(j){which(rownames(clinical)==common.cids[j])}),]
cf.results <- lapply(1:length(common.cids),function(j){cf.results[[which(samples==common.cids[j])]]})
ttp.eligible <- which(clinical$ttp.eligible)

###
ks <- sapply(1:length(cf.results),function(j){length(which(cf.results[[j]]$psi>0))})
columns <- colnames(clinical)

###Treatment (significant):
#vars; note: time variables are in months (365.25/12).
White.blood.count <- clinical$White.blood.count
Rai.Stage <- clinical$Rai.Stage
TimeSample2SigTreat <- clinical$TimeSample2SigTreat
TimeDiag2SigTreat <- clinical$TimeDiagnosis2SigTreat
PFS.time <- clinical$PFS.time
EFS.time <- clinical$EFS.time
mutation.status <- clinical$mutation.status
OSAD <- clinical$OSAfterDiagnosis
OSAS <- clinical$OSAfterSample
TTP <- clinical$TTP.time
TxDSS <- clinical$TxDSS.time

#status (numeric)
NumericSigTreatment <- clinical$NumericSigTreatment
PFS.numericStatus <- clinical$PFS.numericStatus
EFS.numericStatus <- clinical$EFS.numericStatus
NumericVitalStatus <- clinical$NumericVitalStatus
NumericTTPStatus <- clinical$TTP.numericStatus
NumericTxDSSStatus <- clinical$TxDSS.numericStatus

########Cox Proportional Hazards:
colors <- c('red', 'blue', 'purple', 'darkorange')

###Unified segmentation of SNP data
segmentation <- function(dat, bandwidth=10){
  Reduce(rbind,lapply(1:22,function(chrom){
    data <- dat[dat$chrom==chrom,]
    breaks <- c(data$loc.start,data$loc.end)
    #plot(density(breaks,from=0,to=max(na.omit(data$loc.end)),
    #             bw=100),xlim=c(0,max(na.omit(data$loc.end))),main='Segment Breaks',xlab='Position')
    dens <- density(breaks,from=0,to=max(na.omit(data$loc.end)), bw=bandwidth)
    y <- dens$y
    peaks <- findPeaks(y)
    breaks <- round(unique(c(min(dens$x[peaks]),dens$x[peaks])))
    data.frame(chr=chrom,start=breaks[1:(length(breaks)-1)]+1,end=breaks[2:length(breaks)])
  }))
}
segments <- segmentation(rewind)
##

segInfo <- lapply(1:length(common.cids),function(i){
  resdat <- cf.results[[i]]$data$snpdata
  resdat.filt <- cf.results[[i]]$filtered.data$cndata.filt
  A.filt <- cf.results[[i]]$A
  B.filt <- cf.results[[i]]$B
  A <- B <- matrix(1,nrow=nrow(resdat),ncol=length(cf.results[[i]]$psi))
  for(j in 1:length(resdat.filt$seg)){
    A[resdat.filt$seg[j],] <- A.filt[j,]
    B[resdat.filt$seg[j],] <- B.filt[j,]
  }
  psi <- cf.results[[i]]$psi
  dat <- rewind[rewind$SamID==common.cids[i],]
  combined <- Reduce(rbind,lapply(1:22,function(j){
    resdat.chr <- resdat[resdat$chr==j,]
    A.chr <- A[resdat$chr==j,]
    B.chr <- B[resdat$chr==j,]
    CN.chr <- A.chr + B.chr
    if(length(which(resdat$chr==j))==1){
      CN.chr <- t(as.matrix(CN.chr))
    }
    dat.chr <- dat[dat$chrom==j,]
    df <- segments[segments$chr==j,]
    markers.chr <- pos[pos$Chr==j,]
    chrdf <- data.frame(t(sapply(1:nrow(df),function(k){
      overlapping <- dat.chr[which(dat.chr$loc.start<=df$end[k] & dat.chr$loc.end>=df$start[k]),]
      indices <- which(resdat.chr$LRR %in% overlapping$seg.median & resdat.chr$BAF %in% overlapping$AvgBAF)
      indices <- intersect(indices, which(!is.na(resdat.chr$BAF) & which(!is.na(resdat.chr$LRR))))
      if(length(indices)>0){
        index <- indices[which.max(resdat.chr$markers[indices])]
        up <- which(CN.chr[index,]>2)
        down <- which(CN.chr[index,]<2)
        AF <- max(sum(psi[up]), sum(psi[down]))
        if(sum(psi[down]) > sum(psi[up])){
          AF <- -1*AF
        }
        TCN <- resdat.chr$X[index] + resdat.chr$Y[index]
        markers <- length(which(markers.chr$Position <= df$end[k] & markers.chr$Position >= df$start[k]))
        row <- data.frame('chr'=j,'start'=df$start[k],'end'=df$end[k],'AF'=AF,'TCN'=TCN,'markers'=markers)
      }else{
        row <- data.frame('chr'=j,'start'=df$start[k],'end'=df$end[k], 'AF'=NA,'TCN'=NA,'markers'=NA)
      }
      row
    })))
    colnames(chrdf) <- c('chr','start','end','AF','TCN','markers')
    chrdf$chr <- unlist(chrdf$chr)
    chrdf$start <- unlist(chrdf$start)
    chrdf$end <- unlist(chrdf$end)
    chrdf$AF <- unlist(chrdf$AF)
    chrdf$TCN <- unlist(chrdf$TCN)
    chrdf$markers <- unlist(chrdf$markers)
    chrdf
  }))
  combined
})
names(segInfo) <- common.cids
save(segInfo, file='C:/Users/Mark/OneDrive - The Ohio State University/clinical/segInfo.rda')

###Response data:
response.status <- rep(0,length(ttp.eligible))
response.status[which(clinical$Response.1st.sig.treat[ttp.eligible]=='CR')] <- 1


###
load('C:/Users/Mark/OneDrive - The Ohio State University/clinical/segInfo.rda')

#####Clinical vars as a function of allele fraction (or total CN): 
suffixes <- c('OSAS','T2ST','PFS','EFS','TTP','OSAT')
nameset <- c('Overall Survival (from sample)','Time to Treatment','Progression Free Survival',
             'Event Free Survival','Time to Progression','Overall Survival (from treatment)')
filtered <- c(FALSE,FALSE,TRUE,TRUE,TRUE,TRUE)
survtimes <- list(OSAS,TimeDiag2SigTreat,PFS.time,EFS.time,TTP,OSAS)
statuses <- list(NumericVitalStatus,NumericSigTreatment,PFS.numericStatus,EFS.numericStatus,NumericTTPStatus,NumericVitalStatus)
pmat <- pmat.bin <- pmat.cut <- matrix(NA,nrow(segInfo[[1]]),ncol=length(suffixes))
response.pval <- response.oddsRatio <- rep(NA, nrow(segInfo[[1]]))
colnames(pmat) <- colnames(pmat.bin) <- colnames(pmat.cut) <- suffixes
loessfits <- array(NA,dim=c(nrow(segInfo[[1]]),length(suffixes),2,50))

#nMin is the minimum number of non-na alteration fractions for analysis to continue:
nMin <- 20
#kMin is the minimum number of entries for each level necessary to proceed; i.e.,
#if all or almost all samples are unaltered, there's no point in trying to compute
#a p-value
kMin <- 4

#Cox PH Regression (0 vs. nonzero) and Loess fits:
for(seg in 1:nrow(segInfo[[1]])){
  afs <- unlist(sapply(1:length(segInfo),function(j){segInfo[[j]]$AF[seg]}))
  tcns <- unlist(sapply(1:length(segInfo),function(j){segInfo[[j]]$TCN[seg]}))
  afs[which(is.na(tcns))] <- NA
  af.status <- rep(NA, length(afs))
  if(length(na.omit(afs))>=nMin){
    notna <- which(!is.na(afs))
    if(mean(na.omit(afs))>0){
      af.status[notna] <- 0
      af.status[notna][which(afs[notna]>0)] <- 1
    }else{
      af.status[notna] <- 0
      af.status[notna][which(afs[notna]<0)] <- 1
    } 
  }
  L <- length(which(sapply(0:1,function(i){length(which(af.status==i))>kMin})))
  if(L==2 & length(which(!is.na(afs)))>=nMin){
    tbl <- table(af.status[ttp.eligible], response.status)
    test <- fisher.test(tbl)
    response.pval[seg] <- test$p.value
    response.oddsRatio[seg] <- test$estimate
    for(h in 1:length(suffixes)){
      segname <- paste('chr',segInfo[[1]]$chr[seg][[1]],'_',round(segInfo[[1]]$start[seg][[1]]),sep='')
      filename.reg <- paste('C:/Users/Mark/OneDrive - The Ohio State University/clinical/AF_reg/',segname,'-',suffixes[h],'.png',sep='')
      filename.loess <- paste('C:/Users/Mark/OneDrive - The Ohio State University/clinical/AF_loess/',segname,'-',suffixes[h],'.png',sep='')
      if(filtered[h]){
        af.var <- afs[ttp.eligible]
        af.status.var <- af.status[ttp.eligible]
        cph <- coxph(Surv(survtimes[[h]][ttp.eligible],statuses[[h]][ttp.eligible])~af.var)
        cph.bin <- coxph(Surv(survtimes[[h]][ttp.eligible],statuses[[h]][ttp.eligible])~af.status.var)
        fit <- survfit(Surv(survtimes[[h]][ttp.eligible],statuses[[h]][ttp.eligible])~af.status.var)
      }else{
        af.var <- afs
        af.status.var <- af.status
        cph <- coxph(Surv(survtimes[[h]],statuses[[h]])~af.var)
        cph.bin <- coxph(Surv(survtimes[[h]],statuses[[h]])~af.status.var)
        fit <- survfit(Surv(survtimes[[h]],statuses[[h]])~af.status.var) 
      }
      ns <- fit$n
      p <- round(summary(cph)$logtest[3], digits=4)
      p.bin <- round(summary(cph.bin)$logtest[3], digits=4)
      pmat[seg,h] <- p
      pmat.bin[seg,h] <- p.bin
      if(p<.05){
        mainlab <- paste(nameset[h],' - Chrom ',round(segInfo[[1]]$chr[seg][[1]]),':',round(segInfo[[1]]$start[seg][[1]]),
                         '-',round(segInfo[[1]]$end[seg][[1]]),sep='')
        martingales <- residuals(cph,type="martingale")
        loessfit <- loess.smooth(cph$linear.predictors, martingales)
        loessfits[seg,h,1,] <- loessfit$x
        loessfits[seg,h,2,] <- loessfit$y
        legendlabs <- c(paste('AF = 0; n = ',ns[1],sep=''),paste('AF > 0; n = ',ns[2],sep=''))
        xmax <- max(na.omit(survtimes[[h]]))
        png(filename.reg,height=800,width=1000)
        plot(fit,col=colors[1:length(ns)],main=paste(mainlab,'; p-value = ',p,sep=''))
        legend(xmax-.35*xmax,1,legendlabs,col=colors[1:length(ns)],lty=rep(1,3),cex=.9)
        dev.off()
        png(filename.loess,height=800,width=1000)
        scatter.smooth(na.omit(af.var), martingales,xlab='Allele Fraction', ylab='Martingale Residual',main=mainlab)
        dev.off()
      }
    }
  }
}
save(pmat,file='C:/Users/Mark/OneDrive - The Ohio State University/clinical/seg_pvals.rda')
save(pmat.bin,file='C:/Users/Mark/OneDrive - The Ohio State University/clinical/seg_pvals_bin.rda')
save(response.pval,file='C:/Users/Mark/OneDrive - The Ohio State University/clinical/response.pval.rda')
save(response.oddsRatio,file='C:/Users/Mark/OneDrive - The Ohio State University/clinical/response.oddsRatio.rda')
save(loessfits,file='C:/Users/Mark/OneDrive - The Ohio State University/clinical/loess_fits.rda')


###
load('C:/Users/Mark/OneDrive - The Ohio State University/clinical/segInfo.rda')
markers <- segInfo[[1]]$markers
pvals <- get(load('C:/Users/Mark/OneDrive - The Ohio State University/clinical/seg_pvals.rda'))
pvals.bin <- get(load('C:/Users/Mark/OneDrive - The Ohio State University/clinical/seg_pvals_bin.rda'))
response.pval <- get(load('C:/Users/Mark/OneDrive - The Ohio State University/clinical/response.pval.rda'))
response.oddsRatio <- get(load('C:/Users/Mark/OneDrive - The Ohio State University/clinical/response.oddsRatio.rda'))
thresh <- 100
hist(pvals[which(markers>thresh),1],breaks=30,xlab='P-Value',main='Overall Survival')
hist(pvals[which(markers>thresh),2],breaks=30,xlab='P-Value',main='Time to Sig. Treat.')
hist(pvals[which(markers>thresh),3],breaks=30,xlab='P-Value',main='Progression Free Survival')
hist(pvals[which(markers>thresh),4],breaks=30,xlab='P-Value',main='Event Free Survival')
hist(pvals[which(markers>thresh),5],breaks=30,xlab='P-Value',main='Time to Progression')
hist(pvals[which(markers>thresh),6],breaks=30,xlab='P-Value',main='Overall Survival (after treatment)')
hist(response.pval[which(markers>thresh)],breaks=50,xlab='P-Value',main='Response to Treatment')
hist(response.oddsRatio[which(markers>thresh)],breaks=50,xlab='Coefficient',main='Response to Treatment')

#Multiple test correction:
cutoff.t2st <- cutoffSignificant(bt2t <- Bum(pvals[which(markers>thresh),2]), alpha=0.1, by="FDR")
cutoff.ttp <- cutoffSignificant(Bum(pvals[which(markers>thresh),5]), alpha=0.1, by="FDR")

#Looking at 'significant' segs:
segs.sig_t2st <- segInfo[[1]][which(pvals[,2]<cutoff.t2st & markers>thresh),1:3]
segs.sig_ttp <- segInfo[[1]][which(pvals[,5]<cutoff.ttp & markers>thresh),1:3]

#Figure 1: seg 29 (chr 2); seg 346 (chr 12); 476 (chr 13); seg 558 (chr 17); seg 617 (chr 19); seg 627 (chr 21). 
segs <- c(29,346,558,617,627)
ttp.elig <- TTP[ttp.eligible]
numericStatus <- NumericTTPStatus[ttp.eligible]

###Generating prodlim input:
f <- function(seg){
  afs <- unlist(sapply(1:length(segInfo),function(j){segInfo[[j]]$AF[seg]}))
  tcns <- unlist(sapply(1:length(segInfo),function(j){segInfo[[j]]$TCN[seg]}))
  afs[which(is.na(tcns))] <- NA
  af.status <- rep(NA, length(afs))
  notna <- which(!is.na(afs))
  if(mean(na.omit(afs))>0){
    af.status[notna] <- 0
    af.status[notna][which(afs[notna]>0)] <- 1
  }else{
    af.status[notna] <- 0
    af.status[notna][which(afs[notna]<0)] <- 1
  }
  af.status.var <- af.status[ttp.eligible]
  Altered <- af.status.var
  df <- data.frame('TTP'=ttp.elig,'numericStatus'=numericStatus,'Altered'=Altered)
  df
}

dfs <- lapply(c(346,558,29,617,627),f)
prod346 <- prodlim(Surv(TTP,numericStatus)~Altered,data=dfs[[1]])
prod558 <- prodlim(Surv(TTP,numericStatus)~Altered,data=dfs[[2]])
prod29 <- prodlim(Surv(TTP,numericStatus)~Altered,data=dfs[[3]])
prod617 <- prodlim(Surv(TTP,numericStatus)~Altered,data=dfs[[4]])
prod627 <- prodlim(Surv(TTP,numericStatus)~Altered,data=dfs[[5]])

#Figure 1: Known CNVs
par(mfrow=c(1,2))
plot(prod346,logrank=TRUE)
title('A) Trisomy 12',cex=1)
plot(prod558,logrank=TRUE)
title('B) 17p deletion',cex=1)

#Figure 2: Novel CNVs
par(mfrow=c(1,3))
plot(prod29,logrank=TRUE)
title('B) 2p gain',cex=1)
plot(prod617,logrank=TRUE)
title('B) 19p deletion',cex=1)
plot(prod627,logrank=TRUE)
title('B) 21p gain',cex=1)

#Figure 3: Alteration fractions (known and novel)
afmat <- t(sapply(1:length(segInfo),function(i){segInfo[[i]]$AF[segs]}))
par(mfrow=c(2,3))
plot(afmat[,1],ylim=c(0,1),xlab='Patient Index',ylab='Alteration Fraction',main='Alteration: 2P Gain')
plot(afmat[,2],ylim=c(0,1),xlab='Patient Index',ylab='Alteration Fraction',main='Alteration: Trisomy 12')
plot(afmat[,3],ylim=c(0,1),xlab='Patient Index',ylab='Alteration Fraction',main='Alteration: 17P Deletion')
plot(afmat[,4],ylim=c(0,1),xlab='Patient Index',ylab='Alteration Fraction',main='Alteration: 19P Deletion')
plot(afmat[,5],ylim=c(0,1),xlab='Patient Index',ylab='Alteration Fraction',main='Alteration: 21P Gain')

#Figure 4: Is clinical association due to presence vs. absence of CNV, or low vs. high alteration fraction?
#Not sure what to put here yet.

#Figure 5: Response to treatment.
#Not sure why, but R thinks some pvals given by 'fisher.test' that look like 1 are really >1.
response.pval[response.pval>1] <- 1
cutoff.rtt <- cutoffSignificant(Bum(response.pval[which(markers>thresh)]), alpha=0.1, by="FDR")
#Nothing passes?



###I know longer need this for figures 1 and 2, but here's prodlim being used within a function
#and the ensuing issues.
pfun <- function(seg){
  afs <- unlist(sapply(1:length(segInfo),function(j){segInfo[[j]]$AF[seg]}))
  tcns <- unlist(sapply(1:length(segInfo),function(j){segInfo[[j]]$TCN[seg]}))
  afs[which(is.na(tcns))] <- NA
  af.status <- rep(NA, length(afs))
  notna <- which(!is.na(afs))
  if(mean(na.omit(afs))>0){
    af.status[notna] <- 0
    af.status[notna][which(afs[notna]>0)] <- 1
  }else{
    af.status[notna] <- 0
    af.status[notna][which(afs[notna]<0)] <- 1
  }
  af.status.var <- af.status[ttp.eligible]
  Altered <- af.status.var
  df1 <- data.frame('TTP'=ttp.elig,'numericStatus'=numericStatus,'Altered'=Altered)
  prod <- prodlim(Surv(TTP,numericStatus)~Altered,data=df1)
  #plot(prod,logrank=TRUE,xlab='Time to Progression')
  #title(mainlab,cex=1)
  prod
}

par(mfrow=c(1,2))
pfun(346,mainlab='A) Trisomy 12')
pfun(558,mainlab='B) 17p deletion')

par(mfrow=c(1,3))
pfun(29,mainlab='A) 2p gain')
pfun(617,mainlab='B) 19p deletion')
pfun(627,mainlab='C) 21p gain')

temp29 <- pfun(29)
temp617 <- pfun(617)
temp627 <- pfun(627)

temp100 <- pfun(100)
plot(temp100,logrank=TRUE)

par(mfrow=c(1,3))
plot(temp617,logrank=TRUE)
plot(temp29,logrank=TRUE)
plot(temp627,logrank=TRUE)
