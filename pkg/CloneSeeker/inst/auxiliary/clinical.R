library(survival)
library(prodlim)
library(quantmod)
cll.path <- 'C:/Users/Mark/OneDrive - The Ohio State University/clonetools/res/res-CLL'
cll.files <- list.files(cll.path)
samples <- sapply(1:length(cll.files),function(j){strsplit(strsplit(cll.files[j],split='-')[[1]][4],split='.rda')[[1]][1]})
cf.results <- lapply(1:length(cll.files), function(j){get(load(paste(cll.path, '/', cll.files[j], sep='')))$res})
cll.ids <- sapply(1:length(cll.files),function(j){strsplit(strsplit(cll.files[j],split='-')[[1]][4],split='.rda')[[1]][1]})
pars <- cf.results[[1]]$pars
rewind <- get(load('C:/Users/Mark/OneDrive - The Ohio State University/clonetools/SNP_array_data/rewind_both.rda'))
total.markers <- sum(rewind[rewind$SamID==samples[1],]$num.mark)
clinical <- get(load('C:/Users/Mark/Documents/Lab/clinok.rda'))
#clinical <- get(load('E:/CLL/ClinicalOK.rda'))
eligible <- get(load('C:/Users/Mark/OneDrive - The Ohio State University/clinical/eligible.rda'))
clinical$ttp.eligible <- eligible
cids <- rownames(clinical)
common.cids <- intersect(c(cids),samples)
clinical <- clinical[sapply(1:length(common.cids),function(j){which(rownames(clinical)==common.cids[j])}),]
cf.results <- lapply(1:length(common.cids),function(j){cf.results[[which(samples==common.cids[j])]]})
ttp.eligible <- which(clinical$ttp.eligible)

#Table of characteristics: age, ethnicity, sex, Rai stage
race <- table(clinical$Race)
names(race) <- c('Asian','Black','Hispanic','White')
sex <- table(clinical$Sex)
names(sex) <- c('Female','Male')
mutstatus <- table(clinical$mutation.status)
names(mutstatus) <- c('Mutated','Unmutates')
rai <- table(clinical$Rai.Stage)
zap70 <- table(clinical$Zap70Protein)
median.age.dx <- median(clinical$AgeAtDx)
dohner <- table(clinical$Dohner)
timeSam2Treat <- median(clinical$TimeDiagnosis2SigTreat)
df <- data.frame('NumberofPatients'=c(race,sex,mutstatus,rai,zap70,median.age.dx,dohner,timeSam2Treat))
namesvec <- c(names(race), names(sex), names(mutstatus), names(rai), names(zap70), "MedianAgeatDiagnosis", names(dohner), "TimefromSampletoTreatment")
#rownames(df) <- namesvec
#save(df, file='C:/Users/Mark/OneDrive - The Ohio State University/clonetools/df.tsv')

###Common alteration CNV fractions:
#Seperate genome into common segments according to all segment breaks
sams <- unique(rewind$SamID)
chrs <- unique(rewind$chrom)

###
ks <- sapply(1:length(cf.results),function(j){length(which(cf.results[[j]]$psi>0))})
columns <- colnames(clinical)

alt.frac <- sapply(1:length(cf.results),function(j){
  temp <- rewind[rewind$SamID==common.cids[j],]
  notna <- which(!is.na(temp$loc.end) & !is.na(temp$loc.start))
  all.lens <- sum(as.numeric(temp$loc.end[notna]) - as.numeric(temp$loc.start[notna]))
  alt <- which(sapply(1:nrow(cf.results[[j]]$A),function(k){
    unique(all.equal(unname(cf.results[[j]]$A[k,1:length(which(cf.results[[j]]$psi>0))]), rep(1, length(which(cf.results[[j]]$psi>0))))!=TRUE) & 
      unique(all.equal(unname(cf.results[[j]]$B[k,1:length(which(cf.results[[j]]$psi>0))]), rep(1, length(which(cf.results[[j]]$psi>0))))!=TRUE)
  }))
  lens <- unlist(sapply(1:length(alt), function(k){
    index <- which(rewind$seg.median==cf.results[[j]]$filtered.data$cndata.filt$LRR[alt[k]] &
                   rewind$AvgBAF==cf.results[[j]]$filtered.data$cndata.filt$BAF[alt[k]])
    rewind$loc.end[index] - rewind$loc.start[index]
  }))
  sum(as.numeric(lens))/all.lens
})

k.status <- af.status <- rep(0,length(ks))
names(ks) <- names(alt.frac) <- common.cids
k.status[ks>1] <- 1
af.status[alt.frac>median(alt.frac)] <- 1

merge <- function(vector){
  if(length(vector)>1){
    starts <- rep(NA, length(vector))
    ends <- starts
    diffs <- diff(vector)
    stop <- FALSE
    curr <- 1
    index <- 1
    start <- vector[1]
    end <- vector[1]
    while(stop==FALSE) {
      if(diffs[index]==1){
        end <- vector[index+1]
      }else{
        starts[curr] <- start
        ends[curr] <- end
        start <- vector[index+1]
        end <- start
        curr <- curr + 1
      }
      if(index+1>=length(vector)){
        stop <- TRUE
        if(diffs[index]==1){
          starts[curr] <- start
          ends[curr] <- end
        }
      }
      index <- index+1
    }
    if(vector[length(vector)] - vector[length(vector)-1] >1){
      starts[length(na.omit(starts))+1] <- vector[length(vector)]
      ends[length(na.omit(ends))+1] <- vector[length(vector)]
    }else{
      ends[length(na.omit(ends))] <- vector[length(vector)]
    }
    df <- data.frame('start'=na.omit(starts), 'end'=na.omit(ends))
  }else{
    df <- data.frame('start'=vector[1], 'end'=vector[1])
  }
  df
}

marker.thresh <- 3000
frac.thresh <- .2
complex <- sapply(1:length(cf.results),function(j){
  segs <- nrow(cf.results[[j]]$A)
  psi <- cf.results[[j]]$psi
  markers <- cf.results[[j]]$filtered.data$cndata.filt$markers
  cns <- cf.results[[j]]$A + cf.results[[j]]$B
  gainfrac <- sapply(1:segs,function(k){sum(psi[which(cns[k,]>2)])})
  lossfrac <- sapply(1:segs,function(k){sum(psi[which(cns[k,]<2)])})
  merged.gain <- na.omit(merge(which(gainfrac>frac.thresh)))
  merged.loss <- na.omit(merge(which(lossfrac>frac.thresh)))
  gain.sizes <- numeric(0)
  loss.sizes <- numeric(0)
  if(nrow(merged.gain)>0){
    gain.sizes <- sapply(1:nrow(merged.gain),function(k){sum(markers[merged.gain$start[k]:merged.gain$end[k]])})
  }
  if(nrow(merged.loss)>0){
    loss.sizes <- sapply(1:nrow(merged.loss),function(k){sum(markers[merged.loss$start[k]:merged.loss$end[k]])})
  }
  length(which(gain.sizes>marker.thresh)) + length(which(loss.sizes>marker.thresh)) >= 3
})

###Treatment (significant):
#vars; note: time variables are in months (365.25/12).
White.blood.count <- clinical[,which(columns=='White.blood.count')]
Rai.Stage <- clinical[,which(columns=='Rai.Stage')]
TimeSample2SigTreat <- clinical[,which(columns=='TimeSample2SigTreat')]
TimeDiag2SigTreat <- clinical$TimeDiagnosis2SigTreat
PFS.time <- clinical[,which(columns=='PFS.time')]
EFS.time <- clinical[,which(columns=='EFS.time')]
mutation.status <- clinical[,which(columns=='mutation.status')]
OSAD <- clinical[,which(columns=='OSAfterDiagnosis')]
OSAS <- clinical[,which(columns=='OSAfterSample')]
TTP <- clinical[,which(columns=='TTP.time')]
TxDSS <- clinical[,which(columns=='TxDSS.time')]

#status (numeric)
NumericSigTreatment <- clinical[,which(columns=='NumericSigTreatment')]
PFS.numericStatus <- clinical[,which(columns=='PFS.numericStatus')]
EFS.numericStatus <- clinical[,which(columns=='EFS.numericStatus')]
NumericVitalStatus <- clinical[,which(columns=='NumericVitalStatus')]
NumericTTPStatus <- clinical[,which(columns=='TTP.numericStatus')]
NumericTxDSSStatus <- clinical[,which(columns=='TxDSS.numericStatus')]

########Cox Proportional Hazards:
colors <- c('red', 'blue', 'purple', 'darkorange')
Multiclonal <- k.status
Complex <- rep(NA, nrow(clinical))
Complex[which(clinical$CatCyto=='Simple')] <- 0
Complex[which(clinical$CatCyto=='Complex')] <- 1
High.AF <- af.status
Clones <- ks
Clones <- as.character(ks)
Clones[Clones=='3' | Clones=='4'] <- '3+'
Clones <- as.factor(Clones)
K <- Clones

###Time to significant treatment
t2st.karyo <- prodlim(Surv(TimeDiag2SigTreat, NumericSigTreatment) ~ Multiclonal + Complex)
diff <- survdiff(Surv(TimeDiag2SigTreat, NumericSigTreatment) ~ Multiclonal + Complex)
p_t2st.karyo <- round(pchisq(diff$chisq, length(diff$n)-1, lower.tail = FALSE), digits=4)
plot(t2st.karyo, logrank=TRUE)

t2st.AF <- prodlim(Surv(TimeDiag2SigTreat, NumericSigTreatment) ~ Multiclonal + High.AF)
diff <- survdiff(Surv(TimeDiag2SigTreat, NumericSigTreatment) ~ Multiclonal + High.AF)
p_t2st.AF <- round(pchisq(diff$chisq, length(diff$n)-1, lower.tail = FALSE), digits=4)
plot(t2st.AF, logrank=TRUE)

t2st <- prodlim(Surv(TimeDiag2SigTreat, NumericSigTreatment) ~ Multiclonal)
diff <- survdiff(Surv(TimeDiag2SigTreat, NumericSigTreatment) ~ Multiclonal)
p_t2st <- round(pchisq(diff$chisq, length(diff$n)-1, lower.tail = FALSE), digits=4)
plot(t2st, logrank=TRUE)

t2st.int <- prodlim(Surv(TimeDiag2SigTreat, NumericSigTreatment) ~ K)
diff <- survdiff(Surv(TimeDiag2SigTreat, NumericSigTreatment) ~ K)
p_t2st.int <- round(pchisq(diff$chisq, length(diff$n)-1, lower.tail = FALSE), digits=4)

###OSAS:
osas.karyo <- prodlim(Surv(OSAS, NumericVitalStatus) ~ Multiclonal + Complex)
diff <- survdiff(Surv(OSAS, NumericVitalStatus) ~ Multiclonal + Complex)
p_osas.karyo <- round(pchisq(diff$chisq, length(diff$n)-1, lower.tail = FALSE), digits=4)
plot(osas.karyo, logrank=TRUE)

osas.AF <- prodlim(Surv(OSAS, NumericVitalStatus) ~ Multiclonal + High.AF)
diff <- survdiff(Surv(OSAS, NumericVitalStatus) ~ Multiclonal + High.AF)
p_osas.AF <- round(pchisq(diff$chisq, length(diff$n)-1, lower.tail = FALSE), digits=4)
plot(osas.AF, logrank=TRUE)

osas <- prodlim(Surv(OSAS, NumericVitalStatus) ~ Multiclonal)
diff <- survdiff(Surv(OSAS, NumericVitalStatus) ~ Multiclonal)
p_osas <- round(pchisq(diff$chisq, length(diff$n)-1, lower.tail = FALSE), digits=4)
plot(osas, logrank=TRUE)

osas.int <- prodlim(Surv(OSAS, NumericVitalStatus) ~ K)
diff <- survdiff(Surv(OSAS, NumericVitalStatus) ~ K)
p_osas.int <- round(pchisq(diff$chisq, length(diff$n)-1, lower.tail = FALSE), digits=4)

###Progression free survival
ttp.time <- TTP[ttp.eligible]
TTP.status <- NumericTTPStatus[ttp.eligible]
Multiclonal <- Multiclonal[ttp.eligible]
Complex <- Complex[ttp.eligible]
High.AF <- High.AF[ttp.eligible]
K <- K[ttp.eligible]

ttp.karyo <- prodlim(Surv(ttp.time, TTP.status) ~ Multiclonal + Complex)
diff <- survdiff(Surv(ttp.time, TTP.status) ~ Multiclonal + Complex)
p_ttp.karyo <- round(pchisq(diff$chisq, length(diff$n)-1, lower.tail = FALSE), digits=4)
plot(ttp.karyo, logrank=TRUE)

ttp <- prodlim(Surv(ttp.time, TTP.status) ~ Multiclonal)
diff <- survdiff(Surv(pfs.time, TTP.status) ~ Multiclonal)
p_ttp <- round(pchisq(diff$chisq, length(diff$n)-1, lower.tail = FALSE), digits=4)
plot(ttp, logrank=TRUE)

ttp.int <- prodlim(Surv(ttp.time, TTP.status) ~ K)
diff <- survdiff(Surv(ttp.time, TTP.status) ~ K)
p_ttp.int <- round(pchisq(diff$chisq, length(diff$n)-1, lower.tail = FALSE), digits=4)

###OSAT:
osat.time <- TxDSS[ttp.eligible]
vitalStatus <- NumericVitalStatus[ttp.eligible]
osat.karyo <- prodlim(Surv(osat.time, vitalStatus) ~ Multiclonal + Complex)
diff <- survdiff(Surv(osat.time, vitalStatus) ~ Multiclonal + Complex)
p_osat.karyo <- round(pchisq(diff$chisq, length(diff$n)-1, lower.tail = FALSE), digits=4)
plot(osat.karyo, logrank=TRUE)

osat.AF <- prodlim(Surv(osat.time, vitalStatus) ~ Multiclonal + High.AF)
diff <- survdiff(Surv(osat.time, vitalStatus) ~ Multiclonal + High.AF)
p_osat.AF <- round(pchisq(diff$chisq, length(diff$n)-1, lower.tail = FALSE), digits=4)
plot(osat.AF, logrank=TRUE)

osat <- prodlim(Surv(osat.time, vitalStatus) ~ Multiclonal)
diff <- survdiff(Surv(osat.time, vitalStatus) ~ Multiclonal)
p_osat <- round(pchisq(diff$chisq, length(diff$n)-1, lower.tail = FALSE), digits=4)
plot(osat, logrank=TRUE)

osat.int <- prodlim(Surv(osat.time, vitalStatus) ~ K)
diff <- survdiff(Surv(osat.time, vitalStatus) ~ K)
p_osat.int <- round(pchisq(diff$chisq, length(diff$n)-1, lower.tail = FALSE), digits=4)


###Combined plots
K <- Clones
png('C:/Users/Mark/OneDrive - The Ohio State University/clonetools/cf_images/surv_int.png',height=700,width=2100,res=140)
par(mfrow=c(1,3))
plot(t2st.int, logrank=TRUE)
title("Time to Treatment",cex=1.5)
K <- K[ttp.eligible]
plot(pfs.int, logrank=TRUE)
title("Progression-Free Survival",cex=1.5)
plot(osat.int, logrank=TRUE)
title("Overall Survival after Treatment",cex=1.5)
dev.off()

#Among 2-cloners, does the size of the smaller subclone relate to clinical outcome?
indices <- which(ks==2)
psis <- t(sapply(1:length(cf.results),function(j){cf.results[[j]]$psi}))
minor.psis <- psis[,2]
minor.psi.status <- rep(0,length(minor.psis))
cutoff <- .2
minor.psi.status[which(minor.psis>cutoff)] <- 1
minor.psi.status2 <- minor.psi.status[indices]
survosas2 <- OSAS[indices]
survt2st2 <- TimeDiag2SigTreat[indices]
survttp2 <- TTP[intersect(indices,ttp.eligible)]
nvs <- NumericVitalStatus[indices]
nst <- NumericSigTreatment[indices]
nttps <- NumericTTPStatus[intersect(indices,ttp.eligible)]
minor.eligible <- minor.psis[intersect(indices,ttp.eligible)]
minor.eligible.status2 <- minor.psi.status[intersect(indices,ttp.eligible)]

osas.2 <- prodlim(Surv(survosas2, nvs)~minor.psi.status2)
t2st.2 <- prodlim(Surv(survt2st2,nst)~minor.psi.status2)
ttp.2 <- prodlim(Surv(survttp2,nttps)~minor.eligible.status2)
plot(osas.2,logrank =TRUE)
plot(t2st.2,logrank =TRUE)
plot(ttp.2,logrank =TRUE)

