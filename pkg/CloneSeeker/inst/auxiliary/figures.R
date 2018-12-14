library(prodlim)
library(survival)
###Figure 1:
plotSetTraits <- function(setDir, seq=TRUE, cn=TRUE, merge=FALSE, kmax=5){
  simfiles <- list.files(setDir)
  metadata <- simfiles[1]
  simfiles <- simfiles[2:length(simfiles)]
  sims <- lapply(1:length(simfiles),function(j){get(load(paste(setDir, '\\', simfiles[j], sep='')))})
  if(length(sims[[1]])==1){
    sims <- lapply(1:length(sims),function(x){sims[[x]]$tumor})
  }
  ks <- sapply(1:length(sims),function(j){length(which(sims[[j]]$psi>0))})
  psis <- t(sapply(1:length(sims),function(j){
    psi <- sims[[j]]$psi
    if(length(psi)<kmax){
      psi <- c(psi, rep(0, kmax-length(psi)))
    }
    psi
  }))
  frac.altered <- sapply(1:length(sims),function(j){
    sim <- sims[[j]]
    lengths <- sim$clones[[1]]$cn$end + 1 - sim$clones[[1]]$cn$start
    sum(lengths[unique(unlist(lapply(1:length(sim$clones),function(x){
      which(sim$clones[[x]]$cn$A!=1 | sim$clones[[x]]$cn$B!=1)})))])/sum(lengths)
  })
  if(merge){
    sim <- sims[[j]]
    lengths <- sim$clones[[1]]$cn$end + 1 - sim$clones[[1]]$cn$start
    altered <- sapply(1:length(sim$clones),function(x){sim$clones[[x]]$cn$A!=1 | sim$clones[[x]]$cn$B!=1})
    altered.indices <- which(altered)
    diffs <- diff(altered.indices)
    #To be continued?
  }else{
    length.altered <- lapply(1:length(sims),function(j){
      sim <- sims[[j]]
      lengths <- sim$clones[[1]]$cn$end + 1 - sim$clones[[1]]$cn$start
      lengths[unique(unlist(lapply(1:length(sim$clones),function(x){
        which(sim$clones[[x]]$cn$A!=1 | sim$clones[[x]]$cn$B!=1)})))]
    })
    nCNVs <- sapply(1:length(length.altered),function(z){length(length.altered[[z]])}) 
  }
  length.altered <- unlist(length.altered)
  mutations <- sapply(1:length(sims),function(j){
    sim <- sims[[j]]
    if(length(sim)==1){
      sim <- sim$tumor
    }
    length(unique(unlist(lapply(1:length(sim$clones), function(x){sim$clones[[x]]$seq$mut.id}))))
  })
  if(cn==TRUE & seq==FALSE){
    plot(density(unlist(psis),bw=.05), xlab='Psi', xlim=c(0,1), main='F) Density of Psi Values across Simulations')
    plot(density(frac.altered,bw=.01), xlab='Genome Altered Fraction', xlim=c(0,1), main='G) Genome AF Density among Simulations')
  }else if(cn==FALSE & seq==TRUE){
    plot(density(unlist(psis),bw=.05), xlab='Psi', xlim=c(0,1), main='A) Density of Psi Values across Simulations')
    plot(density(mutations,bw=5), xlab='Mutations per Sample',xlim=c(0,max(mutations)), main='B) Density of Mutations per Simulation')
  }else if(cn & seq){
    plot(density(unlist(psis),bw=.05), xlab='Psi', xlim=c(0,1), main='C) Density of Psi Values across Simulations')
    plot(density(frac.altered,bw=.01), xlab='Genome Altered Fraction', xlim=c(0,1), main='D) Genome AF Density among Simulations')
    plot(density(mutations,bw=5), xlab='Mutations per Sample',xlim=c(0,max(mutations)), main='E) Density of Mutations per Simulation')
  }
}

png('C:/Users/Mark/OneDrive - The Ohio State University/clonetools/cf_images/figure1.png',height=500,width=900,res=100)
layout(matrix(c(1,1,1,2,2,2,3,3,4,4,5,5,6,6,6,7,7,7), 3, 6, byrow = TRUE))
plotSetTraits('C:\\Users\\Mark\\OneDrive - The Ohio State University\\clonetools\\sim\\sims-1', cn=FALSE)
plotSetTraits('C:\\Users\\Mark\\OneDrive - The Ohio State University\\clonetools\\sim\\sims-2')
plotSetTraits('C:\\Users\\Mark\\OneDrive - The Ohio State University\\clonetools\\sim\\sims-3', seq=FALSE)
dev.off()


#Figure S1 and S2:
metadata <- get(load('C:\\Users\\Mark\\OneDrive - The Ohio State University\\clonetools\\sim\\sims-3\\metadata-3.rda'))
plotSNP <- function(chr.sim, chr.dat, j, sam, homFrac=.68, mu.0=0, sigma.0=.03){
  dat <- get(load(paste('C:\\Users\\Mark\\OneDrive - The Ohio State University\\clonetools\\dat\\data-3', '\\dat-3-', j,'.rda', sep='')))
  sim <- get(load(paste('C:\\Users\\Mark\\OneDrive - The Ohio State University\\clonetools\\sim\\sims-3', '\\sim-3-', j,'.rda', sep='')))
  segdata <- dat$dat$cn.data
  #sigma0.lrr <- metadata$data.params$sigma0.lrr
  sigma0.lrr <- .12
  sigma0.baf <- metadata$data.params$sigma0.baf
  segs <- segdata[segdata$chr==chr.sim,]
  x.mu <- unlist(sapply(1:nrow(segs),function(k){rep(segs$X[k], segs$markers[k])}))
  y.mu <- unlist(sapply(1:nrow(segs),function(k){rep(segs$Y[k], segs$markers[k])}))
  index <- 1:length(x.mu)
  hom <- sample(1:length(index), round(length(index)*homFrac), replace=FALSE)
  inverted <- sample(1:length(index), round(length(index)*.5), replace=FALSE)
  x.mu[hom] <- x.mu[hom] + y.mu[hom]
  y.mu[hom] <- 0
  x.inverted <- x.mu[inverted]
  x.mu[inverted] <- y.mu[inverted]
  y.mu[inverted] <- x.inverted
  lrr.mu <- log10((x.mu + y.mu)/2)
  baf.mu <- x.mu/(x.mu + y.mu)
  lrr <- rnorm(length(lrr.mu), lrr.mu, sigma0.lrr)
  baf <- rnorm(length(baf.mu), baf.mu, sigma0.baf)
  baf[which(baf<0)] <- -baf[which(baf<0)]
  baf[which(baf>1)] <- 1/baf[which(baf>1)]
  path <- paste('E:\\SNP-HapMap\\MergeByChrom\\Full Data Table_Chr_', chr.dat, '.txt', sep='')
  cndata <- read.table(path, header=T)
  columns <- unique(sapply(4:ncol(cndata),function(x){strsplit(colnames(cndata)[x],split='[.]')[[1]][1]}))
  samn <- which(columns==sam)
  lrr.dat <- cndata[,samn*3+2]
  baf.dat <- cndata[,samn*3+3]
  par(mfrow=c(2,2))
  #maintab <- paste('Chr ', chr, sep='')
  plot(1:length(lrr), lrr, xlab='SNP Index', ylab='LRR', ylim=c(-1.5,1), pch='.', col='gray',
       main='A) SNP Array Data (simulated data)')
  plot(1:length(lrr.dat), lrr.dat, xlab='SNP Index', ylab='LRR', ylim=c(-1.5,1), pch='.', col='gray',
       main='B) SNP Array Data (real data)')
  abline(h=c(log10((1:4)/2)), col=c('gray','black','gray','gray'))
  abline(h=c(log10((1:4)/2)), col=c('gray','black','gray','gray'))
  plot(1:length(baf), baf, main='', xlab='SNP Index', ylab='BAF', ylim=c(0,1), pch='.', col='gray')
  abline(h=.5, col='black')
  plot(1:length(baf.dat), baf.dat, main='', xlab='SNP Index', ylab='BAF', ylim=c(0,1), pch='.', col='gray')
  abline(h=.5, col='black')
  #plot(density(baf), xlim=c(-1.5,1.5), ylim=c(0,5), main=paste('SD=',round(sd(baf),digits=3), sep=''))
}

#Figure S2
#Note: I only plot the first 225 SNVs in the function below because (unlike with the SNP array), the number of observations
#varies between samples, so I picked a number such that whatever samples or simulations I plot, they will both have at least
#This many SNVs. 
plotSEQ <- function(seqdata, seqsim){
  som.dat <- seqdata[seqdata$somatic_status=='SOMATIC',]
  vaf.dat <- na.omit(som.dat$alt_read_count/(som.dat$alt_read_count + som.dat$ref_read_count))
  tc.dat <- som.dat$alt_read_count + som.dat$ref_read_count
  seq.sim <- seqsim$dat$seq.data
  som.sim <- seq.sim[seq.sim$status=='somatic',]
  tc.sim <- som.sim$totalCounts
  vaf.sim <- som.sim$VAF
  par(mfrow=c(2,2))
  plot(tc.sim[1:225], ylim=c(0,250), xlab='Index', ylab='Total Read Counts', main='A) Somatic Variants (Simulation)', pch=16)
  plot(tc.dat[1:225], ylim=c(0, 250), xlab='Index', ylab='Total Read Counts', main='B) Somatic Variants (Data)', pch=16)
  plot(vaf.sim[1:225], ylim=c(0, 1), xlab='Index', ylab='VAF', main='', pch=16)
  plot(vaf.dat[1:225], ylim=c(0, 1), xlab='Index', ylab='VAF', main='', pch=16)
}

png('C:/Users/Mark/OneDrive - The Ohio State University/clonetools/cf_images/figureS1.png',width=1000,height=500)
plotSNP(chr.sim=2, j=68, homFrac=.68, mu.0=0, sigma.0=.03, chr.dat=2, sam='NA18855')
dev.off()
png('C:/Users/Mark/OneDrive - The Ohio State University/clonetools/cf_images/figureS2.png',width=1000,height=500)
plotSEQ(seqdata=get(load('C:\\Users\\Mark\\OneDrive - The Ohio State University\\deepcnv\\GS1115.rda')),
        seqsim=get(load('C:\\Users\\Mark\\OneDrive - The Ohio State University\\clonetools\\dat\\data-1\\dat-1-10.rda')))
dev.off()
#Figure S3:Karyotype plot:
###Not sure whether to include###

#Figure 2: Just mutations
assessments1 <- get(load('C:/Users/Mark/OneDrive - The Ohio State University/clonetools/assess/assess-1/assessments-1.rda'))
kdiff.cf1 <- na.omit(sapply(1:length(assessments1),function(j){assessments1[[j]][[1]]$kdiff}))
psi.dist.cf1 <- na.omit(sapply(1:length(assessments1),function(j){assessments1[[j]][[1]]$psi.dist}))
cn.dist.cf1 <- na.omit(sapply(1:length(assessments1),function(j){assessments1[[j]][[1]]$cn.dist}))
mut.dist.cf1 <- na.omit(sapply(1:length(assessments1),function(j){assessments1[[j]][[1]]$mut.dist}))
kdiff.sc1 <- na.omit(sapply(1:length(assessments1),function(j){assessments1[[j]][[2]]$kdiff}))
psi.dist.sc1 <- na.omit(sapply(1:length(assessments1),function(j){assessments1[[j]][[2]]$psi.dist}))
cn.dist.sc1 <- na.omit(sapply(1:length(assessments1),function(j){assessments1[[j]][[2]]$cn.dist}))
mut.dist.sc1 <- na.omit(sapply(1:length(assessments1),function(j){assessments1[[j]][[2]]$mut.dist}))
kdiff.exp1 <- na.omit(sapply(1:length(assessments1),function(j){assessments1[[j]][[3]]$kdiff}))
psi.dist.exp1 <- na.omit(sapply(1:length(assessments1),function(j){assessments1[[j]][[3]]$psi.dist}))
cn.dist.exp1 <- na.omit(sapply(1:length(assessments1),function(j){assessments1[[j]][[3]]$cn.dist}))
mut.dist.exp1 <- na.omit(sapply(1:length(assessments1),function(j){assessments1[[j]][[3]]$mut.dist}))
kdiff.exp1 <- na.omit(sapply(1:length(assessments1),function(j){assessments1[[j]][[3]]$kdiff}))
psi.cdf.cf1 <- list('x'=seq(from=min(psi.dist.cf1),to=max(psi.dist.cf1),length=150),'y'=sapply(seq(from=min(psi.dist.cf1),
   to=max(psi.dist.cf1),length=150),function(i){length(which(psi.dist.cf1<=i))/length(psi.dist.cf1)}))
cn.cdf.cf1 <- list('x'=seq(from=min(cn.dist.cf1),to=max(cn.dist.cf1),length=150),'y'=sapply(seq(from=min(cn.dist.cf1),
   to=max(cn.dist.cf1),length=150),function(i){length(which(cn.dist.cf1<=i))/length(cn.dist.cf1)}))
mut.cdf.cf1 <- list('x'=seq(from=min(mut.dist.cf1),to=max(mut.dist.cf1),length=150),'y'=sapply(seq(from=min(mut.dist.cf1),
   to=max(mut.dist.cf1),length=150),function(i){length(which(mut.dist.cf1<=i))/length(mut.dist.cf1)}))
psi.cdf.sc1 <- list('x'=seq(from=min(psi.dist.sc1),to=max(psi.dist.sc1),length=150),'y'=sapply(seq(from=min(psi.dist.sc1),
   to=max(psi.dist.sc1),length=150),function(i){length(which(psi.dist.sc1<=i))/length(psi.dist.sc1)}))
cn.cdf.sc1 <- list('x'=seq(from=min(cn.dist.sc1),to=max(cn.dist.sc1),length=150),'y'=sapply(seq(from=min(cn.dist.sc1),
   to=max(cn.dist.sc1),length=150),function(i){length(which(cn.dist.sc1<=i))/length(cn.dist.sc1)}))
mut.cdf.sc1 <- list('x'=seq(from=min(mut.dist.sc1),to=max(mut.dist.sc1),length=150),'y'=sapply(seq(from=min(mut.dist.sc1),
   to=max(mut.dist.sc1),length=150),function(i){length(which(mut.dist.sc1<=i))/length(mut.dist.sc1)}))
psi.cdf.exp1 <- list('x'=seq(from=min(psi.dist.exp1),to=max(psi.dist.exp1),length=150),'y'=sapply(seq(from=min(psi.dist.exp1),
   to=max(psi.dist.exp1),length=150),function(i){length(which(psi.dist.exp1<=i))/length(psi.dist.exp1)}))
cn.cdf.exp1 <- list('x'=seq(from=min(cn.dist.exp1),to=max(cn.dist.exp1),length=150),'y'=sapply(seq(from=min(cn.dist.exp1),
   to=max(cn.dist.exp1),length=150),function(i){length(which(cn.dist.exp1<=i))/length(cn.dist.exp1)}))
mut.cdf.exp1 <- list('x'=seq(from=min(mut.dist.exp1),to=max(mut.dist.exp1),length=150),'y'=sapply(seq(from=min(mut.dist.exp1),
   to=max(mut.dist.exp1),length=150),function(i){length(which(mut.dist.exp1<=i))/length(mut.dist.exp1)}))
xrange.psi1=c(min(psi.cdf.cf1$x,psi.cdf.sc1$x,psi.cdf.exp1$x),max(psi.cdf.cf1$x,psi.cdf.sc1$x,psi.cdf.exp1$x))
interval.psi1 <- max(xrange.psi1)-min(xrange.psi1)
xrange.cn1=c(min(cn.cdf.cf1$x,cn.cdf.sc1$x,cn.cdf.exp1$x),max(cn.cdf.cf1$x,cn.cdf.sc1$x,cn.cdf.exp1$x))
interval.cn1 <- max(xrange.cn1)-min(xrange.cn1)
xrange.mut1=c(min(mut.cdf.cf1$x,mut.cdf.sc1$x,mut.cdf.exp1$x),max(mut.cdf.cf1$x,mut.cdf.sc1$x,mut.cdf.exp1$x))
interval.mut1 <- max(xrange.mut1)-min(xrange.mut1)

png('C:/Users/Mark/OneDrive - The Ohio State University/clonetools/cf_images/figure2_new.png',width=1400,height=2000,res=200)
layout(matrix(c(1,2,3,4,4,4,5,5,5,6,6,6),nrow=6,ncol=2))
hist(kdiff.cf1, breaks=seq(from=-5.5,to=5.5, by=1),xlab='Difference', main='A) Histogram of K difference (Clonefinder)')
hist(kdiff.sc1, breaks=seq(from=-5.5,to=5.5, by=1),xlab='Difference', main='Histogram of K difference (SciClone)')
hist(kdiff.exp1, breaks=seq(from=-5.5,to=5.5, by=1),xlab='Difference', main='Histogram of K difference (Expands)')
plot(psi.cdf.cf1$x,psi.cdf.cf1$y, ylim=c(0,1), main='B)  Psi Distance Cumulative Distribution',type='l',
     xlim=c(min(psi.cdf.cf1$x,psi.cdf.sc1$x,psi.cdf.exp1$x),max(psi.cdf.cf1$x,psi.cdf.sc1$x,psi.cdf.exp1$x)),
     xlab='Psi Distance',ylab='Cumulative Density')
lines(psi.cdf.sc1$x,psi.cdf.sc1$y,col='blue')
lines(psi.cdf.exp1$x,psi.cdf.exp1$y,col='red')
legend('bottomright', c('Clonefinder','Sciclone','Expands'),col=c('black','blue','red'),lty=c(1,1),cex=.7)
plot(mut.cdf.cf1$x,mut.cdf.cf1$y, ylim=c(0,1), main='C) Mutation Distance Cumulative Distribution',type='l',
     xlim=c(min(mut.cdf.cf1$x,mut.cdf.sc1$x,mut.cdf.exp1$x),max(mut.cdf.cf1$x,mut.cdf.sc1$x,mut.cdf.exp1$x)),
     xlab='Mutation Distance',ylab='Cumulative Density')
lines(mut.cdf.sc1$x,mut.cdf.sc1$y,col='blue')
lines(mut.cdf.exp1$x,mut.cdf.exp1$y,col='red')
legend('bottomright', c('Clonefinder','Sciclone','Expands'),col=c('black','blue','red'),lty=c(1,1),cex=1)
dev.off()


#Figure 3 - New (CNV-sim and mixtures)
assessments2 <- get(load('C:/Users/Mark/OneDrive - The Ohio State University/clonetools/assess/assess-2/assessments-2.rda'))
kdists.cf2 <- na.omit(sapply(1:length(assessments2),function(j){assessments2[[j]][[1]]$k.dist}))
psi.dist.cf2 <- na.omit(sapply(1:length(assessments2),function(j){assessments2[[j]][[1]]$psi.dist}))
cn.dist.cf2 <- na.omit(sapply(1:length(assessments2),function(j){assessments2[[j]][[1]]$cn.dist}))
kdiff2 <- na.omit(sapply(1:length(assessments2),function(j){assessments2[[j]][[1]]$kdiff}))
cndiff.avg2 <- na.omit(sapply(1:length(assessments2),function(j){assessments2[[j]][[1]]$avg.cndiff}))
psi.cdf.cf2 <- list('x'=seq(from=min(psi.dist.cf2),to=max(psi.dist.cf2),length=150),'y'=sapply(seq(from=min(psi.dist.cf2),
    to=max(psi.dist.cf2),length=150),function(i){length(which(psi.dist.cf2<=i))/length(psi.dist.cf2)}))
cn.cdf.cf2 <- list('x'=seq(from=min(cn.dist.cf2),to=max(cn.dist.cf2),length=150),'y'=sapply(seq(from=min(cn.dist.cf2),
    to=max(cn.dist.cf2),length=150),function(i){length(which(cn.dist.cf2<=i))/length(cn.dist.cf2)}))

assessmentsMix <- get(load('C:/Users/Mark/OneDrive - The Ohio State University/clonetools/assess/assess-mix/assessments-mix.rda'))
kdists.cfMix <- na.omit(sapply(1:length(assessmentsMix),function(j){assessmentsMix[[j]][[1]]$k.dist}))
psi.dist.cfMix <- na.omit(sapply(1:length(assessmentsMix),function(j){assessmentsMix[[j]][[1]]$psi.dist}))
cn.dist.cfMix <- na.omit(sapply(1:length(assessmentsMix),function(j){assessmentsMix[[j]][[1]]$cn.dist}))
kdiffMix <- na.omit(sapply(1:length(assessmentsMix),function(j){assessmentsMix[[j]][[1]]$kdiff}))
cndiff.avgMix <- na.omit(sapply(1:length(assessmentsMix),function(j){assessmentsMix[[j]][[1]]$avg.cndiff}))
psi.cdf.cfMix <- list('x'=seq(from=min(psi.dist.cfMix),to=max(psi.dist.cfMix),length=150),'y'=sapply(seq(from=min(psi.dist.cfMix),
       to=max(psi.dist.cfMix),length=150),function(i){length(which(psi.dist.cfMix<=i))/length(psi.dist.cfMix)}))
cn.cdf.cfMix <- list('x'=seq(from=min(cn.dist.cfMix),to=max(cn.dist.cfMix),length=150),'y'=sapply(seq(from=min(cn.dist.cfMix),
       to=max(cn.dist.cfMix),length=150),function(i){length(which(cn.dist.cfMix<=i))/length(cn.dist.cfMix)}))

png('C:/Users/Mark/OneDrive - The Ohio State University/clonetools/cf_images/figure3_new.png',width=1500,height=2000,res=150)
par(mfrow=c(3,2))
hist(kdiff2, breaks=seq(from=-5.5,to=5.5, by=1),xlab='Difference', main='A) Histogram of K difference')
hist(kdiffMix, breaks=seq(from=-5.5,to=5.5, by=1),xlab='Difference', main='D) Histogram of K difference')
plot(psi.cdf.cf2$x,psi.cdf.cf2$y, ylim=c(0,1), main='B) Psi Distance Cumulative Distribution',type='l',
     xlim=c(min(psi.cdf.cf2$x),max(psi.cdf.cf2$x)),
     xlab='Psi Distance',ylab='Cumulative Density')
plot(psi.cdf.cfMix$x,psi.cdf.cfMix$y, ylim=c(0,1), main='E) Psi Distance Cumulative Distribution',type='l',
     xlim=c(min(psi.cdf.cfMix$x),max(psi.cdf.cfMix$x)),
     xlab='Psi Distance',ylab='Cumulative Density')
plot(cn.cdf.cf2$x,cn.cdf.cf2$y, ylim=c(0,1), main='C) Copy Number Distance Cumulative Distribution',type='l',
     xlim=c(min(cn.cdf.cf2$x),max(cn.cdf.cf2$x)),
     xlab='Copy Number Distance',ylab='Cumulative Density')
plot(cn.cdf.cfMix$x,cn.cdf.cfMix$y, ylim=c(0,1), main='F) Copy Number Distance Cumulative Distribution',type='l',
     xlim=c(min(cn.cdf.cfMix$x),max(cn.cdf.cfMix$x)),
     xlab='Copy Number Distance',ylab='Cumulative Density')
dev.off()


#Figure 4 New (sim-CNV+Mut, 2X2)
assessments3 <- get(load('C:/Users/Mark/OneDrive - The Ohio State University/clonetools/assess/assess-3/assessments-3.rda'))
cn.concordance.cf <- na.omit(sapply(1:length(assessments3),function(j){assessments3[[j]][[1]]$concord.cn}))
mut.concordance.cf <- na.omit(sapply(1:length(assessments3),function(j){assessments3[[j]][[1]]$concord.mut}))
psis.concordance.cf <- na.omit(sapply(1:length(assessments3),function(j){assessments3[[j]][[1]]$concord.psi}))
kdists.cf <- na.omit(sapply(1:length(assessments3),function(j){assessments3[[j]][[1]]$k.dist}))
psi.dist.cf <- na.omit(sapply(1:length(assessments3),function(j){assessments3[[j]][[1]]$psi.dist}))
cn.dist.cf <- na.omit(sapply(1:length(assessments3),function(j){assessments3[[j]][[1]]$cn.dist}))
mut.dist.cf <- na.omit(sapply(1:length(assessments3),function(j){assessments3[[j]][[1]]$mut.dist}))
kdiff.cf <- na.omit(sapply(1:length(assessments3),function(j){assessments3[[j]][[1]]$kdiff}))
psi.cdf.cf <- list('x'=seq(from=min(psi.dist.cf),to=max(psi.dist.cf),length=150),'y'=sapply(seq(from=min(psi.dist.cf),
                                                                                                to=max(psi.dist.cf),length=150),function(i){length(which(psi.dist.cf<=i))/length(psi.dist.cf)}))
cn.cdf.cf <- list('x'=seq(from=min(cn.dist.cf),to=max(cn.dist.cf),length=150),'y'=sapply(seq(from=min(cn.dist.cf),
                                                                                             to=max(cn.dist.cf),length=150),function(i){length(which(cn.dist.cf<=i))/length(cn.dist.cf)}))
mut.cdf.cf <- list('x'=seq(from=min(mut.dist.cf),to=max(mut.dist.cf),length=150),'y'=sapply(seq(from=min(mut.dist.cf),
                                                                                                to=max(mut.dist.cf),length=150),function(i){length(which(mut.dist.cf<=i))/length(mut.dist.cf)}))

psi.dist.sc <- na.omit(sapply(1:length(assessments3),function(j){assessments3[[j]][[2]]$psi.dist}))
cn.dist.sc <- na.omit(sapply(1:length(assessments3),function(j){assessments3[[j]][[2]]$cn.dist}))
mut.dist.sc <- na.omit(sapply(1:length(assessments3),function(j){assessments3[[j]][[2]]$mut.dist}))
kdiff.sc <- na.omit(sapply(1:length(assessments3),function(j){assessments3[[j]][[2]]$kdiff}))
psi.cdf.sc <- list('x'=seq(from=min(psi.dist.sc),to=max(psi.dist.sc),length=150),'y'=sapply(seq(from=min(psi.dist.sc),
                                                                                                to=max(psi.dist.sc),length=150),function(i){length(which(psi.dist.sc<=i))/length(psi.dist.sc)}))
cn.cdf.sc <- list('x'=seq(from=min(cn.dist.sc),to=max(cn.dist.sc),length=150),'y'=sapply(seq(from=min(cn.dist.sc),
                                                                                             to=max(cn.dist.sc),length=150),function(i){length(which(cn.dist.sc<=i))/length(cn.dist.sc)}))
mut.cdf.sc <- list('x'=seq(from=min(mut.dist.sc),to=max(mut.dist.sc),length=150),'y'=sapply(seq(from=min(mut.dist.sc),
                                                                                                to=max(mut.dist.sc),length=150),function(i){length(which(mut.dist.sc<=i))/length(mut.dist.sc)}))

psi.dist.exp <- na.omit(sapply(1:length(assessments3),function(j){assessments3[[j]][[3]]$psi.dist}))
cn.dist.exp <- na.omit(sapply(1:length(assessments3),function(j){assessments3[[j]][[3]]$cn.dist}))
mut.dist.exp <- na.omit(sapply(1:length(assessments3),function(j){assessments3[[j]][[3]]$mut.dist}))
kdiff.exp <- na.omit(sapply(1:length(assessments3),function(j){assessments3[[j]][[3]]$kdiff}))
psi.cdf.exp <- list('x'=seq(from=min(psi.dist.exp),to=max(psi.dist.exp),length=150),'y'=sapply(seq(from=min(psi.dist.exp),
                                                                                                   to=max(psi.dist.exp),length=150),function(i){length(which(psi.dist.exp<=i))/length(psi.dist.exp)}))
cn.cdf.exp <- list('x'=seq(from=min(cn.dist.exp),to=max(cn.dist.exp),length=150),'y'=sapply(seq(from=min(cn.dist.exp),
                                                                                                to=max(cn.dist.exp),length=150),function(i){length(which(cn.dist.exp<=i))/length(cn.dist.exp)}))
mut.cdf.exp <- list('x'=seq(from=min(mut.dist.exp),to=max(mut.dist.exp),length=150),'y'=sapply(seq(from=min(mut.dist.exp),
                                                                                                   to=max(mut.dist.exp),length=150),function(i){length(which(mut.dist.exp<=i))/length(mut.dist.exp)}))
xrange.psi=c(min(psi.cdf.cf1$x,psi.cdf.sc1$x,psi.cdf.exp1$x),max(psi.cdf.cf1$x,psi.cdf.sc1$x,psi.cdf.exp1$x))
interval.psi <- max(xrange.psi)-min(xrange.psi)
xrange.cn=c(min(cn.cdf.cf$x,cn.cdf.sc$x,cn.cdf.exp$x),max(cn.cdf.cf$x,cn.cdf.sc$x,cn.cdf.exp$x))
interval.cn <- max(xrange.cn)-min(xrange.cn)
xrange.mut=c(min(mut.cdf.cf$x,mut.cdf.sc$x,mut.cdf.exp$x),max(mut.cdf.cf$x,mut.cdf.sc$x,mut.cdf.exp$x))
interval.mut <- max(xrange.mut)-min(xrange.mut)

png('C:/Users/Mark/OneDrive - The Ohio State University/clonetools/cf_images/figure4_new.png',height=1900,width=2000,res=220)
layout(matrix(c(1,2,3,4,4,4,5,5,5,6,6,6),nrow=6,ncol=2))
hist(kdiff.cf, breaks=seq(from=-5.5,to=5.5, by=1),xlab='Difference', main='A) Histogram of K difference')
hist(kdiff.sc, breaks=seq(from=-5.5,to=5.5, by=1),xlab='Difference', main='Histogram of K difference')
hist(kdiff.exp, breaks=seq(from=-5.5,to=5.5, by=1),xlab='Difference', main='Histogram of K difference')
plot(psi.cdf.cf$x,psi.cdf.cf$y, ylim=c(0,1), main='B) Psi Distance Cumulative Distribution',type='l',
     xlim=c(min(psi.cdf.cf$x,psi.cdf.sc$x,psi.cdf.exp$x),max(psi.cdf.cf$x,psi.cdf.sc$x,psi.cdf.exp$x)),
     xlab='Psi Distance',ylab='Cumulative Density')
lines(psi.cdf.sc$x,psi.cdf.sc$y,col='blue')
lines(psi.cdf.exp$x,psi.cdf.exp$y,col='red')
legend('bottomright', c('Clonefinder','Sciclone','Expands'),col=c('black','blue','red'),lty=c(1,1),cex=1)
plot(cn.cdf.cf$x,cn.cdf.cf$y, ylim=c(0,1), main='C) Copy Number Distance Cumulative Distribution',type='l',
     xlim=c(min(cn.cdf.cf$x,cn.cdf.sc$x,cn.cdf.exp$x),max(cn.cdf.cf$x,cn.cdf.sc$x,cn.cdf.exp$x)),
     xlab='Copy Number Distance',ylab='Cumulative Density')
lines(cn.cdf.sc$x,cn.cdf.sc$y,col='blue')
lines(cn.cdf.exp$x,cn.cdf.exp$y,col='red')
legend('bottomright', c('Clonefinder','Sciclone','Expands'),col=c('black','blue','red'),lty=c(1,1),cex=1)
plot(mut.cdf.cf$x,mut.cdf.cf$y, ylim=c(0,1), main='D) Mutation Distance Cumulative Distribution',type='l',
     xlim=c(min(mut.cdf.cf$x,mut.cdf.sc$x,mut.cdf.exp$x),max(mut.cdf.cf$x,mut.cdf.sc$x,mut.cdf.exp$x)),
     xlab='Mutation Distance',ylab='Cumulative Density')
lines(mut.cdf.sc$x,mut.cdf.sc$y,col='blue')
lines(mut.cdf.exp$x,mut.cdf.exp$y,col='red')
legend('bottomright', c('Clonefinder','Sciclone','Expands'),col=c('black','blue','red'),lty=c(1,1),cex=1)
dev.off()


###Figure S4 (?):
res.path <- 'C:/Users/Mark/OneDrive - The Ohio State University/clonetools/res/res-CLL'
rewind <- get(load('C:/Users/Mark/OneDrive - The Ohio State University/clonetools/SNP_array_data/rewind_both.rda'))
samples <- unique(rewind$SamID)
resfiles <- list.files(res.path)
reslist <- lapply(1:length(resfiles),function(i){get(load(paste(res.path, '/', resfiles[i], sep='')))})
psis <- t(sapply(1:length(resfiles),function(i){reslist[[i]]$psi}))
ks <- sapply(1:nrow(psis),function(i){length(which(psis[i,]>0))})
alt.frac <- sapply(1:length(reslist),function(j){
  temp <- rewind[rewind$SamID==samples[j],]
  notna <- which(!is.na(temp$loc.end) & !is.na(temp$loc.start))
  all.lens <- sum(as.numeric(temp$loc.end[notna]) - as.numeric(temp$loc.start[notna]))
  alt <- which(sapply(1:nrow(reslist[[j]]$A),function(k){
    unique(all.equal(unname(reslist[[j]]$A[k,1:length(which(reslist[[j]]$psi>0))]), rep(1, length(which(reslist[[j]]$psi>0))))!=TRUE) & 
      unique(all.equal(unname(reslist[[j]]$B[k,1:length(which(reslist[[j]]$psi>0))]), rep(1, length(which(reslist[[j]]$psi>0))))!=TRUE)
  }))
  lens <- unlist(sapply(1:length(alt), function(k){
    index <- which(rewind$seg.median==reslist[[j]]$filtered.data$cndata.filt$LRR[alt[k]] &
                     rewind$AvgBAF==reslist[[j]]$filtered.data$cndata.filt$BAF[alt[k]])
    rewind$loc.end[index] - rewind$loc.start[index]
  }))
  #alt.marks <- sum(cf.results[[j]]$filtered.data$cndata.filt$markers[alt])
  #alt.marks/total.markers
  sum(as.numeric(lens))/all.lens
})
psivec <- as.vector(psis)

png('C:/Users/Mark/OneDrive - The Ohio State University/clonetools/cf_images/figureS4.png',width=1000,height=700,res=100)
par(mfrow=c(2,1))
hist(ks,breaks=seq(from=.5,to=5.5,by=1), main='Frequency distribution of K',xlab='Number of subclones')
#plot(density(psivec[which(psivec>0)],from=0,to=1,bw=.01),main='Distribution of Psis',xlab='Psi',xlim=c(0,1))
plot(ks, alt.frac, main='', xlab='Number of subclones', ylab='Altered Genomic Fraction')
dev.off()


###Figure 5 new:
data.path <- 'C:/Users/Mark/OneDrive - The Ohio State University/clonetools/dat/data-CLL'
cll.files <- list.files(data.path)
cids.data <- sapply(1:length(cll.files),function(j){strsplit(strsplit(cll.files[j],split='dat-')[[1]][2],split='.rda')[[1]][1]})
clinical <- get(load('C:/Users/Mark/Documents/Lab/clinOk.rda'))
cids.clinical <- rownames(clinical)
common.cids <- intersect(cids.clinical,cids.data)
eligible <- get(load('C:/Users/Mark/OneDrive - The Ohio State University/clinical/eligible.rda'))
clinical$ttp.eligible <- eligible
clinical <- clinical[sapply(1:length(common.cids),function(j){which(rownames(clinical)==common.cids[j])}),]
ttp.eligible <- which(clinical$ttp.eligible)
res.files <- list.files(res.path)
cf.results <- lapply(1:length(res.files), function(j){get(load(paste(res.path, '/', res.files[j], sep='')))})
#cf.results <- lapply(1:length(common.cids),function(j){cf.results[[which(samples==common.cids[j])]]})
ks <- sapply(1:length(cf.results),function(j){length(which(cf.results[[j]]$psi>0))})
ok <- sapply(1:length(common.cids),function(i){which(cids.data==common.cids[i])})
ks <- ks[ok]
k.status <- rep(0,length(ks))
k.status[ks>1] <- 1
columns <- colnames(clinical)

TimeDiag2SigTreat <- clinical$TimeDiagnosis2SigTreat
TTP <- clinical[,which(columns=='TTP.time')]
TTP.numericStatus <- clinical[,which(columns=='TTP.numericStatus')]
osat.time <- clinical$TxDSS.time[ttp.eligible]
osas.time <- clinical$OSAfterSample
numericVitalStatus <- clinical$NumericVitalStatus
numericVitalStatus.elig <- numericVitalStatus[ttp.eligible]
NumericSigTreatment <- clinical$NumericSigTreatment

colors <- c('red', 'blue', 'purple', 'darkorange')
Multiclonal <- k.status
Multiclonal.elig <- k.status[ttp.eligible]
Complex <- rep(NA, nrow(clinical))
Complex[which(clinical$CatCyto=='Simple')] <- 0
Complex[which(clinical$CatCyto=='Complex')] <- 1
Complex.elig <- Complex[ttp.eligible]
df.main <- data.frame('osas.time'=osas.time,'numericVitalStatus'=numericVitalStatus,'Multiclonal'=Multiclonal,
    'Complex'=Complex,'TimeDiag2SigTreat'=TimeDiag2SigTreat,'NumericSigTreatment'=NumericSigTreatment)
df.elig <- data.frame('ttp.time'=TTP[ttp.eligible],'Multiclonal'=Multiclonal.elig,'Complex'=Complex.elig,
    'numericVitalStatus'=numericVitalStatus.elig,'osat.time'=osat.time,'TTP.status'=TTP.numericStatus[ttp.eligible])

osas <- prodlim(Surv(osas.time, numericVitalStatus) ~ Multiclonal,data=df.main)
diff <- survdiff(Surv(osas.time, numericVitalStatus) ~ Multiclonal,data=df.main)
p_osas <- round(pchisq(diff$chisq, length(diff$n)-1, lower.tail = FALSE), digits=4)
osas.karyo <- prodlim(Surv(osas.time, numericVitalStatus) ~ Multiclonal + Complex,data=df.main)
diff <- survdiff(Surv(osas.time, numericVitalStatus) ~ Multiclonal + Complex,data=df.main)
p_osas.karyo <- round(pchisq(diff$chisq, length(diff$n)-1, lower.tail = FALSE), digits=4)

t2st <- prodlim(Surv(TimeDiag2SigTreat, NumericSigTreatment) ~ Multiclonal,data=df.main)
diff <- survdiff(Surv(TimeDiag2SigTreat, NumericSigTreatment) ~ Multiclonal,data=df.main)
p_t2st <- round(pchisq(diff$chisq, length(diff$n)-1, lower.tail = FALSE), digits=4)
t2st.karyo <- prodlim(Surv(TimeDiag2SigTreat, NumericSigTreatment) ~ Multiclonal + Complex,data=df.main)
diff <- survdiff(Surv(TimeDiag2SigTreat, NumericSigTreatment) ~ Multiclonal + Complex,data=df.main)
p_t2st.karyo <- round(pchisq(diff$chisq, length(diff$n)-1, lower.tail = FALSE), digits=4)

osat <- prodlim(Surv(osat.time, numericVitalStatus) ~ Multiclonal,df.elig)
diff <- survdiff(Surv(osat.time, numericVitalStatus) ~ Multiclonal,df.elig)
p_osat <- round(pchisq(diff$chisq, length(diff$n)-1, lower.tail = FALSE), digits=4)
osat.karyo <- prodlim(Surv(osat.time, numericVitalStatus) ~ Multiclonal + Complex,df.elig)
diff <- survdiff(Surv(osat.time, numericVitalStatus) ~ Multiclonal + Complex,df.elig)
p_osat.karyo <- round(pchisq(diff$chisq, length(diff$n)-1, lower.tail = FALSE), digits=4)

ttp <- prodlim(Surv(ttp.time, TTP.status) ~ Multiclonal,df.elig)
diff <- survdiff(Surv(ttp.time, TTP.status) ~ Multiclonal,df.elig)
p_ttp <- round(pchisq(diff$chisq, length(diff$n)-1, lower.tail = FALSE), digits=4)
ttp.karyo <- prodlim(Surv(ttp.time, TTP.status) ~ Multiclonal + Complex,df.elig)
diff <- survdiff(Surv(ttp.time, TTP.status) ~ Multiclonal + Complex,df.elig)
p_ttp.karyo <- round(pchisq(diff$chisq, length(diff$n)-1, lower.tail = FALSE), digits=4)

png('C:/Users/Mark/OneDrive - The Ohio State University/clonetools/cf_images/figure5_new.png',width=1600,height=2600,res=200)
par(mfrow=c(2,1))
plot(osas)
title(paste('A) Overall Survival after Sample; P = ',p_osas,sep=''))
plot(osat)
title(paste('B) Overall Survival after Treatment; P = ',p_osat,sep=''))
dev.off()

png('C:/Users/Mark/OneDrive - The Ohio State University/clonetools/cf_images/figureS5.png',width=1200,height=1600,res=140)
par(mfrow=c(2,1))
plot(t2st)
title(paste('C) Time to Treatment; P = ',p_t2st,sep=''))
plot(ttp)
title(paste('D) Time to Progression; P = ',p_ttp,sep=''))
dev.off()

#Figure S6
png('C:/Users/Mark/OneDrive - The Ohio State University/clonetools/cf_images/figureS6.png',width=1800,height=1800,res=140)
par(mfrow=c(2,2))
plot(osas.karyo)
title(paste('A) Overall Survival after Sample; P = ',p_osas.karyo,sep=''))
plot(osat.karyo)
title(paste('B) Overall Survival after Treatment; P = ',p_osat.karyo,sep=''))
plot(t2st.karyo)
title(paste('C) Time to Treatment; P = ',p_t2st.karyo,sep=''))
plot(ttp.karyo)
title(paste('D) Time to Progression; P = ',p_ttp.karyo,sep=''))
dev.off()
