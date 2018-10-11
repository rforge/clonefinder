#library(CloneFinder)

tumorGen <- function(...) {
  tumor <- Tumor(...)
  as(tumor, "list")
}

snpDataGen <- function(tumor, snps.cgh=600000, sigma0.lrr=.01, sigma0.baf=.01, density.sigma=.1){
  psi <- tumor$psi
  cn.clones <- lapply(1:length(tumor$clones), function(i){tumor$clones[[i]]$cn})
  chr <- cn.clones[[1]]$chr
  eta <- data.frame('A'=Reduce('+', lapply(1:length(which(psi>0)), function(j){psi[j]*cn.clones[[j]]$A})), 
                    'B'=Reduce('+', lapply(1:length(which(psi>0)), function(j){psi[j]*cn.clones[[j]]$B})))
  eta.total <- eta$A + eta$B
  mu.baf <- eta$B/eta.total
  mu.lrr <- log10(eta.total/2)
  snp.dens <- rbeta2(1:nrow(eta), .5, density.sigma)
  markers <- as.vector(rmultinom(1, size=snps.cgh, prob=snp.dens))
  sigmas.lrr <- sigma0.lrr/(markers)^.5
  sigmas.baf <- sigma0.baf/(markers)^.5
  lrr <- rnorm(length(markers), mean=mu.lrr, sd=sigmas.lrr)
  baf <- rnorm(length(markers), mean=mu.baf, sd=sigmas.baf)
  baf[baf<0] <- -baf[baf<0]
  baf[baf>1] <- 1/baf[baf>1]
  total.backcomp <- 2*(10^(lrr))
  X <- total.backcomp*(baf)
  Y <- total.backcomp*(1 - baf)
  segdf <- data.frame('chr'=chr, seg=1:nrow(eta), 'LRR'=lrr, 'BAF'=baf, 'X'=X, 'Y'=Y, 'markers'=markers)
  segdf
}

#snp.info is an option for an already made snp df (like from snpDataGen) to generate
#read counts from those in addition to generating new ones.
seqDataGen <- function(tumor, snps.seq=1000000, density.sigma, mu, sigma.reads){
  psi <- tumor$psi
  cn.clones <- lapply(1:length(tumor$clones), function(i){tumor$clones[[i]]$cn})
  chrs <- cn.clones[[1]]$chr
  psi <- tumor$psi
  eta <- data.frame('A'=Reduce('+', lapply(1:length(which(psi>0)), function(j){psi[j]*cn.clones[[j]]$A})), 
                    'B'=Reduce('+', lapply(1:length(which(psi>0)), function(j){psi[j]*cn.clones[[j]]$B})))
  eta.total <- eta$A + eta$B
  snp.dens <- rbeta2(1:nrow(eta), .5, density.sigma)
  markers <- as.vector(rmultinom(1, size=snps.seq, prob=snp.dens))
  starts <- cn.clones[[1]]$start
  ends <- cn.clones[[1]]$ends
  snpdf <- matrix(NA, nrow=sum(markers), ncol=7)
  colnames(snpdf) <- c('chr', 'seg', 'mut.id', 'refCounts', 'varCounts', 'VAF', 'totalCounts')
  for(j in 1:length(markers)){
    z <- sample(1:2, 2, replace=FALSE)
    pR <- eta[j, z[1]]/(eta[j, z[1]]+eta[j, z[2]])
    totalCounts <- round(rnorm(markers[j], 2*mu, sigma.reads))
    totalCounts[totalCounts<0] <- -totalCounts[totalCounts<0]
    refCounts <- rbinom(markers[j], totalCounts, pR)
    varCounts <- totalCounts - refCounts
    chr <- rep(chrs[j], markers[j])
    VAFs <- varCounts/(totalCounts)
    seg <- rep(j, markers[j])
    for(k in 1:length(refCounts)){
      snpdf[sum(markers[1:j]) - markers[j] +k, ] <- c(chr[k], seg[k], NA, refCounts[k], varCounts[k], VAFs[k], totalCounts[k])
    }
  }
  snpdf <- as.data.frame(snpdf)
  snpdf$status <- rep('germline', nrow(snpdf))
  seq.clones <- lapply(1:length(tumor$clones), function(i){tumor$clones[[i]]$seq})
  cn.clones <- lapply(1:length(tumor$clones), function(i){tumor$clones[[i]]$cn})
  mutids <- unique(unlist(lapply(1:length(which(psi>0)), function(j){seq.clones[[j]]$mut.id})))
  mutdf <- matrix(NA, nrow=length(mutids), ncol=7)
  colnames(mutdf) <- c('chr', 'seg', 'mut.id', 'refCounts', 'varCounts', 'VAF', 'totalCounts')
  if(length(mutids)>0){
    for(j in 1:length(mutids)){
      indices <- which(sapply(1:length(which(psi>0)), function(k){mutids[j] %in% seq.clones[[k]]$mut.id}))
      indices.unmutated <- which(sapply(1:length(which(psi>0)), function(k){!mutids[j] %in% seq.clones[[k]]$mut.id}))
      coefs <- psi[indices]
      mutated <- sapply(indices, function(k){seq.clones[[k]]$mutated.copies[which(seq.clones[[k]]$mut.id==mutids[j])]})
      normal <- sapply(indices, function(k){seq.clones[[k]]$normal.copies[which(seq.clones[[k]]$mut.id==mutids[j])]})
      seg <- seq.clones[[indices[1]]]$seg[which(seq.clones[[indices[1]]]$mut.id==mutids[j])]
      if(length(indices.unmutated)>0){
        coefs.unmutated <- psi[indices.unmutated]
        eta.unmutated <- sum(psi[indices.unmutated]*sapply(indices.unmutated, function(i){
          cn.clones[[i]]$A[seg]+cn.clones[[i]]$B[seg]}))
      }else{
        eta.unmutated <- 0
      }
      eta <- c(sum(coefs*mutated), sum(coefs*normal) + eta.unmutated)
      pR <- eta[2]/(eta[2]+eta[1])
      totalCount <- round(rnorm(1, 2*mu, sigma.reads))
      totalCount[totalCount<0] <- -totalCount[totalCount<0]
      refCount <- rbinom(1, totalCount, pR)
      varCount <- totalCount - refCount
      VAF <- varCount/(varCount+refCount)
      chr <- seq.clones[[indices[1]]]$chr[which(seq.clones[[indices[1]]]$mut.id==mutids[j])]
      mutdf[j, ] <- c(chr, seg, mutids[j], refCount, varCount, VAF, totalCount)
    }
  }
  mutdf <- as.data.frame(mutdf)
  mutdf$status <- rep('somatic', nrow(mutdf))
  vardf <- rbind(snpdf, mutdf)
  vardf <- vardf[with(vardf, order(seg)), ]
  vardf
}

dataGen <- function(tumor, snps.seq, snps.cgh, mu, sigma.reads, sigma0.lrr, sigma0.baf, density.sigma){
  if(snps.cgh>0){
    cndat <- snpDataGen(tumor, snps.cgh, sigma0.lrr, sigma0.baf, density.sigma)
  }else{
    cndat <- NA
  }
  if(snps.seq){
    seqdat <- seqDataGen(tumor, snps.seq, density.sigma, mu, sigma.reads)
  }else{
    seqdat <- NA
  }
  list('cn.data'=cndat, 'seq.data'=seqdat)
}

#A plot function to visualize data and verify that the simulation is working.
plot.data <- function(tumor, data, snp=TRUE, somatic=TRUE, germline=FALSE){
  snpdata <- data$cn.data
  seqdata <- data$seq.data
  psi <- tumor$psi
  cn.clones <- lapply(1:length(tumor$clones), function(i){tumor$clones[[i]]$cn})
  psi <- tumor$psi
  markers <- snpdata$markers
  starts <- sapply(1:length(markers), function(j){sum(markers[1:j])-markers[j]})
  ends <- sapply(1:length(markers), function(j){sum(markers[1:j])})
  eta <- data.frame('A'=Reduce('+', lapply(1:length(psi), function(j){psi[j]*cn.clones[[j]]$A})), 
                    'B'=Reduce('+', lapply(1:length(psi), function(j){psi[j]*cn.clones[[j]]$B})))
  lesser <- sapply(1:nrow(eta), function(j){min(eta$A[j], eta$B[j])})
  lesser.data <- sapply(1:nrow(eta), function(j){min(snpdata$X[j], snpdata$Y[j])})
  greater <- sapply(1:nrow(eta), function(j){max(eta$A[j], eta$B[j])})
  greater.data <- sapply(1:nrow(eta), function(j){max(snpdata$X[j], snpdata$Y[j])})
  errors <- c(lesser.data - lesser, greater.data - greater)
  par(mfrow=c(3, 1))
  plot(0, type='n', xlim=c(0, sum(markers)), ylim=c(0, 5), main='Number of Lesser Allele Copies', 
       xlab='Marker Index', ylab='Expected Copy Number')
  sapply(1:nrow(eta), function(j){
    segments(x0=starts[j], x1=ends[j], y0=lesser[j], y1=lesser[j], col='black')
    segments(x0=starts[j], x1=ends[j], y0=lesser.data[j], y1=lesser.data[j], col='red')
  })
  plot(0, type='n', xlim=c(0, sum(markers)), ylim=c(0, 5), main='Number of Greater Allele Copies', 
       xlab='Marker Index', ylab='Expected Copy Number')
  sapply(1:nrow(eta), function(j){
    segments(x0=starts[j], x1=ends[j], y0=greater[j], y1=greater[j], col='black')
    segments(x0=starts[j], x1=ends[j], y0=greater.data[j], y1=greater.data[j], col='blue')
  })
  plot(density(errors), xlab='Error Index', main='Density of Segment Errors')
  mutations <- seqdata[seqdata$status=='somatic', ]
}
