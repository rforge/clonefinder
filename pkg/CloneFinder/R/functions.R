# Method of moments to reparametrize beta distribution.
estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  c('alpha' = abs(alpha), 'beta' = abs(beta)) # shouldn't negatives be an error ?! 
}

rbeta2 <- function(n, mu, sigma){
  params <- estBetaParams(mu, sigma^2)
  rbeta(n, params[1], params[2])
}

dbeta2 <- function(x, mu, sigma, log=FALSE){
  params <- estBetaParams(mu, sigma^2)
  dbeta(x, params[1], params[2], log=log)
}

psiPrior <- function(psi, alpha=.5, kmax=5, minim=.001){
  psivec <- psi
  psivec[which(psivec<minim)] <- 0
  psivec <- sort(psivec/sum(psivec), decreasing=TRUE)
  if(length(which(psivec>0))==1){
    psivec[1] <- 1 - minim
    psivec[2] <- minim
  }
  p <- log(ddirichlet(psivec[c(which(psivec>0))], rep(alpha, length(which(psivec>0)))))
  p
}


filterCN <- function(data, threshold){
  indices <- which(abs(data$X - 1) >= threshold | abs(data$Y - 1) >= threshold)
  if (length(indices) < 1) { #always return at least mone segment, even if it's normal
    indices <- 1
  }
  list('mat'=data[indices,], 'indices'=indices)
}

filterMutations <- function(mutdata, mu, threshold){
  indices <- which(abs(mutdata$refCounts - mu)>=threshold | abs(mutdata$varCounts - mu)>=threshold)
  ids <- mutdata$mut.id[indices]
  list('mat'=mutdata[indices,], 'ids'=ids)
}

seqSeg <- function(seqsnps, len, thresh){
  sigma0 <- sd(c(seqsnps$refCounts, seqsnps$varCounts))
  dens <- density(seqsnps$refCounts)
  mu <- dens$x[which.max(dens$y)]
  chrs <- unique(seqsnps$chr)
  dflist <- lapply(1:length(chrs),function(i){
    seqdata.chr <- seqsnps[seqsnps$chr == chrs[i],]
    meandiffs <- t(sapply(1:(nrow(seqdata.chr) - len), function(j) {
      mu1.ref <- mean(seqdata.chr$refCounts[((j - 1)*len+1):(j*len)])
      mu2.ref <- mean(seqdata.chr$refCounts[(j*len+1):((j+1)*len)])
      mu1.var <- mean(seqdata.chr$varCounts[((j - 1)*len+1):(j*len)])
      mu2.var <- mean(seqdata.chr$varCounts[(j*len + 1):((j + 1)*len)])
      c(abs(mu2.ref - mu1.ref), abs(mu2.var - mu1.var))
    }))
    peaks.ref <- findPeaks(meandiffs[,1], thresh=.25)
    peaks.var <- findPeaks(meandiffs[,2], thresh=.25)
    peaks <- sort(unique(c(peaks.ref, peaks.var)), decreasing=FALSE)
    if(length(peaks) > 0){
      start <- c(1, 1 + len*sapply(1:length(peaks), function(j) {peaks[j]} ))
      end <- c(start[-1] - 1, nrow(seqdata.chr))
    }else{
      start <- 1
      end <- nrow(seqdata.chr)
    }
    markers <- end - start + 1
    X <- sapply(1:length(start), function(j) {median(seqdata.chr$refCounts[start[j]:end[j]])})/mu
    Y <- sapply(1:length(start), function(j) {median(seqdata.chr$varCounts[start[j]:end[j]])})/mu
    LRR <- log10((X + Y)/2)
    BAF <- X/(X + Y)
    data.frame('chr'=rep(chrs[i],length(start)), 'start'=start, 'end'=end, 'LRR'=LRR, 
               'BAF'=BAF, 'X'=X, 'Y'=Y, 'markers'=markers)
  })
  df <- Reduce(rbind, dflist)
  df$seg <- 1:nrow(df)
  df
}
