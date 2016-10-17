source('C:\\Users\\Mark\\Dropbox\\R Files\\cll\\functions-mcmc.r')
library(gtools)
library(quantmod)
#load('D:\\sdtab.rda')
chroms <- read.table('C:\\Users\\Mark\\Dropbox\\R Files\\cll\\chroms.txt', header=T)
chrset <- chrset <- c(1:22, 'X', 'Y')
chlens <- sapply(chrset, function(i){length(which(chroms[,1]==i))})
positions <- read.table('C:\\Users\\Mark\\Dropbox\\R Files\\cll\\positions.txt', header=T)
posdf <- data.frame('chr'=chroms[,1], 'pos'=positions[,1])
posdf <- posdf[which(posdf$chr %in% chrset),]
posdf$chr.n <- sapply(1:nrow(posdf), function(i){which(chrset==posdf$chr[i])})
df <- posdf[with(posdf, order(chr.n, pos)),]

simulate <- function(psilist, theta, sigma, min, snprate){
  psi <- unlist(psilist)
  K <- length(psi)
  nseg <- round(runif(1, 200, 400))
  segchroms <- c(1:24, sample(1:24, nseg-24, replace=TRUE, prob=chlens))
  segs <- unlist(sapply(1:24, function(i){
    starts <- sort(c(1, round(runif(length(which(segchroms==i))-1, 2, chlens[i]))))
    if(length(starts)>1){
      ends <- c(starts[2:length(starts)]-1, chlens[i])
    }else{
      ends <- chlens[i]
    }
    temp <- unlist(sapply(1:length(starts), function(j){rep(paste(i, '.', j, sep=''), 
                                                            ends[j]-starts[j] + 1)}))
  }))
  df$seg <- segs
  status <- sample(c('het', 'homA', 'homB'), nrow(df), replace=TRUE, prob=c(snprate, (1-snprate)/2, (1-snprate)/2))
  df$status <- status
  segs <- unique(df$seg)
  seglens <- sapply(1:length(segs), function(i){length(which(df$seg==segs[i]))})
  base <- data.frame('A'=rep(1, length(segs)), 'B'=rep(1, length(segs)))
  altmat <- t(combn(0:3, 2))
  altmat <- rbind(altmat, cbind(altmat[,2], altmat[,1]))
  altmat <- rbind(altmat, t(matrix(c(0, 0, 1, 1, 3, 3), ncol=3, nrow=2)))
  #altprobs <- dexp(rowSums(abs(altmat-2)), theta)
  altprobs <- 1/rowSums(abs(altmat-2))^4
  clones <- list(base)
  i <- 1
  while(i <= K){
    normal <- which(base$A==1 & base$B==1)
    N <- round(runif(1, 1, 8))
    normal.big <- intersect(normal, which(seglens >= median(seglens)))
    tochange <- sample(normal.big, runif(1, 1, N))
    tochange <- c(tochange, sample(normal, N - length(tochange)))
    alts <- sample(1:nrow(altmat), tochange, replace=TRUE, prob=altprobs)
    base <- as.matrix(base)
    for(j in 1:length(tochange)){
      base[tochange[j],] <- altmat[alts[j],]
    }
    newclone <- as.data.frame(base)
    clones <- c(clones, list(newclone))
    base <- clones[[sample(1:length(clones), 1)]]
    i <- i + 1
  }
  clones <- clones[2:length(clones)]
  clonemat <- as.matrix(Reduce(cbind, clones))
  #colnames(clonemat) <- c(sapply(1:k, function(i){c(paste('A.',i , sep=''), paste('B.',i , sep=''))}))
  rownames(clonemat) <- segs
  xmat <- Reduce(cbind, lapply(1:length(psi), function(i){psi[i]*clonemat[,2*(i-1)+1]}))
  xmean <- rowSums(xmat)
  ymat <- Reduce(cbind, lapply(1:length(psi), function(i){psi[i]*clonemat[,2*(i-1)+2]}))
  ymean <- rowSums(ymat)
  mu.x <- unlist(sapply(1:length(seglens), function(i){rep(xmean[i], seglens[i])}))
  mu.y <- unlist(sapply(1:length(seglens), function(i){rep(ymean[i], seglens[i])}))
  mu.x[which(df$status=='homA')] <- mu.x[which(df$status=='homA')]+mu.y[which(df$status=='homA')]
  mu.y[which(df$status=='homA')] <- min
  mu.y[which(df$status=='homB')] <- mu.y[which(df$status=='homB')]+mu.x[which(df$status=='homB')]
  mu.x[which(df$status=='homB')] <- min
  df$X <- rlnorm(nrow(df), log(mu.x), rep(sigma, nrow(df)))
  df$Y <- rlnorm(nrow(df), log(mu.y), rep(sigma, nrow(df)))
  df$lrr <- log10((df$X+df$Y)/2)
  df$baf <- df$X/(df$X+df$Y)
  df$total <- df$X + df$Y
  df <- df[with(df, order(chr.n, pos)),]
  lens <- sapply(1:22, function(i){max(df[df$chr.n==i,]$pos)})
  #metapos <- sapply(1:q,function(i){
  #  if(df$chr.n[i]>1){
  #    output <- sum(lens[1:(df$chr.n[i]-1)]) + df$pos[i]}
  #  else{
  #    output <- df$pos[i]
  #  }
  #  output
  #})
  #df$metapos <- metapos
  rownames(df) <- NULL
  #plot(df[df$chr==2,]$pos, df[df$chr==2,]$lrr, pch='.')
  segmented <- as.data.frame(t(sapply(1:length(segs),function(i){
    colMeans(df[which(df$seg==segs[i]),6:9])})))
  #seg.baf <- sapply(1:length(segs), function(i){
  #  den <- density(df[df$seg==segs[i],]$baf, bw=.1)
  #  denx <- den$x
  #  low <- den$x[min(findValleys(den$y))]
  #  high <- den$x[max(findValleys(den$y))]
  #  baf <- df[df$seg==segs[i],]$baf
  #  subset <- baf[which(baf>low & baf<high)]
  #  subset[which(subset>.5)] <- 1 - subset[which(subset>.5)]
  #  if(length(subset)>2){
  #    newden <- density(subset)
  #    peaks <- findPeaks(newden$y)
  #    peaks2x <- newden$x[peaks]
  #    peaks2y <- newden$y[peaks]
  #    output <- peaks2x[which.max(peaks2y)]
  #  }else{
  #    output <- median(baf)
  #    if(output>.5){
  #      output <- 1 - .5
  #    }
  #  }
  #  output
  #})
  seg.baf <- sapply(1:length(segs),function(i){
    med <- median(df[df$seg==segs[i],]$baf)
    min(med, 1-med)
  })
  segmented$markers <- seglens
  segmented$total <- segmented$X + segmented$Y
  segmented$seg <- segs
  segmented$baf <- seg.baf
  #X.new <- sapply(1:nrow(segmented), function(i){
  #  if(segmented$X[i] < segmented$Y[i]){
  #    out <- segmented$total[i]*segmented$baf[i]
  #  }else{
  #    out <- segmented$total[i]*(1 - segmented$baf[i])
  #  }
  #  out
  #})
  Y.new <- segmented$total*segmented$baf
  X.new <- segmented$total - Y.new
  segmented$Y <- Y.new
  segmented$X <- X.new
  
  #filtered <- filter(segmented, 1.2, .9, .002)
  #guess <- as.numeric(rownames(filtered))
  true <- which(sapply(1:nrow(clonemat), function(i){length(which(clonemat[i,]==1))==ncol(clonemat)})=='FALSE')
  #true.filt <- clonemat[guess,]
  filtered <- segmented[true,]
  ####
  condensed <- t(sapply(1:nrow(clonemat), function(i){sapply(1:length(psi),
                                                             function(j){sum(clonemat[i,((j-1)*2+1):(2*j)])})}))
  
  filter2 <- function(data, sigma, threshold, min){
    dat <- data
    logged <- log(dat$total)
    log.mu <- log(round(dat$total))
    log.mu[which(is.infinite(log.mu))] <- min
    dif <- abs(logged - log.mu)
    pvec <- 1- (pnorm(log.mu + dif, log.mu, sigma)-pnorm(log.mu - dif,
                                                         log.mu, sigma))
    #hist(pvec, breaks=100)
    dat[which(pvec <= threshold),]
  }
  #filtered <- filter2(segmented, .6, .9, .002)
  #newdf <- df
  #cols <- sapply(1:nrow(newdf),function(i){
  #  if(newdf$seg %in% filtered$seg){
  #    col <- 'red'
  #  }else{
  #    col <- 'gray'
  #  }
  #  col
  #})
  #newdf$col <- cols
  guess <- as.numeric(rownames(filtered))
  true <- which(sapply(1:nrow(condensed), function(i){length(unique(condensed[i,]))>1}))
  true.filt <- clonemat[guess,] 
  data <- segmented[true,]
  expected <- data.frame('X'=xmean, 'Y'=ymean)
  expected2 <- as.data.frame(t(sapply(1:nrow(expected),function(i){sort(expected[i,], decreasing=TRUE)})))
  colnames(expected2) <- c('X', 'Y')
  truexy.unordered <- expected[true,]
  truexy <- expected2[true,]
  
  truecn <- clonemat[true,]
  truecn <- truecn[,c(seq(from=1, to=ncol(truecn), by=2), seq(from=2, to=ncol(truecn), by=2))]
  colnames(truecn) <- c(sapply(1:(ncol(truecn)/2),function(i){paste('A',i,sep='')}), 
                        sapply(1:(ncol(truecn)/2),function(i){paste('B',i,sep='')}))
  truemat <- truecn
  truemat <- t(sapply(1:nrow(truemat),function(i){
    K <- ncol(truemat)/2
    order.ab <- unname(c(which.max(unlist(truexy.unordered[i,])), which.min(unlist(truexy.unordered[i,]))))
    truemat[i,c(sapply(order.ab, function(j){((j-1)*K+1):(j*K)}))]
  }))
  list(segmented, filtered, data, clonemat, truemat, truexy, psilist, df)
}

#test1 <- simulate(psilist=list(.8, .2), theta=1, sigma=.2, min=.005, snprate=.6)
#filt <- test1[[2]]
#truth <- test1[[5]]

#dat <- c(filt[,1], filt[,2])
#AB <- rbind(truth[,1:2], truth[,3:4])


###Simulating just the altered segments:
#s <- 5
#Max <- 4
#psilist <- c(.3, .7)
#theta <- .5
#sd <- .1

simulate2 <- function(s, Max, psilist, theta, sd){
  psi <- sort(unlist(psilist))
  markers <- round(rlnorm(s, log(500), .5))
  markers <- rep(markers, 2)
  markers[markers < 50] <- 50
  K <- length(psi)
  pool <- c(0:Max)
  probs <- dgeom(abs(pool - 1), prob=theta)
  vec <- sample(pool, 2*K*s, prob=probs, replace=TRUE)
  mat <- matrix(vec, ncol=K, nrow=2*s)
  normal <- which(sapply(1:nrow(mat), function(i){length(which(mat[i,]==1))==ncol(mat)}))
  if(length(normal>0)){
    newvals <- sample(c(0, 2:Max), length(normal), replace=TRUE, prob=probs[c(1, 3:(Max+1))])
    rows <- sample(1:ncol(mat), length(normal), replace=TRUE)
    for(j in 1:length(normal)){
      mat[normal[j],rows[j]] <- newvals[j]
    }
  }
  expected <- rowSums(mat %*% psi)
  error <- sapply(1:length(expected), function(i){rnorm(1, 0, sd/(markers[i]^.5))})
  neg <- which(sapply(1:length(error), function(i){error[i] <= -expected[i]}))
  error[neg] <- -error[neg]
  dat <- expected + error
  colnames(mat) <- c(1:K)
  df <- as.data.frame(mat)
  df$Expected <- expected
  df$Data <- dat
  list(df, psi)
}

