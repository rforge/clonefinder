library(gtools)
library(CloneFinder)

estBetaParams <- function(mu, s2) {
  temp <- mu*(1-mu)/s2 - 1
  alpha <- mu*temp
  beta <- (1 - mu)*temp
  abs(c(alpha = alpha, beta = beta))
}

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

##Loading CLL data
cll <- read.table('D:\\SNP Array Data\\Data_CLL.txt', header=T)
cids <- unique(cll$SamID)

#'Analyze' computes the most likely sets of psi values and Z matrices for 
#each value of 'k' (number of subclones)
analyze <- function(n, sigma2, data){
  cid <- cids[n]
  #sub <- cll[cll$SamID==cid,]
  sub <- as.data.frame(data[which(data[,1]==cid),])
  markers <- sub$num.mark
  dataset <- data.frame('lrr'=sub$seg.median, 'baf'=sub$AvgBAF)
  colnames(dataset) <- c("x", "y")
  dataset <- as.data.frame(dataset)
  compModel <- CompartmentModel(markers, xy, sigma2)  
  # good guess at psi-vector
  pcm <- PrefitCloneModel(dataset, compModel)
  upd <- updatePhiVectors(pcm, compModel)
  mfun <- function(m){
    set.seed(n+3)
    # temp <- upd@phipick
    # for(i in 1:nrow(temp)){
    #   inds <- which(temp[i,] %in% sort(temp[i,], decreasing=TRUE)[1:n])
    #   exc <- which(!1:5 %in% inds)
    #   vals <- temp[i,inds]
    #   temp[i,inds] <- vals/sum(vals)
    #   temp[i,exc] <- 0
    # }
    #upd@phipick <- temp
    estpsi <- guessPsi(upd, m) # 3 is number of clones we are trying to fit
    #final <- runEMalg(estpsi, dataset, compModel, upd@phipick)
    #not sure why upd@phipick is there?
    final <- runEMalg(estpsi, dataset, compModel)
    list(final, upd@phipick)
  }
  results <- sapply(1:6, function(z){mfun(z)[[1]]})
  output <- list(results, dataset, lapply(1:6, function(z){mfun(z)[[2]]}))
  #  assign(paste("sim.", n, sep=""),output)
  #save(list=paste("sim.", n, sep=""), file=paste("D:\\Mark\\Simulations5\\sim.", n, ".rda", sep=""))
  #  save(output, file=paste("D:\\Mark\\Simulations5\\sim.", n, ".rda", sep=""))
  output
}

temp <- analyze(1, sigma2=1, cll, cids)
results <- lapply(1:length(cids), analyze, sigma2=1, data=cll)
#results2 <- results
#save(results2, file='G:\\results-cll2.rda')


###Looking at results:
#load('D:\\results-cll2.rda')
#results <- results2

###Now, we apply the prior on k, the number of subclones, and take the value
#resulting from the plurality of thetas considered:
thetas <- seq(from=1, to=500, length=500)
AuerGervini <- function(i){
  loglikes <- unlist(results[[i]][[1]][3,])
  thetafun <- function(theta){
    which.max(log(dexp(1:6, theta)) + loglikes)
  }
  picks <- sapply(thetas, thetafun)
  bin <- sapply(2:length(picks), function(j){picks[j]!=picks[j-1]})
  breaks <- which(bin) + 1
  plot(1, type='n',xlim=c(thetas[1], thetas[length(thetas)]), ylim=c(0, 7), xlab='theta', 
       ylab='K', main=cids[i])
  segments(x0=c(0, breaks-1), x1=c(breaks, length(picks)), y0=unique(picks), 
           y1=unique(picks))
  Mode(picks)  
}
Ks <- sapply(1:length(results), AuerGervini)
length(which(Ks==1))

#Plotting the LRR-BAF graphs:
plotfun <- function(i){
 sizes <- cll[cll$SamID==cids[i],]$num.mark
 dmat <- sapply(1:5, function(j){((results[[i]][[2]][,1] - xy$x[j])^2 + 
                         (results[[i]][[2]][,2] - xy$y[j])^2)^.5})
 options <- c('gray47', 'purple', 'darkorange', 'seagreen3', 'seagreen4')
 colors <- options[sapply(1:nrow(dmat), function(k){which.min(dmat[k,])})]
 pureColors <- c('darkorange', 'purple', 'gray47', 'seagreen3', 'seagreen4')
 plot(results[[i]][[2]][,1], results[[i]][[2]][,2], pch=16, cex=(sizes^.333)/10, 
      col=colors,
      xlim=c(-.5,.5), ylim=c(0, .5), xlab='LRR', ylab='BAF', main=cids[i])
 points(x=log10(c(1, 2, 2, 3, 4)/2), y=c(0, 0, .5, 1/3, 1/4), pch=1, cex=8, col=
         pureColors)
}

plotfun(8)
#monoclonal samples (maybe) : 3, 4, 9, 11, 42, 57, 67, 74, 76

#psi <- results[[i]][[1]][1,3]

####Looking at plots:
AuerGervini(11)
plotfun(11)

##Figures:
hist(Ks, breaks=5, xlab='K = number of subclones', main='Histogram of number of clones per patient')
plotfun(1)

##Alterations and clonality: do alterations in particular chrom arms 
##correlate with clonal heterogeneity?
#Chlens is a file with chromosome lengths
chlens <- read.table('F:\\chlens.txt')
#cb <- read.table('C:\\Users\\Mark\\Documents\\Temp\\cytoBand.txt')
cb <- read.table('F:\\cytoBand.txt')
cens <- cb[cb$V5=='acen',]
chrs <- as.character(cens$V1)
chr.n <- as.numeric(sapply(1:nrow(cens), function(i){strsplit(chrs[i], split='r')[[1]][2]}))
cens$chr <- chr.n
cens <- cens[with(cens, order(chr)),]
cens <- cens[1:44,]
multi <- which(Ks>1)
mono <- which(Ks==1)

foo <- function(chunk, thresh1, thresh2){
  chr <- ceiling(chunk/2)
  arm <- chunk - chr*2 + 1
  sub <- cll[cll$chrom==chr & cll$loc.start>=cens$V2[chunk] & cll$loc.end <=chlens$V2[chr] & cll$num.mark>=thresh1 & 
               abs(cll$seg.median)>=thresh2,]
  cids.sub <- unique(sub$SamID)
  cids.sub
}
abers <- lapply(1:44, foo, thresh1=100, thresh2=.05)

abfrac.mono <- sapply(1:length(abers), function(i){length(intersect(abers[[i]], unique(cll$SamID)[mono]))/length(mono)})
abfrac.multi <- sapply(1:length(abers), function(i){length(intersect(abers[[i]], unique(cll$SamID)[multi]))/length(multi)})
ch <- as.vector(sapply(1:22, function(i){rep(i, 2)}))
arm <- as.vector(t(sapply(c('p','q'), function(i){rep(i, 222)})))
ch.arm <- sapply(1:44, function(i){paste(ch[i], arm[i], sep='')})
df <- data.frame('Chr.arm'=ch.arm, 'Single'=abfrac.mono, 'Multiple'=abfrac.multi)
write.csv(df, file='F:\\Manuscripts\\Clonefinder-CNV\\abfreq.csv')


###Clinical: Regression of survival on clone size (with a given alteration)
clinical <- read.csv('E:\\CLL\\Clinical\\SNPClinicalCall_revised.csv')
diagnosis <- clinical$Date.of.diag0sis
treatment <- clinical$Date.of.1st.treatment

clinical <- read.csv('E:\\CLL\\Clinical\\SNPClinicalCall_revised.csv')
arms <- rep(c('p', 'q'), 22)
cens$arms <- sapply(1:nrow(cens), function(i){paste(cens$chr[i], arms[i], sep='')})

###Attempt at assessing aberration frequency by chromosome arm:
g <- function(arm){
  ind <- which(cens$arms==arm)
  chrom.arm <- ceiling(ind/2)
  start.arm <- cens$V2[ind]
  end.arm <- chlens$V2[chrom.arm]
  psis.arm <- sapply(1:length(cids), function(i){
    segs <- cll[cll$SamID==cids[i],]
    segs.arm <- intersect(intersect(which(segs$chrom==chrom.arm), which(segs$loc.start>=start.arm)), 
              which(segs$loc.end<=end.arm))
    Zs <- results2[[i]][[1]][2,Ks[i]]$Zmats
    inner <- function(j){
      z.sub <- Zs[segs.arm,,j]
      if(class(z.sub)=='matrix'){
        any(z.sub[,1]!=1)
      }else{
        z.sub[1]!=1
      }
    }
    altered <- which(sapply(1:dim(Zs)[3], inner))
    psi.altered <- unlist(results2[[i]][[1]][1,Ks[i]])[altered]
    if(length(altered)>0){
      output <- max(psi.altered)
    }else{
      output <- 0
    }
    output
  })

  ttt <- sapply(1:length(treatment), function(i){
    tr <- treatment[i]
    di <- diagnosis[i]
    if(is.na(tr)){
      output <- 0
    }else{
      end <- unlist(strsplit(as.character(tr), split='/'))
      end <- end[c(3, 1, 2)]
      start <- unlist(strsplit(as.character(di), split='/'))
      start <- start[c(3, 1, 2)]
      diff <- as.Date(paste(end, collapse='-')) - as.Date(paste(start, collapse='-'))
      output <- as.numeric(diff)
    }
    output
  })

  x <- sapply(1:nrow(clinical), function(i){psis.arm[which(cids %in% clinical[i,1])]})
  y <- ttt
  include <- which(y!=0)
  y <- y[include]
  x <- x[include]
  linmod <- lm(y~x)
  intercept <- summary(linmod)[4]$coefficients[1,1]
  slope <- round(as.numeric(unlist(summary(linmod)[4]$coefficients[2,1])), digits=3)
  r2 <- round(as.numeric(unlist(summary(linmod)[8])), digits=3)
  maintext <- paste(arm, '; R^2=', r2, '; coef=', slope, sep='')
  plot(x, y, main=maintext, xlab='psi', ylab='time to treatment (days)')
  abline(a=intercept, b=slope, col='red')
  data.frame('SamID'=cids[include], 'Freq'=x)
}

temp <- g('13q')

###Generating simulations:
psiList <- unlist(lapply(1:5, function(i){lapply(1:30, function(j){
  as.vector(rdirichlet(1, rep(1, i)))})}), recursive=FALSE)
nsegs <- round(runif(length(psiList), 100, 250))
s2.lrr <- .05
s2.baf <- .01
theta1 <- .5
probs <-c(.3, .3, .3, .1)
L <- 600000
x <- log10(c(2, 2, 1, 3, 4)/2)
y <- c(.5, 0.0001, 0.0001, .25, .333)

sim <- function(i){
  npick <- round(rexp(length(psiList[[i]]), theta1))
  npick[npick==0] <- 1
  list1 <- lapply(1:length(npick), function(j){sample(1:nsegs[i], npick[j])})
  markers <- runif(nsegs[i], 1, 20)
  markers <- markers/sum(markers)
  markers <- round(markers*L)
  alpha <- t(sapply(1:length(markers),function(j){sapply(1:5,function(i){
    estBetaParams(y[i], s2.baf)[1]})/markers[j]}))
  beta <- t(sapply(1:length(markers), function(j){sapply(1:5, function(i){
    estBetaParams(y[i], s2.baf)[2]})/markers[j]}))
  A <- lapply(1:length(psiList[[i]]), function(j){
    mat <- matrix(rep(0, nsegs[i]*5), nrow=nsegs[i], ncol=5)
    mat[,1] <- rep(1, nrow(mat))
    colpick <- sample(2:5, length(list1[[j]]), prob=probs, replace=TRUE)
    for(k in list1[[j]]){
      mat[k,colpick[k]] <- 1
      mat[k,1] <- 0
    }
    baf <- rbeta(nsegs[i], alpha[1], beta[1])
    lrr <- rnorm(nsegs[i], x[1], s2.lrr)
    baf[list1[[j]]] <- rbeta(length(list1[[j]]), alpha[colpick], beta[colpick])
    lrr[list1[[j]]] <- rnorm(length(list1[[j]]), x[colpick], s2.lrr)
    baf[baf>.5] <- 1 - baf[baf>.5]
    df <- data.frame('lrr'=lrr, 'baf'=baf)
    list(mat, df)
  })
  B <- list(psiList[[i]], lapply(1:length(psiList[[i]]), function(k){A[[k]][[1]]}))
  C <- list(Reduce('+', lapply(1:length(A), 
                                function(k){A[[k]][[2]]*psiList[[i]][[k]]})), markers)
  list(B, C)
}

sims <- lapply(1:length(psiList), sim)
cids.sim <- sapply(1:length(sims), function(i){paste('sim-', i, sep='')})
stuff <- sims
simtab <- data.frame('SamID'=unlist(lapply(1:length(cids.sim),function(j){rep(cids.sim[j], 
  nrow(stuff[[j]][[1]][[2]][[1]]))})), 'seg.median'=unlist(lapply(1:length(cids.sim), 
  function(j){stuff[[j]][[2]][[1]][,1]})), 'AvgBAF'=unlist(lapply(1:length(cids.sim),
    function(j){stuff[[j]][[2]][[1]][,2]})), 'num.mark'=unlist(lapply(1:length(cids.sim), 
      function(j){stuff[[j]][[2]][[2]]}))
)
stuff <- simtab
res <- lapply(1:length(cids.sim) ,analyze, sigma2=.05, data=stuff, cids=cids.sim)

###Checking simulatiosn:
plot(density(sims[[31]]))

###Running on already made simulations:
load('E:\\Mark\\res.rda')
results <- res
cids.sim <- sapply(1:length(res), function(i){paste('sim-', i, sep='')})
cids <- cids.sim
thetas <- seq(from=1, to=10, length=500)
Ks <- sapply(1:length(results), AuerGervini)
true <- sapply(1:5, function(i){rep(i, 30)})
length(which(Ks==true))/length(Ks)
