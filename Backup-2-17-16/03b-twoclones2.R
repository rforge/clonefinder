### R code from vignette source '03b-twoclones.Rnw'

###################################################
### code chunk number 1: setOptions
###################################################
setwd('F:\\H&N')
set.seed(2692476)
options(width=88)
options(SweaveHooks = list(fig = function() par(bg='white')))


###################################################
### code chunk number 2: makeFiguresDirectory
###################################################
if (!file.exists("Figures")) {
  dir.create("Figures")
}


###################################################
### code chunk number 3: load
###################################################
load("hn2.Rda")
ls()


###################################################
### code chunk number 4: 03b-twoclones.Rnw:172-179
###################################################
goo <- function(cn, theta) exp(-theta*(abs(cn-2)))
opar <- par(mfrow=c(2,1))
pts <- barplot(goo(0:15, 0.01), yaxt='n', main="Prior When Theta = 0.01")
mtext(0:15, side=1, at=pts, line=1)
pts <- barplot(goo(0:15, 0.9), yaxt='n', main="Prior When Theta = 0.9")
mtext(0:15, side=1, at=pts, line=1)
par(opar)


###################################################
### code chunk number 5: cnMAP
###################################################
source("00fracture.R")
cnMAP


###################################################
### code chunk number 6: LL
###################################################
LL


###################################################
### code chunk number 7: LLfun
###################################################
LLfun


###################################################
### code chunk number 8: fracture
###################################################
Fracture


###################################################
### code chunk number 9: fits
###################################################
gl <- GenomicLocation(locns)
pg <- setupGrid()
f <- "newfits.Rda"
if (file.exists(f)) {
  load(f)
} else {
  fracs <- list()
  for (j in 1:ncol(imputedCNV2)) {
    fracs[[j]] <- Fracture(imputedCNV2[,j], gl, pg, 
                           theta=0.9, sigma=0.04, LABEL=j)
  }
  save(fracs, file=f)
}
rm(f)
xnew <- sapply(fracs, function(f) f@fit$par[1])
psi <- sapply(fracs, function(f) f@fit$par[2])


###################################################
### code chunk number 10: makeFigs (eval = FALSE)
###################################################
## if (!file.exists("ColoredFigs")) dir.create("ColoredFigs")
## for (j in 1:length(fracs)) {
##   frac <- fracs[[j]]
##   png(paste("ColoredFigs/fig", j, ".png", sep=''), 
##       width=1800, height=1200, bg="white")
##   tag <- paste(colnames(imputedCNV2)[j], "  (", 
##                sampleInfo$platform[j], ", ", 
##                sampleInfo$gender[j], ", HPV ",
##                sampleInfo$hpv.status[j], ")", sep='')
##   opar <- par(mfrow=c(2,1))
##   hist(frac, main=colnames(imputedCNV2)[j])
##   plot(frac, main=tag)
##   par(opar)
##   dev.off()
## }


###################################################
### code chunk number 11: 03b-twoclones.Rnw:262-266
###################################################
opar <- par(mai=c(0.7, 0.1, 0.7, 0.1))
barplot(rep(1,11), names.arg=0:10, col=cncolor, yaxt='n')
title("Color Code for Copy Number")
par(opar)


###################################################
### code chunk number 12: 03b-twoclones.Rnw:280-293 (eval = FALSE)
###################################################
## par(bg="white")
## plot(xnew, pch=16, col=pcol, ylab="x0")
## abline(h=2)
## identify(1:132, xnew, colnames(imputedCNV2))
## dev.copy(png, "Figures/x0.png", width=800, height=600, bg="white")
## dev.off()
## 
## par(bg="white")
## plot(psi, pch=16, col=pcol)
## abline(h=2)
## identify(1:132, psi, colnames(imputedCNV2))
## dev.copy(png, "Figures/psi.png", width=800, height=600, bg="white")
## dev.off()


###################################################
### code chunk number 13: tilt
###################################################
source("01tilt.R")
DOT

dot

tilt

showtilt


###################################################
### code chunk number 14: jeffreys
###################################################
library(mc2d)
library(fields)
set.seed(398699)
jeffreys <- rdirichlet(2000, rep(1/2, 3))
uniform <- rdirichlet(2000, rep(1, 3))


###################################################
### code chunk number 15: 03b-twoclones.Rnw:370-374
###################################################
opar <- par(mfrow=c(2,1), mai=c(0.8, 0.8, 0.1, 0.1))
showtilt(uniform)
showtilt(jeffreys)
par(opar)


###################################################
### code chunk number 16: folded
###################################################
m <- apply(jeffreys[, 1:2], 1, which.max)
folded <- jeffreys
folded[m==1,1] <- jeffreys[m==1,2]
folded[m==1,2] <- jeffreys[m==1,1]


###################################################
### code chunk number 17: 03b-twoclones.Rnw:398-399
###################################################
showtilt(folded, xlim=c(-1,1), ylim=c(-0.5, 1))


###################################################
### code chunk number 18: hexgrid
###################################################
v <- function(angle) c(cos(angle), sin(angle))
v(pi/3)

gr1 <- matrix(rep(0:30, each=31), ncol=1) %*% matrix(c(1,0), nrow=1)
gr2 <- matrix(rep(0:30, times=31), ncol=1) %*% matrix(v(pi/3), nrow=1)
colnames(gr1) <- colnames(gr2) <- c("x", "y")
gr <- gr1+gr2
gr <- gr[gr[,'x'] < 15.3,]
gr <- gr/(10*sqrt(3))
gr <- sweep(gr, 2, c(sqrt(3)/2, 1/2), "-")
# squeeze in from the edges...
eps <- 0.001
ding <- (1-2*eps/sqrt(3))*gr[,'x']
dong <- (1 - 4/3*eps)*gr[,'y'] + eps/3
gr[,'x'] <- ding
gr[,'y'] <- dong


###################################################
### code chunk number 19: 03b-twoclones.Rnw:438-440
###################################################
showtilt(folded, xlim=c(-1,1), ylim=c(-0.5, 1))
points(gr, pch=16)


###################################################
### code chunk number 20: invert.tilt
###################################################
# still need to invert the "tilt" mapping.
verse

VERSE
twopsis <- t(apply(gr, 1, VERSE))
colnames(twopsis) <- paste("psi", 1:3, sep='')
dink <- ddirichlet(twopsis, rep(1/2, 3))
krig <- Krig(gr, log(dink))
surface(krig)


###################################################
### code chunk number 21: twowayPriors
###################################################
library(gplots)
# need a better prior on PAIRS of tumor copy numbers
setClass("CNPrior", representation=list(
  cn1 = "numeric",
  cn2 = "numeric",
  logprior = "numeric",
  pmat = "matrix",
  cn = "numeric",
  weight = "numeric"
  ))

CNPrior <- function(theta, wt=4) {
  cn <- 0:9
  cn1 <- rep(cn, each=10)
  cn2 <- rep(cn, times=10)
  logprior <- - theta*( abs(cn1-2) + abs(cn2-2) + wt*abs(cn1-cn2) )
  m <- 1.2*min(logprior)
  logprior[cn1 > 2 & cn2 < 2] <- m
  logprior[cn1 < 2 & cn2 > 2] <- m
  normalizer <- log(sum(exp(logprior)))
  logprior <- logprior - normalizer
  pmat <- matrix(logprior, ncol=10)
  new("CNPrior", cn1=cn1, cn2=cn2, logprior=logprior, pmat=pmat, 
      cn=cn, weight=wt)
}

if (!isGeneric("image"))
  setGeneric("image",
             function(x, ...) standardGeneric("image"))

setMethod('image', signature(x='CNPrior'), function(x, ...) {
  CN1 <- CN2 <- x@cn
  image(CN1, CN2, x@pmat, col=rich.colors(32))
  invisible(x)
})
cnp <- CNPrior(0.9, 30)


###################################################
### code chunk number 22: 03b-twoclones.Rnw:503-504
###################################################
image(cnp, col=rich.colors(32))


###################################################
### code chunk number 23: cnMAP2
###################################################
cnMAP2 <- function(x, x0, psis, cnp, sigma) {
  mu <- x0 + psis[1] * (cnp@cn1 - 2) + psis[2] * (cnp@cn2 - 2)
  post <- dnorm(x, mu, sigma, log=TRUE) + cnp@logprior
  index <- which.max(post)
  c(ll=post[index], cn1=cnp@cn1[index], cn2=cnp@cn2[index])
}


###################################################
### code chunk number 24: LL2
###################################################
LL2 <- function(datavector, x0, psis, cnp, sigma) {
  res <- sapply(datavector, cnMAP2, x0=x0, psis=psis, cnp=cnp, sigma=sigma)
  genell <- res["ll",]
  sum(genell)
}


###################################################
### code chunk number 25: moreColors
###################################################
moreColors <- matrix("black", nrow=6, ncol=6)
diag(moreColors) <- c("red", "orange", "gray", cncolor[4], "blue", "purple")
moreColors[1:3, 1] <- moreColors[1, 1:3] <- 
                      colorRampPalette(c("red", "white"))(5)[1:3]
moreColors[2:3, 2] <- moreColors[2, 2:3] <- 
                      colorRampPalette(c("orange", "white"))(5)[c(1,3)]
moreColors[3:4, 4] <- moreColors[4, 3:4] <- 
                      colorRampPalette(c("white", cncolor[4]))(5)[c(3, 5)]
moreColors[3:5, 5] <- moreColors[5, 3:5] <- 
                      colorRampPalette(c("white", "blue"))(5)[3:5]
moreColors[3:6, 6] <- moreColors[6, 3:6] <- 
                      colorRampPalette(c("white", "purple"))(7)[4:7]
image(0:5, 0:5, matrix(1:36, ncol=6), col=moreColors)


###################################################
### code chunk number 26: newrefits (eval = FALSE)
###################################################
source('mixtures2.R')
mixtureList <- list(mixfun(list('GS1445','GS1461'),list(.5,.5),),
                  mixfun(list('GS1138','GS1096'),list(.5,.5),),
                  mixfun(list('GS1059','GS1115'),list(.5,.5),))
mixtures <- cbind(mixtureList[[1]][[2]],mixtureList[[2]][[2]],
                  mixtureList[[3]][[2]])
rownames(mixtures) <- rownames(imputedCNV2)
mixtures <- as.data.frame(mixtures)
truepsis <- rbind(mixtureList[[1]][[3]],mixtureList[[2]][[3]],
                  mixtureList[[3]][[3]])

verbose <- TRUE
mysig <- 0.02
f <- "newrefits.Rda"
if (file.exists(f)) {
 load(f)
} else {
  refits <- list()
   for (j in 1:ncol(mixtures)) {
     x0 <- xnew[j]
     xg <- mixtures[!(locns$Chr %in% c("chrX", "chrY")),j]
     xg <- xg[xg < 10] # greater than 10 copies is definitely an outlier
     loglikes <- rep(NA, nrow(twopsis))
     for (i in 1:nrow(twopsis)) {
       if (i%%10==0 & verbose) cat(i, "\n", file=stderr())
       loglikes[i] <- LL2(xg, x0, twopsis[i,], cnp, sigma=mysig)
     }
     dex <- which.max(loglikes)
     round(psis <- twopsis[dex,], 3)
     obj <- Krig(gr, loglikes)
     refits[[j]] <- list(obj=obj, loglikes=loglikes, psis=psis, x0=x0)
   }
   save(refits, file=f)
 }
rm(f) 

plot(density(mixtures[,1]),xlim=c(1,3))
med <- median(mixtures[,1])
abline(v=med)
abline(v=c(med+truepsis[1,],med-truepsis[1,]),col='red')
abline(v=c(med+refits[[1]]$psis[2:3],med-refits[[1]]$psis[2:3]), col='blue')
###################################################
### code chunk number 27: foo
###################################################
foo <- function(cn1, cn2) {
  if(cn1 > 5 | cn2 > 5) return ("pink")
  moreColors[1+cn1, 1+cn2]
}


###################################################
### code chunk number 28: 03b-twoclones.Rnw:580-626 (eval = FALSE)
###################################################
## 
## sap <- sapply(bap <- seq(0, 10000, by=100), function(k) {
##   which.max(loglikes + k*log(dink))
## })
## plot(sap)
## sap[1:15]
## bap[1:15]
## 
## 
## for (j in 1:ncol(imoputedCNV2)) {
##   dex <- which.max(refits[[j]]$loglikes)
##   round(psis <- twopsis[dex,], 3)
##   surface(obj)
##   points(gr[dex,1], gr[dex,2], pch=16)
##   xg <- imputedCNV2[,j]  
##   psis <- refits[[j]]$psis
##   res <- sapply(xg, cnMAP2, x0=x0, psis=psis, cnp=cnp, sigma=mysig)
##   genell <- res["ll",]
##   sum(genell)
##   tab <- table(res["cn1",], res["cn2",])
##   tab
##   foolish <- rep(NA, ncol(res))
##   for (nc in 1:ncol(res)) {
##     foolish[nc] <- foo(res["cn1", nc], res["cn2", nc])
##   }
## 
##   plot(xg, col=foolish, ylim=c(1, 3), pch=16,
##        xlab="Gene Index", ylab="Copy Number",
##        main=paste(colnames(imputedCNV2)[j], "  (", 
##                   sampleInfo$platform[j], ", ", 
##                   sampleInfo$gender[j], ", HPV ",
##                   sampleInfo$hpv.status[j], ")", sep=''))
##   abline(h=0:6, col='gray75', lty='dashed')
##   abline(v=gl@revert)
##   mtext(levels(gl@locns$Chr), side=3, at=gl@xpts)
##   for (nr in 1:min(6, nrow(tab))) {
##     for (nc in 1:min(6, ncol(tab))) {
##       if (tab[nr, nc] > 50) {
##         loc <- x0 + (nr-3)*psis[1] + (nc-3)*psis[2]
##         col <- moreColors[nr, nc]
##         abline(h=loc, col=col)
##         mtext(paste("(", nr-1, ',', nc-1, ")", sep=''), side=4, at=loc, col=col, las=2)
##       }
##     }
##   }
## }


###################################################
### code chunk number 29: 03b-twoclones.Rnw:629-689 (eval = FALSE)
###################################################
## j <- 7; xo <- 1.843; psi <- 0.245;
## 
## 
## pho <- sapply(1:132, function(j) {
##   xnought[[j]] - dens[[j]]$x[uds[[j]][1]]
## })
## gag <- data.frame(x0=xnought, y0=ynought, pho=pho, x1=xnew, psi=psi)
## rownames(gag) <- colnames(imputedCNV2)
## 
## 
## plot(factor(pcol, levels=tricol), xnought, col=tricol)
## plot(factor(pcol, levels=tricol), xnew, col=tricol)
## lot(xnought, xnew, pch=16, col=pcol)
## abline(0,1)
## plot(psi, pho, pch=16, col=pcol)
## abline(0,1)
## plot(psi, pch=16, col=pcol)
## plot(factor(pcol, levels=tricol), psi, col=tricol)
## 
## 
## flagger <- sapply(1:132, function(j) {
##   xall <- dens[[j]]$x[ fps[[j]] ]
##   xud <- dens[[j]]$x[ uds[[j]] ]
##   x0 <- xnought[j]
##   wick <- which(xall %in% c(x0, xud))
##   xleft <- xall[1:max(wick)][-wick]
##   ifelse(length(xleft) > 0, abs(x0 - xleft[1]), NA)
## })
## library(quantmod)
## 
## msg <- unlist(lapply(refits, function(x) x$message))
## 
## parms <- t(as.data.frame(lapply(refits, function(x) x$par)))
## dimnames(parms) <- list(colnames(imputedCNV2), 
##                         c("x2", "psi2", "zeta"))
## gag <- data.frame(gag, parms)
## pairs(gag, pch=16, col=pcol)
## 
## pairs(gag[msg=="ERROR: ABNORMAL_TERMINATION_IN_LNSRCH",], pch=16, col=pcol)
## pairs(gag[msg=="CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTOL",], pch=16, col=pcol)
## pairs(gag[msg=="CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH",], pch=16, col=pcol)
## 
## 
## 
## plot(0,0, type='n', xlab="Log Intensity", ylab="Frequency",
##      xlim=c(-1,20), ylim=c(0, 0.3), main="Align Big Mode")
## for (j in 1:ncol(imputedCNV2)) {
##   lines(density(log.rna[,j] + 10 -modal[j]), col=j)
## }
## 
## 
## windows()
## plot(apply(adjustedRNA, 2, mean), col=pcol, pch=16)
## plot(apply(adjustedRNA, 2, sd), col=pcol, pch=16)
## 
## f <- function(j) plot(adjustedRNA[j,], pch=16, col=pcol, main=j)
## g <- function(j) plot(log.rna[j,], pch=16, col=pcol, main=j)
## h <- function(j) plot(log.cnv2[j,], adjustedRNA[j,], pch=16, col=pcol, main=j)
## h0 <- function(j) plot(log.cnv2[j,], log.rna[j,], pch=16, col=pcol, main=j)
## 


###################################################
### code chunk number 30: getwd
###################################################
getwd()


###################################################
### code chunk number 31: sessionInfo
###################################################
sessionInfo()


