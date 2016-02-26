setwd('D:\\Mark\\Simulations7')
load('stuff3.rda')
load('results-11-6-15.rda')
real <- stuff3[[1]]
thetas <- seq(from=.01, to=.99, length=99)
ones <- list(list(1), list(1), list(1), list(1), list(1), list(1), list(1), list(1), list(1), list(1), list(1),
             list(1), list(1), list(1), list(1), list(1), list(1), list(1), list(1), list(1), list(1), list(1)
             , list(1), list(1), list(1), list(1), list(1), list(1), list(1), list(1), list(1), list(1), list(1))
twos <- list(list(.1, .9), list(.2, .8), list(.25, .75), list(.3,.7), list(.4, .6), list(.45, .55), list(.5, .5),
             list(.08, .92), list(.61, .39), list(.67, .33), list(.53, .47), list(.38, .62), list(.59, .41),
             list(.06, .94), list(.05, .95), list(.035, .965), list(.17, .83), list(.48, .52)
             , list(.18, .82), list(.14, .86), list(.11, .89), list(.23, .77), list(.26, .74), list(.21, .79)
             , list(.44, .56), list(.41, .59), list(.33, .67), list(.31, .69), list(.37, .63), list(.49, .51), list(.45, .55)
             , list(.04, .96), list(.07, .93), list(.9, .11), list(.2, .8), list(.2, .8))
threes <- list(list(.1,.2,.7),list(.1,.3,.6),list(.2,.3,.5),list(.25,.3,.45),list(.28,.34,.38), list(.4, .33, .27),
               list(.1, .14, .77), list(.44, .39, .17), list(.58, .23, .19), list(.07, .33, .6), list(.05, .08, .87)
               ,list(.12,.29,.59),list(.14,.28,.58),list(.19,.25,.64),list(.04,.44,.52),list(.07,.1,.83),list(.16,.35,.49)
               ,list(.32,.27,.41),list(.64,.2,.16),list(.19,.23,.68),list(.05,.41,.54),list(.07,.17,.76),list(.16,.39,.45))
fours <- list(list(.1, .15, .23, .52), list(.07, .15, .37, .41), list(.1, .27, .3, .33), list(.07, .13, .165, .635)
              , list(.12, .18, .26, .44), list(.18, .15, .23, .44), list(.22, .13, .26, .39), list(.06, .21, .29, .44)
              , list(.10, .19, .37, .34), list(.32, .08, .21, .39))
psiList <- c(ones, twos, threes, fours)
thetatest <- function(theta){
  trueNs <- sapply(1:length(real),function(i){length(real[[i]][[1]][[1]])})
  inferredNs <- sapply(1:length(res3),function(i){
    which.max(unlist(res3[[i]][[1]][3,])*dexp(1:6,theta))})
  list(length(which(trueNs==inferredNs)), data.frame('true'=trueNs,'inferred'=inferredNs))
}
L <- lapply(thetas, thetatest)
i.best <- which.max(sapply(1:length(L),function(i){L[[i]][[1]]}))
tab <- L[[i.best]][[2]]
best <- thetas[i.best]
acc <- length(which(tab[,1]==tab[,2]))/nrow(tab)
accs <- sapply(1:length(thetas),function(i){
  length(which(L[[i]][[2]][,1]==L[[i]][[2]][,2]))/nrow(L[[i]][[2]])})
plot(thetas, accs,type='l')

nset <- c(rep(1, length(ones)),rep(2, length(twos)),rep(3, length(threes)),
   rep(4, length(fours)))
f <- function(i){
thetafun <- function(theta){
  foo <- function(n){
    unlist(res3[[i]][[1]][3,])[n]*dexp(n,theta)
  }
  sapply(1:6, foo)
}
tab <- t(sapply(thetas, thetafun))
plot(thetas,tab[,1], type='l',col='black',ylim=c(min(tab),max(tab)),main=
       paste('sample = ',i, ', nclone = ', nset[i], sep=''))
lines(thetas,tab[,2], type='l',col='gray')
lines(thetas,tab[,3], type='l',col='red')
lines(thetas,tab[,4], type='l',col='blue')
lines(thetas,tab[,5], type='l',col='darkorange')
lines(thetas,tab[,6], type='l',col='purple')
sapply(1:6,function(x){length(which(sapply(nrow(tab), function(j){
  which.max(tab[j,])})==x))})
}

results <- data.frame(t(sapply(1:length(res3),f)))
colnames(results) <- c("1 clone","2 clones","3 clones","4 clones","5 clones",
                       "6 clones")
results$true.nclone <- nset