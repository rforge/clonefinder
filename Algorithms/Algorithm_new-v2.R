library(gtools)
library(combinat)
names <- as.vector(read.table("D:\\Mark\\Simulations2\\namelist.txt")[,1])

callfun <- function(i){
  string <- names[i]
  #dat <- get(load(paste("F:\\Simulations\\", string, ".rda", sep="")))
  dat <- get(load(paste("D:\\Mark\\Simulations2\\", string, ".rda", sep="")))
  phis <- dat[[3]]
  phiv <- as.vector(phis)
  zeds <- dat[[4]]
  dens <- density(phiv, bw="SJ", n=1024)$y
  truepsis <- dat[[1]]
  combos <- unlist(sapply(1:length(truepsis), function(n){combn(truepsis, n, sum)}))
  means <- c(0, combos)
  den <- density(phiv, bw="SJ")$y
  things <- list(truepsis, means, dens, phiv, phis, zeds)
  things
}

f <- function(sim){
set <- callfun(sim)
truepsis <- set[[1]]
nclones <- length(truepsis)
combos <- unlist(sapply(1:nclones, function(n){combn(truepsis, n, sum)}))
vec <- set[[3]]
phiv <- set[[4]]
phiset <- set[[5]]
zeds <- set[[6]]
under <- phiv[phiv<.5]
over <- 1-phiv[phiv>.5]
newphiv <- sort(c(under, over), decreasing=FALSE)
fold <- sort(c(vec[1:(length(vec)/2)], rev(vec[(length(vec)/2+1):length(vec)])))
final <- FALSE

foo <- function(start){ 
postfun <- function(v){
  v - sort(v)
  p1 <- v[1]
  p2 <- v[2]
  p3 <- v[3]
  a <- p1
  b <- p2
  c <- p3
  d <- p1+p2
  e <- p1+p3
  f <- p2+p3  
  letters <- c(a, b, c, d, e, f)
  pnames <- c("p1", "p2", "p3", "p1+p2", "p1+p3", "p2+p3")
  weights <-c(sum(phiset[,1]),sum(phiset[,2]),sum(phiset[,3]),sum(phiset[,4]),
              sum(phiset[,5]))/
    sum(phiv)
  combins <- unlist(sapply(1:nclones, function(n){combn(v, n, sum)}))
  upper <- 1-(1-max(combins))/2
  phiset.trim <- phiset[which(sapply(1:1000,function(i){max(phiset[i,]) < upper})==TRUE),]
  ones <- phiset[which(sapply(1:1000, function(i){max(phiset[i,]) >= upper})==TRUE),]
  i.ones <- which(sapply(1:1000, function(i){max(phiset[i,]) >= upper})==TRUE)
  columns <- sapply(1:nrow(ones), function(i){which.max(ones[i,])})
  fun <- function(j){
    v <- rep(0, 5)
    v[j] <- 1
    v
  }
  oneset <- t(sapply(columns, fun))
  phiset <- phiset.trim
  
  findBestMatch <- function(row){
    vec <- phiset[row,]
    c1 <- c(a, b, c, 2, 3)
    c2 <- c(a, f, 2, 3, 4)
    c3 <- c(b, e, 2, 3, 4)
    c4 <- c(c, d, 2, 3, 4)
    c5 <- c(1, 2, 3, 4, 5)
    L <- list(c1, c2, c3, c4, c5)
    permfun <- function(i){
      perms <- permutations(5, 5, L[[i]], set=FALSE)
      perms[perms==5 | perms==4 | perms==3 | perms==2]<- 0
      unique(perms)
    }
    allperms <- sapply(1:length(L), permfun)
    allperms <- rbind(allperms[[1]],allperms[[2]],allperms[[3]],allperms[[4]],
                      allperms[[5]])
    probfun <- function(permrow){
      sum(abs(allperms[permrow,]-vec))
    }
    probs <- sapply(1:nrow(allperms), probfun)
    prob.min <- probfun(which.min(probs))
    row.min <- allperms[which.min(probs),]
    matchfun <- function(q){
      x <- L[[q]]
      x[x==5 | x==4 | x==3 | x==2]<- 0
      bin <- sort(x)==sort(row.min)
      length(which(bin %in% TRUE))==5
    }
    elements <- L[[which(sapply(1:length(L), matchfun)==TRUE)]]
    elements <- sort(elements[elements < 1])
    if(length(elements)==0){
      combins <- c("p1+p2+p3")
      elements <- L[[which(sapply(1:length(L), matchfun)==TRUE)]]
    }else{
    combins <- pnames[which(letters %in% elements)]
    }
    assign <- rep(0, 5)
    whichfun <- function(ind){
      which(row.min %in% elements[ind])
    }
    order <- unlist(sapply(1:length(elements), whichfun))
    assign[order] <- combins
    #prod1 <- prod(weights[which(row.min %in% c(a, b, c))])
    #prod2 <- prod(weights[which(row.min %in% c(d, e, f))]^2)
    #prod3 <- prod(weights[which(row.min %in% 1)]^3)
    #prod1*prod2*prod3
    if(final==TRUE){
     output <- list(prob.min, assign)
    }else{
      output <- prob.min
    }
    output
  }
  if(final==TRUE){
  picks <- sapply(1:nrow(phiset), function(z){findBestMatch(z)[[1]]})
  assignments <- t(sapply(1:nrow(phiset), function(z){findBestMatch(z)[[2]]}))
  output <- list(prod(picks^.01), assignments)
  }else{
    picks <- sapply(1:nrow(phiset), function(z){findBestMatch(z)})
    output <- prod(picks^.01)
  }
  output
}
if(final==TRUE){
posteriors <- unlist(sapply(1:nrow(start), function(i){postfun(start[i,])[[1]]}))
assignments <- sapply(1:nrow(start), function(i){postfun(start[i,])[2]})
assignments <- assignments[which.min(posteriors)]
pick <- sort(start[which.min(posteriors),])
output <- list(pick, assignments)
}else{
  posteriors <- unlist(sapply(1:nrow(start), function(i){postfun(start[i,])}))
  output <- sort(start[which.min(posteriors),])
}
output
}

cset <- t(combn(1:8, 3))
start1 <- t(sapply(1:nrow(cset), function(i){sort(cset[i,]/sum(cset[i,]))}))
resids <- sapply(1:nrow(start1), function(n){sum(abs(sort(start1[n,])-truepsis))})
truebest <- start1[which.min(resids),]

#note: currently assignments must be reconciled with use of multiple iterations.
iter1 <- foo(start1)
newset <- t(combn(1:10, 3))
start2 <- t(sapply(1:nrow(newset), function(i){sort(newset[i,]/sum(newset[i,]))}))
rangefun <- function(i){
  innerfun <- function(col){
  start2[i,][col] > (iter1[col])-.06 & start2[i,][col] < (iter1[col])+.06
  }
  bin <- sapply(1:length(iter1), innerfun)
  length(which(bin==TRUE))==length(iter1)
}
bin <- sapply(1:nrow(newset), rangefun)
start.trim <- start2[which(bin==TRUE),]

iter2 <- foo(start.trim)
iter2 <- iter2[[1]]
newset2 <- t(combn(1:18, 3))
start3 <- t(sapply(1:nrow(newset2), function(i){sort(newset2[i,]/sum(newset2[i,]))}))
rangefun <- function(i){
  innerfun <- function(col){
    start3[i,][col] > (iter2[col])-.01 & start3[i,][col] < (iter2[col])+.01
  }
  bin <- sapply(1:length(iter2), innerfun)
  length(which(bin==TRUE))==length(iter2)
}
bin <- sapply(1:nrow(newset2), rangefun)
start.trim2 <- start3[which(bin==TRUE),]
final <- TRUE
iter3 <- foo(start.trim2)
preZ <- iter3[2][[1]][[1]]
psis <- iter3[1][[1]]

clonelist <- c("p1", "p2", "p3")
zGen <- function(clone){
  clone <- clonelist[clone]
  outerfun <- function(i){
    innerfun <- function(j){
      grepl(clone, toString(preZ[i,j]))
    }
    zeroes <- rep(0, 5)
    zeroes[which((sapply(1:5, innerfun))==TRUE)] <- 1
    zeroes
  }
  z <- t(sapply(1:nrow(preZ), outerfun))
  zed <- matrix(rep(0, 5000), nrow=1000, ncol=5)
  zed[i.ones,] <- oneset
  sample <- 1:1000
  zed[sample[is.na(pmatch(sample,i.ones))],] <- z
  zed
}
zedList <- lapply(1:length(clonelist), zGen)
accuracy <- c(1-sum(zedList[[1]]!=zeds[[1]])/(2*1000), 1-sum(zedList[[2]]!=zeds[[2]])/
                (2*1000), 1-sum(zedList[[3]]!=zeds[[3]])/(2*1000))
obj <- list(zedList, accuracy, psis, truepsis)
#save(obj, file=paste("D:\\Mark\\Simulations2\\Results\\sim2results.",sim, ".rda",
#sep=""))
obj
}

#problem line: bin <- sapply(1:nrow(newset2), rangefun)

l <- c(105, 106, 108, 115, 117, 119, 122, 123, 156, 157, 162, 173, 174)
sapply(l, f)
#start time: 4:37
#end time: 