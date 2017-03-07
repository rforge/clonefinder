load('C:\\Users\\Mark\\OneDrive - The Ohio State University\\psis3.rda')
load('C:\\Users\\Mark\\Dropbox\\Lab\\cll\\Simulations\\grid1.rda')
grid1 <- grid
lens.psi <- sapply(1:nrow(psis3), function(i){length(which(psis3[i,]>0))})
vecpool <- grid1
vecpool <- lapply(1:length(vecpool), function(i){
  l <- length(vecpool[[i]])
  c(vecpool[[i]], rep(0, 5 - l))
})
vecpool <- unname(Reduce(rbind, vecpool))
psipool <- psis3

###Simulation:
load('C:\\Users\\Mark\\Dropbox\\Lab\\cll\\Simulations\\sims1.rda')
x <- 3
datamat <- sims1[[x]][[2]]
datavec <- c(datamat$X, datamat$Y)
markers <- rep(datamat$markers, 2)
truemat <- sims1[[x]][[3]]
truepsi <- as.vector(sims1[[x]][[4]])