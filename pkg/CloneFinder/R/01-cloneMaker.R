########################################################################
### Part 1: Simplices
###
### The d-dimensional simplex is the set of nonnegative points in R^d
### that satisfy  x_1 + .. + x_d = 1. These points represent partitions
### of data into subsets, such as taking the fraction of cells that are
### part of each clone.

### Sample n points from the d-simplex
sampleSimplex <- function(n, d=5) {
    result <- matrix(NA, nrow=n, ncol=d)
    for (i in 1:n) {
        result[i,] <- diff(sort( c(0, 1, runif(d-1, 0, 1)) ))
    }
    result
}

### Generate a regularly spaced lattice of points in the d-simplex,
### with k points along each edge.
### Use symmetry to reduce the number we need to try.
generateSimplex <- function(k, d, reps=1){
  simplex <- t(xsimplex(d, k))
  for(i in 1:nrow(simplex)){
    simplex[i,] <- sort(simplex[i,], decreasing=TRUE)
  }
  simplex <- unique(simplex)
  psis <- t(sapply(rep(1:nrow(simplex), reps), function(i){simplex[i,]/k}))
  psis
}

########################################################################
### Part 2: WeightVectors
###
### We aloso refer to points in a simplex as "weight vectors" (since
### "vector of clonal fractions" is both too long and too specialized).
### We define these as an S4 class so we can build validitiy checking
### into the constructor, and so we can more reliably use them when we
### generate "Tumor" objects later.

setClass("WeightVector", slots=c(psi = "numeric"))
setValidity("WeightVector", function(object) {
  all( object@psi >= 0 ) & sum(object@psi) == 1
})
setMethod("initialize", "WeightVector", function(.Object, psi = 1, ...) {
  .Object <- callNextMethod(.Object) # in case this gets inherited
  if (any(is.na(psi))) {
    stop("Psi should not contain missing values.")
  }
  if (any(psi < 0)) {
    stop("Psi should not contain negative values.")
  }
  if (all(psi == 0)) {
    stop("At least one psi value must be larger than zero.")
  }
  .Object@psi = psi/sum(psi)
  .Object
})
WeightVector <- function(phi) {
  new("WeightVector", psi = phi)
}
setAs("WeightVector", "numeric", function(from) from@psi)


########################################################################
# Part 3: Simulating Data

############ CLONE ############
# A clone is a list of segments, where each segment represents exactly
# one pure compartment.

setClass("Clone", representation=list(
                    segments = "integer",
                    weights = "numeric"
                    ))

# We should be able to simulate a clone if we know the number of segments,
# the number of compartments, and the relative frequency/prevalence of
# each compartment.
Clone <- function(nSegments, weights=rep(1/5, 5), segnames=NULL) {
  # nSegments = integer, the number of segments
  # weights   = vector, the prevalence of each compartment

  if (nSegments < 1) stop("Number of segments must be positive.")
  
  # start with sanity checks on the weights
  if(any(is.na(weights))) stop("Missing weights.")
  if(any(weights < 0)) stop("Negative weights.")
  weights <- weights/sum(weights) # ensure they sum to 1

  # now sample the compartments to generate a clone
  segs <- sample(length(weights), nSegments, replace=TRUE, prob=weights)
  if (is.null(segnames)) segnames <- paste("Segment", 1:nSegments, sep='')
  if (length(segnames) != length(segs)) stop("Wrong number of segment names.")
  names(segs) <- segnames
  new("Clone", segments=segs, weights=weights)
}

############ TUMOR ############
# A tumor is a set of clones, each of which is associated with a
# fraction, subject to the constraint that the sum of the fractions
# equals one.

setClass("Tumor",
         slots = c(psi = "WeightVector",
                   clones = "list"))
setAs("Tumor", "list", function(from) {
  list(psi = as(from@psi, "numeric"),
       clones = from@clones)
})
setMethod("summary", signature("Tumor"), function(object, ...) {
    tdata <- object@data
    table(A=tdata[,1], B=tdata[,2])
})

getClone <- function(tumor, i) {
  if (inherits(tumor, "Tumor")) {
    tumor <- as(tumor, "list")
  }
  tumor$clones[[i]] # should we do error checking on 'i'?
}

### TODO: Wrap everything except psi into some sort of parameter class.
Tumor <- function(psi, rounds, nu=100, pcnv=0.5, norm.contam=FALSE, cnmax=4) {
  K <- length(which(psi > 0)) # number of clones
  ## Choose the number of copy number (CN) segments.
  total.segs <- round(runif(1, 250, 500))
  ## Distribute segments across the chromosomes.
  segsperchr <- as.vector(rdirichlet(1, chlens/1000000))
  segsperchr <- round(segsperchr*total.segs)
  segsperchr[segsperchr < 1] <- 1
  ## Start building a fake DNA copy data frame. Each segments
  ## needs a chromosome with start and end positions.
  chr <- unlist(lapply(1:24, function(i) {
    rep(i, segsperchr[i])
  }))
  ends <- lapply(1:24, function(i) {
    lens <- c(as.vector(rdirichlet(1, rep(1, segsperchr[i]))))
    sapply(1:length(lens), function(j) {sum(lens[1:j])})
  })
  ends <- lapply(1:length(ends), function(i) {
    round(chlens[chr[i]]*ends[[i]])
  })
  starts <- lapply(1:length(ends), function(i) {
    c(1, ends[[i]][1:(length(ends[[i]])-1)]+1)
  })
  ends <- unlist(ends)
  starts <- unlist(starts)
  ## Make the data frame.
  cnmat <- cbind(chr, unlist(starts), unlist(ends), rep(1, length(starts)), 
                 rep(1, length(starts)), 1:length(starts), rep(NA, length(starts)))
  colnames(cnmat) <- c('chr', 'start', 'end', 'A', 'B', 'seg', 'parent.index')
  ## Evolve some clones, starting with the currently normal one.
  startclone <- list('cn'=cnmat, 'seq'=NA)
  id.end <- 0
  clones <- list(startclone)
  while(length(clones) <= rounds) {
    parent.index <- sample(1:length(clones), 1)
    parent <- clones[[parent.index]]
    if(nu > 0 ) { # Add some mutations. Possibly CNV as well
      cnv <- sample(c(TRUE, FALSE), 1, prob = c(pcnv, 1 - pcnv))
    } else {      # Must be a CNV
      cnv <- TRUE
    }
    newclone.cn <- parent$cn
    newclone.cn[, which(colnames(newclone.cn)=='parent.index')] <- parent.index
    newclone.seq <- parent$seq

    ## handle copy number variants
    if(cnv) {
      allele <- sample(c('A', 'B'), 1)
      pool <- which(parent$cn[, which(colnames(cnmat)==allele)] > 0 &
                      parent$cn[, which(colnames(cnmat)==allele)]<cnmax)
      tochange <- sample(pool, 1)
      other.clones <- sapply(1:length(clones), function(j) {
        clones[[j]]$cn[tochange, which(colnames(cnmat)==allele)]
      })
      if(max(other.clones) > 1) {
        delta <- 1
      } else if(min(other.clones) < 1) {
        delta <- -1
      } else {
        delta <- sample(c(-1, 1), 1)
      }
      if(parent.index>1 & length(which(!is.na(unlist(newclone.seq)))) > 0) {
        muts.tochange <- intersect(which(newclone.seq[, i.seg]==tochange),
                                   which(newclone.seq[, i.allele]==allele))
        if(length(muts.tochange) > 0) {
          newclone.seq[muts.tochange, i.cn] <-
            as.numeric(newclone.seq[muts.tochange, i.cn]) + delta
        }
      }
      newclone.cn[tochange, which(colnames(cnmat)==allele)] <-
         as.numeric(newclone.cn[tochange, which(colnames(cnmat)==allele)]) + delta
    } #END: if(cnv)

    ## handle sequence mutations
    nmuts <- round(runif(1, nu/2, 1.5*nu))
    if(nmuts > 0){
      mutsperchr <- as.vector(rmultinom(1, nmuts, chlens/1000000))
      starts <- lapply(1:length(mutsperchr), function(i) {
        unique(round(sort(runif(mutsperchr[i], 1, chlens[i]))))
      })
      starts <- unlist(starts)
      mut.chrs <- unlist(lapply(1:length(mutsperchr), function(i) {
        rep(i, mutsperchr[i])
      }))
      mut.segs <- unlist(sapply(1:length(starts), function(i) {
        col.index <- which(colnames(cnmat)=='chr')
        tab <- cnmat[which(cnmat[, col.index]==mut.chrs[i]), ]
        bin1 <- which(tab[, which(colnames(tab)=='start')]<=starts[[i]])
        index <- bin1[which(tab[bin1, which(colnames(tab)=='end')]>=starts[[i]])]
        seg <- unname(tab[index, which(colnames(tab)=='seg')])
        seg
      }))
      alleles <- sample(c('A', 'B'), length(starts), replace=TRUE)
      mut.ids <- (id.end+1):(id.end+length(starts))
      seqmat <- cbind(mut.chrs, starts, mut.segs, mut.ids, rep(1, length(starts)), alleles)
      id.end <- max(mut.ids)
      colnames(seqmat) <- c('chr', 'start', 'seg', 'mut.id', 'mutated.copies', 'allele')
      newclone.seq <- rbind(newclone.seq, seqmat)
      rownames(newclone.seq) <- NULL
    } else {
      newclone.seq <- NA # No new mutations
    } #END: if(nmuts > 0)-else
    ## Make the new clone, remembewring what changed.
    if(length(clones) == 1 & nmuts > 0) {
      i.seg <- which(colnames(seqmat)=='seg')
      i.cn <- which(colnames(seqmat)=='mutated.copies')
      i.allele <- which(colnames(seqmat)=='allele')
    }
    newclone <- list('cn'=newclone.cn, 'seq'=newclone.seq)
    clones <- c(clones, list(newclone))
  } #END: while(length(clones) <= rounds)
  
  if(norm.contam==TRUE) { # first "clone" consists of the normal cells.
    sampled <- c(1, sample(2:length(clones), K-1, replace=FALSE))
  }else{
    sampled <- sample(2:length(clones), K, replace=FALSE)
  }
  
  clones.final <- lapply(1:length(sampled), function(i) {
    cndf <- as.data.frame(clones[[sampled[i]]]$cn)
    if(length(which(!is.na(unlist(clones[[sampled[i]]]$seq))))){
      seqdf <- as.data.frame(clones[[sampled[i]]]$seq)
      for(j in which(colnames(seqdf)!='allele')){
        seqdf[, j] <- as.numeric(as.character(seqdf[, j]))
        seqdf <- na.omit(seqdf[with(seqdf, order(seg, start)), ])
        rownames(seqdf) <- NULL
        seqdf
      }
      seqdf$normal.copies <- na.omit(unlist(lapply(1:nrow(cndf), function(j){
        total.cn <- cndf$A[j] + cndf$B[j]
        if(length(which(seqdf$seg==j))>0){
          normal.cns <- total.cn - seqdf$mutated.copies[which(seqdf$seg==j)]
        } else {
          normal.cns <- NA
        }
        normal.cns
      })))
    } else {
      seqdf <- NA
    }
    if(nu>0){
      output <- list('cn'=cndf, 'seq'=seqdf)
    }else{
      output <- list('cn'=cndf)
    }
    output
  })
  list('clones'=clones.final, 'psi'=psi)
  new("Tumor", clones = clones.final, psi = WeightVector(psi))
}


trueZ <- function(tumor) {
  Zmat <- array(0, dim=c(nrow(tumor@centers), # segments
                     nrow(tumor@pureCenters), # compartments
                     ncol(tumor@data)))       # clones
  for (wclone in 1:ncol(tumor@data)) {
    for(wseg in 1:nrow(tumor@centers)) {
      wcomp <- tumor@data[wseg, wclone]
      Zmat[wseg,wcomp,wclone] <- 1
    }
  }
  Zmat
}

############ SIMULATING DATASET ############

# input:
#   'centers; is the result of calling "findCenters"
#   'markers' is the vector with the number of markers per segment
generateData <- function(object) {
  if (!inherits(object, "Tumor"))
    stop(paste("Incorrect class of 'object':", class(object)))
  centers <- object@centers
  markers <- object@markers
  xy <- object@pureCenters
  # now we can genrate the data
  xvec <- rnorm(length(markers), centers$x, object@SEM[,1])
  yvec <- rnorm(length(markers), centers$y, object@SEM[,2])
  data.frame(x=xvec, y=yvec)
}

sizeplot <- function(simdata, tumor) {
    size <- 1 + round(sqrt(tumor@markers)/15)/2
    filt <- size > log(128)/3
    x <- tumor@pureCenters$x
    y <- tumor@pureCenters$y
    eps <- 0.1*diff(range(x))
    plot(x, y, type='n',
         xlim=c(min(x)-eps, max(x)+eps),
         ylim=c(min(y)-eps, max(y)+eps))
    text(x, y, 1:nrow(tumor@pureCenters))
    points(tumor@pureCenters, cex=15, col='gray25')
    points(simdata$x[filt], simdata$y[filt],
           cex=size[filt], col="gray35", pch=16)
    invisible(simdata)
}


####################################################
# from GlobalMaxLike
#  KRC: DO we need these any more? Or are they part of an
# old approach that didn't generalize?

makeLatent <- function(vector){
  clone <- matrix(0, ncol=5, nrow=length(vector))
  for (i in 1:length(vector)){
    clone[i, vector[i]] <- 1
  }
  clone
}

estBetaParams <- function(vec) {
  mu <- mean(vec)
  s2 <- var(vec)
  temp <- mu*(1-mu)/s2 - 1
  alpha <- mu*temp
  beta <- (1 - mu)*temp
  c(alpha = alpha, beta = beta)
}

# this depended on lots of global variables.
# I hope we don't really use it or need it.
maxfun <- function(segment, newphiset, paramtable, temp){
  d <- temp[segment,]
  range <- 1:1000
  range <- range[!range == segment]
  optimfun <- function(n){
    distfun <- function(segment){
      log(dbeta(d[n], paramtable[segment, 1], paramtable[segment, 2]))
    }
    sum(sapply(range, distfun))+log(d[n]/sum(d))
  }
  d <- sapply(1:10000, optimfun)
  i <- which.max(d)
  newphiset[i,]
}
