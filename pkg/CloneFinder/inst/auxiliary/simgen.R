
###Generate a sample simulation population:
#set: 1 (pcnv=0, nu=30, sigma.lrr=.15, sigma.baf=.03, sigma.reads=25); k = 1:5
#set: 2 (pcnv=1, nu=30, sigma.lrr=.15, sigma.baf=.03, sigma.reads=25); k = 1:5
#set: 3 (pcnv=1, nu=0, sigma.lrr=.15, sigma.baf=.03, sigma.reads=25); k = 1:5
gen.sims <- FALSE

#path1 <- 'C:\\Users\\Mark\\OneDrive - The Ohio State University\\Integrated\\parOpt\\sim'
#path2 <- 'C:\\Users\\Mark\\OneDrive - The Ohio State University\\Integrated\\parOpt\\data'
path1 <- 'C:\\Users\\Mark\\OneDrive - The Ohio State University\\clonetools\\sim2'
path2 <- 'C:\\Users\\Mark\\OneDrive - The Ohio State University\\clonetools\\dat2'
if(gen.sims){
  pars.default <- list('rounds'=400, 'nu'=0, 'pcnv'=1, 'norm.contam'=FALSE)
  dataPars <- list('snps.seq'=1000000, 'snps.cgh'=600000, 'mu'=70, 'sigma.reads'=25, 'sigma0.lrr'=.15, 
                   'sigma0.baf'=.03, 'density.sigma'=.1)
  threshold <- .04
  dirs <- list.dirs(path1)
  dirs <- dirs[-1]
  if(length(dirs)>0){
    setIDs <- as.numeric(sapply(1:length(dirs), function(j){strsplit(dirs[j], split='sims-')[[1]][2]}))
    setID <- max(setIDs) + 1
  }else{
    setID <- 1
  }
  dir.create(paste(path1, '\\sims-', setID, sep=''))
  dir.create(paste(path2, '\\data-', setID, sep=''))
  simpath <- paste(path1, '\\sims-', setID, sep='')
  datapath <- paste(path2, '\\data-', setID, sep='')
  #psis <- genSimplex(20, 5, reps=1)
  psis <- get(load('C:\\Users\\Mark\\OneDrive - The Ohio State University\\Integrated\\psis.100.rda'))
  psis <- psis[which(sapply(1:nrow(psis), function(j){length(which(psis[j, ]==0 | psis[j, ]>threshold))==ncol(psis)})), ]
  ks <- sapply(1:nrow(psis), function(j){length(which(psis[j, ]>0))})
  seed <- 123
  set.seed(seed)
  sampled <- unlist(lapply(1:5, function(x){sample(which(ks==x), 60, replace=TRUE)}))
  psis <- psis[sampled, ]
  n <- nrow(psis)
  for(j in 1:nrow(psis)){
    go <- FALSE
    while(go==FALSE){
      tumor <- try(tumorGen(psis[j, ], pars.default$rounds, pars.default$nu, pars.default$pcnv, pars.default$norm.contam), silent=TRUE)
      if(class(tumor)!='try-error'){
        data <- dataGen(tumor, dataPars$snps.seq, dataPars$snps.cgh, dataPars$mu, dataPars$sigma.reads, 
                        dataPars$sigma0.lrr, dataPars$sigma0.baf, dataPars$density.sigma)
      }
      if(class(tumor)!='try-error' & class(data)!='try-error'){
        go <- TRUE
      }
    }
    filename.tumor <- paste(simpath, '\\', 'sim-', setID, '-', j, '.rda', sep='')
    filename.data <- paste(datapath, '\\', 'dat-', setID, '-', j, '.rda', sep='')
    sim <- list('tumor'=tumor)
    save(sim, file=filename.tumor)
    dat <- list('dat'=data)
    save(dat, file=filename.data)
  }
  if(pars.default$nu > 0 & pars.default$pcnv > 0){
    alt.types <- 'both'
  }else if(pars.default$nu > 0 & pars.default$pcnv == 0){
    alt.types <- 'mut'
  }else if(pars.default$nu == 0 & pars.default$pcnv > 0){
    alt.types <- 'cnv'
  }else{
    alt.types <- 'none'
  }
  metadata <- list('tumor.params'=pars.default, 'data.params'=dataPars, 'alt.types'=alt.types, 'date'=Sys.Date(), 
                  'setID'=setID, 'seed'=seed)
  save(metadata, file=paste(simpath, '\\metadata-', setID, '.rda', sep=''))
}


###Mixing real data:
mix <- function(datList, psi, datapath, simpath, index, version, alter=FALSE){
  breaks <- lapply(1:22, function(j){
    brks <- sort(unique(unlist(sapply(1:length(datList), function(k){c(datList[[k]][datList[[k]]$chrom==j, ]$loc.start, 
      datList[[k]][datList[[k]]$chrom==j, ]$loc.end)}))))
    starts <- brks[1:(length(brks)-1)]
    ends <- c(starts[2:length(starts)]-1, brks[length(brks)])
    data.frame('start'=starts, 'end'=ends, 'chr'=rep(j, length(starts)))
  })
  coords <- Reduce(rbind, breaks)
  columns <- sample(1:2, length(psi), replace=TRUE)
  shifts <- sapply(1:length(datList), function(j){
    lrrs <- unlist(sapply(1:nrow(datList[[j]]), function(k){rep(datList[[j]]$seg.median[k], datList[[j]]$num.mark[k])}))
    dens <- density(na.omit(2*10^lrrs))
    dens$x[which.max(dens$y)] - 2
  })
  d2 <- lapply(1:length(datList), function(j){
    temp <- datList[[j]]
    res <- t(sapply(1:nrow(coords), function(k){
      chrdat <- temp[as.character(temp$chrom)==as.character(coords$chr[k]), ]
      seg <- chrdat[chrdat$loc.end>=coords$start[k] & chrdat$loc.start<=coords$end[k], ]
      #Afreq <- seg$AA*2 + seg$AB
      #Bfreq <- seg$BB*2 + seg$AB
      #baf <- Bfreq/(Afreq + Bfreq)
      baf <- seg$AvgBAF
      tcn <- 2*10^(seg$seg.median)
      tcn <- tcn - shifts[j]
      if(length(which(tcn<0))>0){
        tcn[which(tcn<0)] <- -tcn[which(tcn<0)]
      }
      X <- tcn*(1-baf)
      Y <- tcn*baf
      if(nrow(seg)>0){
        output <- c(X, Y)
      }else{
        output <- rep(NA, 2)
      }
      output
    }))
    colnames(res) <- c('X', 'Y')
    res
  })
  data.mixed <- as.data.frame(Reduce('+', lapply(1:length(d2), function(k){psi[k]*d2[[k]]})))
  notna <- which(sapply(1:nrow(data.mixed), function(x){length(which(is.na(data.mixed[x, ])))==0}))
  for(j in 1:length(d2)){
    d2[[j]] <- d2[[j]][notna, ]
  }
  data.mixed <- data.mixed[notna, ]
  starts <- coords$start[notna]
  ends <- coords$end[notna]
  chrs <- coords$chr[notna]
  data.mixed$chr <- chrs
  data.mixed$seg <- 1:length(notna)
  data.mixed$LRR <- log10((data.mixed$X + data.mixed$Y)/2)
  data.mixed$BAF <- data.mixed$X/(data.mixed$X + data.mixed$Y)
  data.mixed$markers <- unlist(sapply(1:22, function(i){
    pos.chr <- pos[pos$Chr==i, ]
    data.chr <- data.mixed[data.mixed$chr==i, ]
    sapply(1:nrow(data.chr), function(j){
      length(which(pos.chr >= starts[which(data.mixed$chr==i)[j]] & pos.chr <= ends[which(data.mixed$chr==i)[j]]))
    })
  }))
  data.mixed$start <- starts
  data.mixed$end <- ends
  data.mixed <- data.mixed[, c(3, 4, 5, 6, 1, 2, 7, 8, 9)]
  clones <- lapply(1:length(d2), function(k){
    #tcn <- d2[[k]][, 1] + d2[[k]][, 2]
    #BAF <- d2[[k]][, 1]/tcn
    #BAF[which(tcn==0)] <- 0
    #df <- as.data.frame(round(t(sapply(1:length(tcn), function(l){tcn[l]*c(BAF[l], 1-BAF[l])}))))
    df <- data.frame('A'=rep(1, nrow(d2[[1]])), 'B'=rep(1, nrow(d2[[1]])))
    colnames(df) <- c('A', 'B')
    df$chr <- data.mixed$chr
    df$start <- starts
    df$end <- ends
    df$seg <- 1:length(starts)
    df$parent.index <- rep(0, length(starts))
    df$markers <- data.mixed$markers
    df <- df[, c(3, 4, 5, 1, 2, 6, 7, 8)]
    list('cn'=df, 'seq'=NULL)
  })
  for(k in 2:nrow(clones[[1]]$cn)){
    bin1 <- data.mixed$X[k]==data.mixed$X[k-1] & data.mixed$Y[k]==data.mixed$Y[k-1]
    bin2 <- clones[[1]]$cn$chr[k]==clones[[1]]$cn$chr[k-1]
    if(is.na(bin1)){
      bin1 <- FALSE
    }
    if(bin1==TRUE & bin2==TRUE){
      for(l in 1:length(clones)){
        clones[[l]]$cn$start[k] <- clones[[l]]$cn$start[k-1]
        clones[[l]]$cn$markers[k] <- clones[[l]]$cn$markers[k] + clones[[l]]$cn$markers[k-1]
        clones[[l]]$cn[k-1, ] <- rep(NA, ncol(clones[[1]]$cn))
      }
      data.mixed$start[k] <- data.mixed$start[k-1]
      data.mixed$markers[k] <- data.mixed$markers[k-1]
      data.mixed[k-1, ] <- rep(NA, ncol(data.mixed))
    }
  }
  notna <- unique(unlist(lapply(1:length(clones), function(j){which(!is.na(clones[[j]]$cn$seg))})))
  for(k in 1:length(clones)){
    clones[[k]]$cn <- clones[[k]]$cn[notna, ]
    clones[[k]]$cn$seg <- rownames(clones[[k]]$cn) <- 1:nrow(clones[[k]]$cn)
    clones[[k]]$cn
  }
  data.mixed <- data.mixed[notna, ]
  dat <- list('cn.data'=data.mixed, 'seq.data'=NULL)
  dat$cn.data <- na.omit(dat$cn.data)
  dat$cn.data$seg <- rownames(dat$cn.data) <- 1:nrow(dat$cn.data)
  if(alter){
    alt.pools <- lapply(1:length(psi), function(j){
      x <- dat$cn.data$X
      y <- dat$cn.data$Y
      markers <- dat$cn.data$markers
      which(markers>1000 & x-shifts[j]/2>.95 & x-shifts[j]/2<1.05 & y-shifts[j]/2>.95 & y-shifts[j]/2<1.05)
    })
    altN <- sample(1:3, length(psi), prob=c(.6, .3, .1), replace=TRUE)
    altered <- lapply(1:length(psi), function(j){sample(alt.pools[[j]], altN[j])})
    change <- lapply(1:length(psi), function(j){sample(c(-1, 1), altN[j], replace=TRUE)})
    for(l in 1:length(psi)){
      if(length(altered[[l]])>0){
        for(m in 1:length(altered[[l]])){
          allele <- sample(c(1, 2), 1)
          clones[[l]]$cn[altered[[l]][[m]], allele+3] <- 
            clones[[l]]$cn[altered[[l]][[m]], allele+3] + change[[l]][[m]]
          dat$cn.data[altered[[l]][[m]], allele+4] <- dat$cn.data[altered[[l]][[m]], allele+4] + 
            change[[l]][[m]]*psi[l]
          if(dat$cn.data[altered[[l]][[m]], allele+4] < 0){
            dat$cn.data[altered[[l]][[m]], allele+4] <- -dat$cn.data[altered[[l]][[m]], allele+4]
          }
        }
      }
    }
  }
  tumor <- list('clones'=clones, 'psi'=psi, 'altered'=altered, 'change'=change)
  save(dat, file=paste(datapath, '/mixdat', version, '-', index, '.rda', sep=''))
  save(tumor, file=paste(simpath, '/mixsim', version, '-', index, '.rda', sep=''))
}

gen.mixtures <- FALSE
if(gen.mixtures){
  version <- 1
  threshold <- .04
  if(version==1){
    version <- ''
  }
  seed <- 12345
  set.seed(seed)
  path <- 'C:/Users/Mark/OneDrive - The Ohio State University/clonetools'
  rewind.hapmap <- get(load(paste(path, '\\SNP_array_data\\rewind.hapmap.rda', sep='')))
  rewind.cll <- get(load(paste(path, '\\', 'SNP_array_data\\rewind_both.rda', sep='')))
  rewind.cll <- rewind.cll[, c(1:6, 11)]
  rewind <- rbind(rewind.cll, rewind.hapmap)
  datapath <- paste('C:/Users/Mark/OneDrive - The Ohio State University/clonetools/dat/data-mix2', version, sep='')
  simpath <- paste('C:/Users/Mark/OneDrive - The Ohio State University/clonetools/sim/sims-mix2', version, sep='')
  psipool <- get(load('C:/Users/Mark/OneDrive - The Ohio State University/clonetools/clonefinder/psis.100.rda'))
  psipool <- psipool[which(sapply(1:nrow(psipool), function(j){length(which(psipool[j, ]==0 | psipool[j, ]>threshold))==ncol(psipool)})), ]
  ks <- sapply(1:nrow(psipool), function(i){length(which(psipool[i, ]>0))})
  N <- 300 #Let's to 60 apiece for k=1, 2, 3, 4, and 5
  indices <- c(sample(which(ks==1), 60, replace=TRUE), sample(which(ks==2), 60, replace=TRUE), sample(which(ks==3), 60, replace=FALSE), 
               sample(which(ks==4), 60, replace=FALSE), sample(which(ks==5), 60, replace=FALSE))
  psis <- lapply(1:N, function(j){psipool[indices[j], ][which(psipool[indices[j], ]>0)]})
  #IDpool <- c('CL001', 'CL060', 'CL082', 'CL084', 'CL102', 'CL185', 'NA19193')
  IDpool <- c('NA07019', 'NA12234', 'NA12249', 'NA12753', 'NA12761', 'NA18545', 'NA18975', 'NA18999', 'NA18517')
  IDs <- lapply(1:N, function(j){sample(IDpool, length(which(psis[[j]]>0)), replace=FALSE)})
  metadata <- list('IDs'=IDs, 'psis'=psis, 'seed'=seed)
  save(metadata, file=paste(simpath, '\\metadata.rda', sep=''))
  pos <- get(load('C:/Users/Mark/OneDrive - The Ohio State University/clonetools/pos.rda'))
  sapply(1:length(psis), function(i){mix(lapply(1:length(IDs[[i]]), function(j){
    rewind[rewind$SamID==IDs[[i]][j], ]}), psis[[i]], datapath, simpath, i, version, alter=TRUE)})
}

check <- FALSE
if(check){
  load('C:\\Users\\Mark\\OneDrive - The Ohio State University\\clonetools\\sim\\sims-mix\\metadata-mix.rda')
  index <- 196
  sams <- metadata$IDs[[index]]
  first.two1 <- paste(strsplit(sams[1], split='')[[1]][1:2], collapse='')
  first.two2 <- paste(strsplit(sams[2], split='')[[1]][1:2], collapse='')
  folder1 <- folder2 <- 'data-CLL'
  if(first.two1=='NA'){
    folder1 <- 'data-hapmap'
  }
  if(first.two2=='NA'){
    folder2 <- 'data-hapmap'
  }
  dat1 <- get(load(paste('C:\\Users\\Mark\\OneDrive - The Ohio State University\\clonetools\\dat\\', 
                         folder1, '\\dat-', sams[1], '.rda', sep='')))
  dat2 <- get(load(paste('C:\\Users\\Mark\\OneDrive - The Ohio State University\\clonetools\\dat\\', 
                         folder2, '\\dat-', sams[2], '.rda', sep='')))
  mix.dat <- get(load(paste('C:\\Users\\Mark\\OneDrive - The Ohio State University\\clonetools\\dat\\data-mix\\mixdat-', index, '.rda', sep='')))
  dat1 <- dat1$dat$cn.data
  dat2 <- dat2$dat$cn.data
  mix.dat <- mix.dat$cn.data
  psi <- metadata$samplePsis[index, ]
  cnplotfun(dat1, 2, 'dat1')
  cnplotfun(dat2, 2, 'dat1')
  cnplotfun(mix.dat, 2, 'mixture')
  mixplotfun(mix.dat, dat1, dat2, 2)
}

###
plotSetTraits <- function(setDir, seq=TRUE, cn=TRUE, merge=FALSE, kmax=5){
  simfiles <- list.files(setDir)
  metadata <- simfiles[1]
  simfiles <- simfiles[2:length(simfiles)]
  sims <- lapply(1:length(simfiles), function(j){get(load(paste(setDir, '\\', simfiles[j], sep='')))})
  if(length(sims[[1]])==1){
    sims <- lapply(1:length(sims), function(x){sims[[x]]$tumor})
  }
  ks <- sapply(1:length(sims), function(j){length(which(sims[[j]]$psi>0))})
  psis <- t(sapply(1:length(sims), function(j){
    psi <- sims[[j]]$psi
    if(length(psi)<kmax){
      psi <- c(psi, rep(0, kmax-length(psi)))
    }
    psi
  }))
  frac.altered <- sapply(1:length(sims), function(j){
    sim <- sims[[j]]
    lengths <- sim$clones[[1]]$cn$end + 1 - sim$clones[[1]]$cn$start
    sum(lengths[unique(unlist(lapply(1:length(sim$clones), function(x){
      which(sim$clones[[x]]$cn$A!=1 | sim$clones[[x]]$cn$B!=1)})))])/sum(lengths)
  })
  if(merge){
    sim <- sims[[j]]
    lengths <- sim$clones[[1]]$cn$end + 1 - sim$clones[[1]]$cn$start
    altered <- sapply(1:length(sim$clones), function(x){sim$clones[[x]]$cn$A!=1 | sim$clones[[x]]$cn$B!=1})
    altered.indices <- which(altered)
    diffs <- diff(altered.indices)
    #To be continued?
  }else{
    length.altered <- lapply(1:length(sims), function(j){
      sim <- sims[[j]]
      lengths <- sim$clones[[1]]$cn$end + 1 - sim$clones[[1]]$cn$start
      lengths[unique(unlist(lapply(1:length(sim$clones), function(x){
        which(sim$clones[[x]]$cn$A!=1 | sim$clones[[x]]$cn$B!=1)})))]
    })
    nCNVs <- sapply(1:length(length.altered), function(z){length(length.altered[[z]])}) 
  }
  length.altered <- unlist(length.altered)
  mutations <- sapply(1:length(sims), function(j){
    sim <- sims[[j]]
    if(length(sim)==1){
      sim <- sim$tumor
    }
    length(unique(unlist(lapply(1:length(sim$clones), function(x){sim$clones[[x]]$seq$mut.id}))))
  })
  if(cn==TRUE & seq==FALSE){
    #par(mfrow=c(1, 2))
    #hist(ks, breaks=c(.5, 1.5, 2.5, 3.5, 4.5, 5.5), xlab='Number of Clones', main='')
    plot(density(unlist(psis), bw=.05), xlab='Psi', xlim=c(0, 1), main='F) Density of Psi Values across Simulations')
    plot(density(frac.altered, bw=.01), xlab='Genome Altered Fraction', xlim=c(0, 1), main='G) Genome AF Density among Simulations')
    #plot(density(length.altered, bw=.01), xlab='CNV Lengths', xlim=c(0, max(length.altered)), main='')
    #hist(nCNVs, xlab='CNVs per Sample', main='')
  }else if(cn==FALSE & seq==TRUE){
    #par(mfrow=c(1, 2))
    #hist(ks, breaks=c(.5, 1.5, 2.5, 3.5, 4.5, 5.5), xlab='Number of Clones', main='')
    plot(density(unlist(psis), bw=.05), xlab='Psi', xlim=c(0, 1), main='A) Density of Psi Values across Simulations')
    plot(density(mutations, bw=5), xlab='Mutations per Sample', xlim=c(0, max(mutations)), main='B) Density of Mutations per Simulation')
  }else if(cn & seq){
    #par(mfrow=c(1, 3))
    #hist(ks, breaks=c(.5, 1.5, 2.5, 3.5, 4.5, 5.5), xlab='Number of Clones', main='')
    plot(density(unlist(psis), bw=.05), xlab='Psi', xlim=c(0, 1), main='C) Density of Psi Values across Simulations')
    plot(density(frac.altered, bw=.01), xlab='Genome Altered Fraction', xlim=c(0, 1), main='D) Genome AF Density among Simulations')
    #plot(density(length.altered, bw=.01), xlab='CNV Lengths', xlim=c(0, max(length.altered)), main='')
    #hist(nCNVs, xlab='CNVs per Sample', main='')
    plot(density(mutations, bw=5), xlab='Mutations per Sample', xlim=c(0, max(mutations)), main='E) Density of Mutations per Simulation')
  }
}

png('C:/Users/Mark/OneDrive - The Ohio State University/clonetools/cf_images/fig_2.png', height=500, width=900, res=100)
layout(matrix(c(1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 6, 7, 7, 7), 3, 6, byrow = TRUE))
plotSetTraits('C:\\Users\\Mark\\OneDrive - The Ohio State University\\clonetools\\sim\\sims-1', cn=FALSE)
plotSetTraits('C:\\Users\\Mark\\OneDrive - The Ohio State University\\clonetools\\sim\\sims-2')
plotSetTraits('C:\\Users\\Mark\\OneDrive - The Ohio State University\\clonetools\\sim\\sims-3', seq=FALSE)
dev.off()

#png('C:\\Users\\Mark\\OneDrive - The Ohio State University\\clonetools\\images\\sims-1.png', width=10, height=3.5, units="in", res=300)
#plotSetTraits('C:\\Users\\Mark\\OneDrive - The Ohio State University\\clonetools\\sim\\sims-1', cn=FALSE)
#dev.off()
#png('C:\\Users\\Mark\\OneDrive - The Ohio State University\\clonetools\\images\\sims-2.png', width=15, height=3.5, units="in", res=300)
#plotSetTraits('C:\\Users\\Mark\\OneDrive - The Ohio State University\\clonetools\\sim\\sims-2')
#dev.off()
#png('C:\\Users\\Mark\\OneDrive - The Ohio State University\\clonetools\\images\\sims-3.png', width=10, height=3.5, units="in", res=300)
#plotSetTraits('C:\\Users\\Mark\\OneDrive - The Ohio State University\\clonetools\\sim\\sims-3', seq=FALSE)
#dev.off()
#png('C:\\Users\\Mark\\OneDrive - The Ohio State University\\clonetools\\images\\sims-mix.png', width=10, height=3.5, units="in", res=300)
#plotSetTraits('C:\\Users\\Mark\\OneDrive - The Ohio State University\\clonetools\\sim\\sims-mix', seq=FALSE, merge=FALSE)
#dev.off()

###
imgpath <- 'C:\\Users\\Mark\\OneDrive - The Ohio State University\\clonetools\\images'
metadata <- get(load('C:\\Users\\Mark\\OneDrive - The Ohio State University\\clonetools\\sim\\sims-3\\metadata-3.rda'))

plotSNP <- function(chr.sim, chr.dat, j, sam, homFrac=.68, mu.0=0, sigma.0=.03){
  dat <- get(load(paste('C:\\Users\\Mark\\OneDrive - The Ohio State University\\clonetools\\dat\\data-3', '\\dat-3-', j, '.rda', sep='')))
  sim <- get(load(paste('C:\\Users\\Mark\\OneDrive - The Ohio State University\\clonetools\\sim\\sims-3', '\\sim-3-', j, '.rda', sep='')))
  segdata <- dat$dat$cn.data
  #sigma0.lrr <- metadata$data.params$sigma0.lrr
  sigma0.lrr <- .12
  sigma0.baf <- metadata$data.params$sigma0.baf
  segs <- segdata[segdata$chr==chr.sim, ]
  x.mu <- unlist(sapply(1:nrow(segs), function(k){rep(segs$X[k], segs$markers[k])}))
  y.mu <- unlist(sapply(1:nrow(segs), function(k){rep(segs$Y[k], segs$markers[k])}))
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
  columns <- unique(sapply(4:ncol(cndata), function(x){strsplit(colnames(cndata)[x], split='[.]')[[1]][1]}))
  samn <- which(columns==sam)
  lrr.dat <- cndata[, samn*3+2]
  baf.dat <- cndata[, samn*3+3]
  par(mfrow=c(2, 2))
  #maintab <- paste('Chr ', chr, sep='')
  plot(1:length(lrr), lrr, xlab='SNP Index', ylab='LRR', ylim=c(-1.5, 1), pch='.', col='gray', 
       main='A) SNP Array Data (simulated data)')
  plot(1:length(lrr.dat), lrr.dat, xlab='SNP Index', ylab='LRR', ylim=c(-1.5, 1), pch='.', col='gray', 
       main='B) SNP Array Data (real data)')
  abline(h=c(log10((1:4)/2)), col=c('gray', 'black', 'gray', 'gray'))
  abline(h=c(log10((1:4)/2)), col=c('gray', 'black', 'gray', 'gray'))
  plot(1:length(baf), baf, main='', xlab='SNP Index', ylab='BAF', ylim=c(0, 1), pch='.', col='gray')
  abline(h=.5, col='black')
  plot(1:length(baf.dat), baf.dat, main='', xlab='SNP Index', ylab='BAF', ylim=c(0, 1), pch='.', col='gray')
  abline(h=.5, col='black')
  #plot(density(baf), xlim=c(-1.5, 1.5), ylim=c(0, 5), main=paste('SD=', round(sd(baf), digits=3), sep=''))
}

plotSEQ <- function(seqdata, seqsim){
  som.dat <- seqdata[seqdata$somatic_status=='SOMATIC', ]
  vaf.dat <- na.omit(som.dat$alt_read_count/(som.dat$alt_read_count + som.dat$ref_read_count))
  tc.dat <- som.dat$alt_read_count + som.dat$ref_read_count
  seq.sim <- seqsim$dat$seq.data
  som.sim <- seq.sim[seq.sim$status=='somatic', ]
  tc.sim <- som.sim$totalCounts
  vaf.sim <- som.sim$VAF
  par(mfrow=c(2, 2))
  plot(tc.sim[1:225], ylim=c(0, 250), xlab='Index', ylab='Total Read Counts', main='A) Somatic Variants (Simulation)', pch=16)
  plot(tc.dat[1:225], ylim=c(0, 250), xlab='Index', ylab='Total Read Counts', main='B) Somatic Variants (Data)', pch=16)
  plot(vaf.sim[1:225], ylim=c(0, 1), xlab='Index', ylab='VAF', main='', pch=16)
  plot(vaf.dat[1:225], ylim=c(0, 1), xlab='Index', ylab='VAF', main='', pch=16)
}

plotSNP(chr.sim=2, j=68, homFrac=.68, mu.0=0, sigma.0=.03, chr.dat=2, sam='NA18855')
plotSEQ(seqdata=get(load('C:\\Users\\Mark\\OneDrive - The Ohio State University\\deepcnv\\GS1115.rda')), 
    seqsim=get(load('C:\\Users\\Mark\\OneDrive - The Ohio State University\\clonetools\\dat\\data-1\\dat-1-10.rda')))

