library(CloneSeeker)

dataPars <- list(snps.seq = 1000000,
                 snps.cgh = 600000,
                 mu = 70,
                 sigma.reads = 25,
                 sigma0.lrr = 0.15, 
                 sigma0.baf = 0.03,
                 density.sigma = 0.1)
psis <- c(0.6, 0.3, 0.1)

RNGversion("3.5.3")
set.seed(412634)
### Mutations only? without contamination
tumor <- Tumor(psis, rounds = 400, nu = 100, pcnv = 0, norm.contam = FALSE)
clone <- getClone(tumor, 1)
summary(clone$cn)  # why is the parent index constant?
summary(clone$seq)

### Mutations only? with normal contamination
tumor <- Tumor(psis, rounds = 400, nu = 100, pcnv = 0, norm.contam = TRUE)
clone <- getClone(tumor, 1)
summary(clone$cn)  # why is the parent index missing?
summary(clone$seq) # this looks like a bug to me.

### CNV and Mutations, without normal contamination
tumor <- Tumor(psis, rounds = 400, nu = 100, pcnv = 0.5, norm.contam = FALSE)
clone <- getClone(tumor, 1)
summary(clone$cn)
summary(clone$seq)

### CNV and Mutations, with normal contamination
tumor <- Tumor(psis, rounds = 400, nu = 100, pcnv = 0.5, norm.contam = TRUE)
clone <- getClone(tumor, 1)
summary(clone$cn)
summary(clone$seq) # this looks like a bug to me.

### CNV-only, with normal contamination
tumor <- Tumor(psis, rounds = 400, nu = 0, pcnv = 1, norm.contam = TRUE) 
clone <- getClone(tumor, 1)
summary(clone$cn)
summary(clone$seq) # this is correct

### CNV-only, without normal contamination
tumor <- Tumor(psis, rounds = 400, nu = 0, pcnv = 1, norm.contam = FALSE)
clone <- getClone(tumor, 1)
summary(clone$cn)
summary(clone$seq) # this is correct

### Test coercion routines
tumor <- as(tumor, "list")
object <- as(tumor, "Tumor")

### Test data generation (i.e., simulation)
dataset <- generateTumorData(object,
                dataPars$snps.seq, dataPars$snps.cgh, dataPars$mu,
                dataPars$sigma.reads, dataPars$sigma0.lrr,
                dataPars$sigma0.baf, dataPars$density.sigma)
class(dataset)
length(dataset)
names(dataset)
lapply(dataset, class)
lapply(dataset, dim)
summary(dataset$cn.data)
summary(dataset$seq.data) # note that all mut.id are NA's since there were no mutations

plotTumorData(object, dataset)

