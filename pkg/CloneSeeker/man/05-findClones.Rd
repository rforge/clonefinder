\name{Seeking Clones}
\alias{seekClones}
\alias{runAlg}
\title{Seeking Tumor Clones From Data}
\description{
  Starting with copy number segmentation data and/or sequencing mutation
  data for a tumor, seek the number of clones, the fraction of cells
  belonging to each clone, and the likely set of abnormalities in each
  clone.
}
\usage{
seekClones(cndata, vardata, cnmodels, psiset, pars, imputedCN = NULL)
runAlg(...)
}
\arguments{
  \item{cndata}{ A data frame with seven columns; can also be NULL. The
    names of the required columns are enumerated in the  man page for
    \code{generateTumorData}. These are the same as the output typically
    produced by the \code{DNAcopy} algorithm, where each row represents a
    segment (a contiguous region of a chromosome) where the copy number is
    believed to be constant.}
  \item{vardata}{ A data frame with eight columns; can also be NULL. The
    names of the required columns are enumerated in the man page for
    \code{generateTumorData}. These are typical of the output of a DNA
    sequencing experiment that has been processed to identify variants,
    which may be either germline or somatic.}
  \item{cnmodels}{ A matrix. Each row represents a model to be considered;
    each column represents a clone.  The entries are integers specifying
    the number of DNA copies present in that clone. See details.}
  \item{psiset}{ A matrix. Each column represents a clone, and each row
    represents a different possible model of the fraction of cells per
    clone. See details. }
  \item{pars}{ A list of algorithm parameters; see details.}
  \item{imputedCN}{a logical value; if missing, should the copy number be
    imputed from the mutation data.}
  \item{...}{additional variables}
}
\details{
  The algorithm starts with an initial set of 'psi' parameters
  (representing the fraction of tumor belonging to each clone). It
  computes the best (maximum a posteriori) clonal copy number and/or
  number of mutated alleles for each clone for each segment/mutation,
  conditional on the data and each of the initial psi vectors. It then
  computes the posterior probability for each psi-vector and its
  computed copy number and mutation parameters. It uses these posterior
  probabilities to resample new possible psi-vectors. The process
  repeats iteratively, and with each iterations obtains a better
  estimate of psi and the clonal segment copy number and mutation
  assignments until it terminates.
   
  The set of copy number models that we use is typically
  generated using the following command:
  \code{as.matrix(expand.grid(lapply(1:5, function(i){0:5})))}
  This setup considers all (7776) possible models with up to five clones,
  where the copy number for each clone ranges from 0 to 5. (In the future,
  we are likely to make this the default; right now, you have to generate
  these models yourself.)

  The set of possible psi-vectors (that is, the fraction of cells allocated
  to each clone) that we use is typically generated using the following
  command:
  \code{psis.20 <- generateSimplex(20,5)}
  This setup considers all (192) possible divisions of the tumor into up to
  five clones, where the fraction of cells per clone is any possible
  multiple of 0.05. Each row is sorted to put the most abundant clones first,
  which makes it easier to identify specific clones, except in the rare case
  when two clones contain exactly the same fraction of cells. (In the future,
  we are likely to make this the default; right now, you have to generate
  these models yourself.)

  The object \code{pars} is a list of numerical algorithm parameters. The
  elements are:

  \describe{
    \item{sigma0}{The standard deviation of measured allelic copy
      number at the SNP level.}
    \item{ktheta}{The probability parameter of the geometric prior
      distribution on \code{K}, the number of clones.} 
    \item{theta}{The probability parameter of the geometric prior
      distribution on genomic copy number.}
    \item{mtheta}{The probability parameter of the geometric prior
      distribution on the occurence of point mutations.}
    \item{alpha}{The (repeated) alpha parameter of a symmetric Dirichlet
      distributed prior on the  fractions of cells belong to each clone;
      default value is 0.5, giving a Jeffreys Prior.}
    \item{thresh}{The threshold determining the smallest possible
      detectable clone.}
    \item{cutoff}{SNP array segments with fewer markers than this are
      excluded.} 
    \item{Q}{Determines the number of new \code{psi} vectors resampled
      from the estimated posterior probability distribution at each
      iteration of the algorithm } 
    \item{iters}{The number of iterations in the algorithm.}
  }
  The default settings we used are from
  commonly used unfinformative priors (e.g., alpha=0.5 for the Dirichlet
  distribution is the Jeffreys Prior) or based on empirical assessments
  of the variation in data (sigma0, for example, which describes
  variation in SNP array data).
  
  Note that \code{runAlg} (an alias for \code{seekClones}) is DEPRECATED.
}
\value{
  The \code{seekClones} function returns a (rather long) list containing:

  \item{psi}{The most likely posterior psi-vector, given the data. The
    number of non-zero entries is the number of clones found, and the non-zero
    entries are the fraction of cells per clone}
  \item{A}{The most likely copy numbers for the A allele in each segment in each clone.}
  \item{B}{The most likely copy numbers for the B allele in each segment in each clone.}
  \item{psibank}{A matrix, where each row is one of the psi-vectors
    considered during the analysis.}
  \item{psiPosts}{A numeric vector, the (marginal) posterior probability of
    each psi-vector considered during the analysis.}
  \item{indices}{???}
  \item{data}{a list with two data-frame components containing the data
    used during the analysis.}
  \item{filtered.data}{a list with two data-frame components containing the
    filtered data used during the analysis. Filtering removes
    non-informative segments that have normal copy number or contain only
    germline mutations.}
  \item{etaA}{A vector of the weighted average allelic copy number for the 'A-Allele' at each segment (that is, the sum of the clonal A-allelic copy number values multiplied by the fraction of the tumor made up by each clone)}
  \item{etaB}{A vector of the weighted average allelic copy number for the 'B-Allele' at each segment}
  \item{etaM}{A vector of the weighted average number of copies of the mutated allele at each mutation}
  \item{mutated}{A matrix of the number of mutated alleles at each locus in each clone, where the number of rows is the number of somatic mutations in the data and the number of columns is the number of clones}
}
\author{
  Kevin R. Coombes \email{krc@silicovore.com},
  Mark Zucker \email{zucker.64@buckeyemail.osu.edu}
}
\examples{
# set up models
psis.20 <- generateSimplex(20,5)
cnmodels <- as.matrix(expand.grid(lapply(1:5, function(i){ 0:5 })))
# set up algortihm parameters
pars <- list(sigma0=5, theta = 0.9, ktheta = 0.3, mtheta = 0.9,
             alpha = 0.5, thresh = 0.04, cutoff = 100, Q = 100, iters = 4)
# create a tumor
psis <- c(0.6, 0.3, 0.1) # three clones
tumor <- Tumor(psis, rounds = 100, nu = 0, pcnv = 1, norm.contam = FALSE)
# simulate a dataset
dataset <- generateTumorData(tumor, 10000, 600000, 70, 25, 0.15, 0.03, 0.1)
\donttest{result <- seekClones(dset$cn.data, dset$seq.data,
             cnmodels, psis.20, pars = pars, imputedCN = NULL) }
}
\keyword{array}
