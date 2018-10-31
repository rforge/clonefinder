\name{Finding Clones}
\alias{findClones}
\alias{runAlg}
\title{ Finding Tumor Clones From Data}
\description{
  Starting with copy number segmentation data and/or sequencing mutation
  data for a tumor, find the number of clones, the fraction of cells
  belonging to each clone, and the likely set of abnormalities in each
  clone.
}
\usage{
findClones(cndata, vardata, cnmodels, psiset, pars, imputedCN = NULL)
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
  TODO: MARK TO EXPLAIN THE BASIC ALGORITHM.
  
  The set of different copy number models that we use is typically
  generated using the following command:
  \code{as.matrix(expand.grid(lapply(1:5, function(i){0:5})))}
  This setup considers all (7776) possible models with up to five clones,
  where the copy number for each cloine ranges from 0 to 5. (In the future,
  we are likely to make this the default; right now, you have to generate
  these models yourself.)

  The set of possible psi-vectors (that is, the fraction of cells allocated
  to each clone) that we use is typically generated using the following
  command:
  \code{psis.20 <- generateSimplex(20,5)}
  This setup considers all (192) possible divisiosn of the tumor into up to
  five clones, where the fraction of cells per clone is any posisble
  multiple of 0.5. Each row is sorted to utt he most abundant clones first,
  which makes it easier to identify specific clones (except in the rare case
  when two clones contain exactly the same fraction of cells). (In the future,
  we are likely to make this the default; right now, you have to generate
  these models yourself.)

  TODO: DESCRIBE OTHER PARS.
  
  Note that \code{runAlg} (an alias for \code{findClones}) is DEPRECATED.
}
\value{
  The \code{fidClones} function returns a (rather long) list containing:

  \item{psi}{The most likely posterior psi-vector, given the data. The
    number of non-zero entris is the number of clones found, and the onzeo
    entries are the fraction of cells per clone}
  \item{A}{The most likely copy numbers for the A allele in each segment in each clone.}
  \item{B}{The most likely copy numbers for the B allele in each segment in each clone.}
  \item{psibank}{a matrix, where each row is one of the psi-vectors
    considered during the analysis}
  \item{psiPosts}{a numeric vector, the (marginal) posterior probability of
    each psi-vector considered during the analysis.}
  \item{indices}{}
  \item{data}{a list with two data-frame components containing the data
    used during the analysis.}
  \item{filtered.data}{a list with two data-frame components containing the
    filtered data used during the analysis. Filtering removes
    non-informative segments that have normal copy number or contain only
    germline mutations.}
  \item{etaA}{TODO}
  \item{etaB}{TODO}
  \item{etaM}{TODO}
  \item{mutated}{TODO}
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
\donttest{ra <- runAlg(dset$cn.data, dset$seq.data,
             cnmodels, psis.20, pars = pars, imputedCN = NULL) }
}
\keyword{array}
