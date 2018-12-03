\name{psiOptim}
\alias{psiOptim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Computes log posterior probabilities for 'psi' vectors, given data, and resamples new 'psi' vectors 
}
\description{
Computes posterior probabilities of input psi vectors, given copy number and/or mutation data, then uses the probabilities to resample new psi-vectors, which it returns as output in addition to the log posterior probabilties  and the most probably clonal segment copy number information associated with each input psi vector.
}
\usage{
psiOptim(cndata.filt, mutdata.filt, psis, cnmodels, pars, cnmax=5, kPriors)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cndata.filt}{Filtered copy number data, with columns: }
  \item{mutdata.filt}{Filtered mutation data, with columns:}
  \item{psis}{A matrix made up of possible of psi vectorsof length 'cnmax' (the maximum number of clones)}
  \item{cnmodels}{A matrix made up of possible clonal segment copy number state vectors of length 'cnmax'}
  \item{pars}{A list of input parameters: 'thresh' (the threshold determining the smalles possible detectable clone); 'theta' (the probability parameter of the geometric prior distribution on genomic copy number); 'ktheta' (the probability parameter of the geometric prior distribution on 'K', the number of clones); 'sigma0' (standard deviation of measured allelic copy number at the SNP level); 'mtheta' (the probability parameter of the geometric prior distribution on the occurence of point mutations)}
  \item{cnmax}{The maximum number of clones possible; default value is 5}
  \item{kPriors}{Prior probabilities of each possible (integer) number of clones; can be set to NULL, in which case, priors will be computed using a geometric distribution with parameter 'ktheta' from the list 'pars'}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{A}{Data frame of the (greater) allelic copy number state for each clone at each segment}
  \item{B}{Data frame of the (lesser) allelic copy number state for each clone at each segment}
  \item{Mutated}{Data frame of the number of mutated alleles for each clone at each mutation}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Mark Zucker
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
cndata.filt <- get(load())
mutdata.filt <- get(load())
psis <- get(load())
cnmodels <-  expand.grid()
pars <- list()
obj <- psiOptim(cndata.filt, mutdata.filt, psis, cnmodels, pars, cnmax=5, kPriors=NULL)
newPsis <- obj$psis
logProbs <- obj$posts
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~CloneFinder }% use one of  RShowDoc("KEYWORDS")
