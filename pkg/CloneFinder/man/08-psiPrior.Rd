\name{psiPrior}
\alias{psiPrior}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Determine prior probability of psi vector for CloneFinder package}
\description{
Computes the Dirichlet prior probability of the 'psi' vector (or list) put into it.
}
\usage{
psiPrior(psi, alpha=.5, kmax=5, minim=.001)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{psi}{A numerical vector or list of values from 0 to 1 that sums to 1 (otherwise it will be divided by the sum in the function) denoting the fraction of cells beling to each clone in a heterogeneous tumor}
  \item{alpha}{The alpha parameter in the Dirichlet distribution used to compute the prior}
  \item{kmax}{The maximum number of clones assumed possible in a tumor}
  \item{minim}{This value is substituted for 0 in the psi vector, since all values supplied to a Dirichlet distribution must be greater than 0.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{The prior probability density value of the psi vector used as input}
\references{
%% ~put references to the literature/web site here ~
}
\author{Mark Zucker}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
priorProb <- psiPrior(psi=c(.7,.3),alpha=.5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~Dirichlet }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~CloneFinder }% __ONLY ONE__ keyword per line
