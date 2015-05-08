\name{updatePhiVectors}
\alias{updatePhiVectors}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
updatePhiVectors(object, nPhi = 10000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
%%     ~~Describe \code{object} here~~
}
  \item{nPhi}{
%%     ~~Describe \code{nPhi} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (object, nPhi = 10000) 
{
    if (!inherits(object, "PrefitCloneModel")) {
        object <- PrefitCloneModel(object, nPhi)
    }
    multiplier <- round(nPhi/1000)
    newphiset <- matrix(NA, ncol = 5, nrow = nPhi)
    for (i in 1:1000) {
        index <- 1 + multiplier * (i - 1)
        iset <- index:(index + multiplier - 1)
        newphiset[iset, ] <- rdirichlet(multiplier, 2 * object@phiset[object@maxLikeIndex[i], 
            ])
    }
    newphiset <- .reorderVectors(newphiset)
    likelihoods <- .computeLikelihoods(newphiset, object@data, 
        tumor, TRUE)
    maxLikeIndex <- apply(likelihoods, 1, which.max)
    phipick <- newphiset[maxLikeIndex, ]
    new("PrefitCloneModel", data = object@data, phiset = newphiset, 
        likelihoods = likelihoods, maxLikeIndex = maxLikeIndex, 
        phipick = phipick)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
