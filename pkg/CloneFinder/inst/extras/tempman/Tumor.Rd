\name{Tumor}
\alias{Tumor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Tumor(object, fracs, weights)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
%%     ~~Describe \code{object} here~~
}
  \item{fracs}{
%%     ~~Describe \code{fracs} here~~
}
  \item{weights}{
%%     ~~Describe \code{weights} here~~
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
function (object, fracs, weights) 
{
    if (any(is.na(fracs))) 
        stop("Missing clone fraction.")
    if (any(fracs < 0)) 
        stop("Negative clone fraction.")
    fracs <- fracs/sum(fracs)
    L <- length(fracs)
    if (is.null(names(fracs))) 
        names(fracs) <- paste("Clone", 1:L, sep = "")
    if (is.null(names(weights))) 
        names(weights) <- paste("Compartment", 1:length(weights), 
            sep = "")
    nSegments <- length(object@markers)
    segnames <- names(object@markers)
    dataset <- matrix(NA, nrow = nSegments, ncol = L)
    for (i in 1:L) {
        dataset[, i] <- s <- Clone(nSegments, weights, segnames)@segments
    }
    dimnames(dataset) <- list(names(s), names(fracs))
    nCompartments <- length(weights)
    fvec <- matrix(fracs, ncol = 1)
    repr <- matrix(NA, nrow = nSegments, ncol = nCompartments)
    for (comp in 1:nCompartments) {
        temp <- 1 * (dataset == comp)
        repr[, comp] <- as.vector(temp \%*\% fvec)
    }
    dimnames(repr) <- list(rownames(dataset), names(weights))
    xy <- as.matrix(object@pureCenters)
    centers <- as.data.frame(repr \%*\% xy)
    secondMoment <- repr \%*\% (xy^2 + object@sigma0^2)
    sigma <- sqrt(sweep(secondMoment - centers^2, 1, object@markers, 
        "/"))
    new("Tumor", object, data = dataset, fraction = fracs, weights = weights, 
        compartments = repr, centers = centers, SEM = sigma)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
