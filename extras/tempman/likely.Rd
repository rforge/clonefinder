\name{likely}
\alias{likely}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
likely(dataset, phi, compartmentModel, log = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
%%     ~~Describe \code{dataset} here~~
}
  \item{phi}{
%%     ~~Describe \code{phi} here~~
}
  \item{compartmentModel}{
%%     ~~Describe \code{compartmentModel} here~~
}
  \item{log}{
%%     ~~Describe \code{log} here~~
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
function (dataset, phi, compartmentModel, log = FALSE) 
{
    xy <- compartmentModel@pureCenters
    markers <- compartmentModel@markers
    sigma0 <- compartmentModel@sigma0
    if (length(phi) < nrow(xy) - 1) 
        stop("You did not supply enough entries in 'phi'")
    if (length(phi) > nrow(xy)) 
        stop("You supplied too many entries in 'phi'")
    if (length(phi) < nrow(xy)) {
        lastphi <- 1 - sum(phi)
        phi <- c(phi, lastphi)
    }
    if (any(phi < 0)) 
        stop("Negative probabilities are not allowed")
    phi <- matrix(phi/sum(phi), nrow = 1)
    center <- as.data.frame(phi \%*\% as.matrix(xy))
    secondMoment <- phi \%*\% (xy^2 + sigma0^2)
    sigma <- sqrt(sweep(secondMoment - center^2, 1, markers, 
        "/"))
    px <- dnorm(dataset$x, center$x, sigma[, 1], log)
    py <- dnorm(dataset$y, center$y, sigma[, 2], log)
    if (log) {
        value <- px + py
    }
    else {
        value <- px * py
    }
    value
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
