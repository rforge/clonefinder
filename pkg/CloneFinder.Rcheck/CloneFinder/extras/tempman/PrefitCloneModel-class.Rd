\name{PrefitCloneModel-class}
\Rdversion{1.1}
\docType{class}
\alias{PrefitCloneModel-class}
\alias{hist,PrefitCloneModel-method}
\alias{plot,PrefitCloneModel,missing-method}
\alias{summary,PrefitCloneModel-method}

\title{Class \code{"PrefitCloneModel"}}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("PrefitCloneModel", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{data}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{phiset}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{likelihoods}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{phipick}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{maxLikeIndex}:}{Object of class \code{"vector"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{hist}{\code{signature(x = "PrefitCloneModel")}: ... }
    \item{plot}{\code{signature(x = "PrefitCloneModel", y = "missing")}: ... }
    \item{summary}{\code{signature(object = "PrefitCloneModel")}: ... }
	 }
}
\references{
%%  ~~put references to the literature/web site here~~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
}
\examples{
showClass("PrefitCloneModel")
}
\keyword{classes}
