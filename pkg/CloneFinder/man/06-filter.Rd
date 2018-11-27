\name{filter}
\alias{filter}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Filters segmented SNP array data to remove segments don't contribute useful information for CloneFinder package.
}
\description{
Filters segmented SNP array data to get rid of segments that appear to be normal and leaving only aberrant segments useful for CloneFinder analysis in order to cut down computation time. 
}
\usage{
filter(data, threshold, cutoff=0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{}
  \item{threshold}{The minimum deviation from a normal, homogeneous copy number state necessary for a segment not to be filtered out}
  \item{cutoff}{The minimum number of markers a segment must have to pass filtering}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{'mat' }{The filtered data}
\item{'indices' }{The indices of rows included in the filtered data}
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
cndata <- get(load('path/cndata.rda'))
filtered <- filter(cndata)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~CloneFinder }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~SNP Array }% __ONLY ONE__ keyword per line
