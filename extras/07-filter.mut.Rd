\name{filter.mut}
\alias{filter.mut}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Filters mutation data to remove single mutations don't contribute useful information for CloneFinder.
}
\description{
Filters segmented SNP array data to get rid of mutations that appear to be clonal and therefore contribute no information pertaining to clonal heterogeneity for CloneFinder analysis in order to cut down computation time. }
\usage{
filter.mut(data, mu, threshold=0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mutdata}{A matrix or data frame of mutations containing the following columns: }
  \item{mu}{The average number of referents counts (or variant counts) one expects to see in a clonal, heterozygous mutations}
  \item{threshold}{Determines the minimum deviation of the observed read counts from mu necessary to pass filtering (if the observed number of both reference counts and variant counts is clone to mu, it suggests the mutation is clonal rather than subclonal, and isn't informative about clonal heterogeneity)}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{'mat' }{The filtered mutation data}
\item{'ids' }{The unique ids of the mutations from the input data included in the filtered data}
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
mutdata <- get(load('path/mutdata.rda'))
filtered <- filter.mut(mutdata,mu=100,threshold=5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~CloneFinder }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~SNP Array }% __ONLY ONE__ keyword per line
