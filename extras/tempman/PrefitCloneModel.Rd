\name{PrefitCloneModel}
\alias{PrefitCloneModel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
PrefitCloneModel(segmentdata, nPhi = 10000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{segmentdata}{
%%     ~~Describe \code{segmentdata} here~~
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
function (segmentdata, nPhi = 10000) 
{
    phiset <- .reorderVectors(sampleSimplex(nPhi, 5))
    likelihoods <- .computeLikelihoods(phiset, segmentdata, tumor, 
        TRUE)
    maxLikeIndex <- apply(likelihoods, 1, which.max)
    phipick <- phiset[maxLikeIndex, ]
    new("PrefitCloneModel", data = segmentdata, phiset = phiset, 
        likelihoods = likelihoods, maxLikeIndex = maxLikeIndex, 
        phipick = phipick)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
