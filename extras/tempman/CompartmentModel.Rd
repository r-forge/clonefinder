\name{CompartmentModel}
\alias{CompartmentModel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
CompartmentModel(markers, pureCenters, sigma0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{markers}{
%%     ~~Describe \code{markers} here~~
}
  \item{pureCenters}{
%%     ~~Describe \code{pureCenters} here~~
}
  \item{sigma0}{
%%     ~~Describe \code{sigma0} here~~
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
function (markers, pureCenters, sigma0) 
{
    if (length(markers) == 1) {
        markers <- round(runif(markers, 25, 2000))
    }
    nSegments <- length(markers)
    if (is.null(names(markers))) 
        names(markers) <- paste("Segment", 1:nSegments, sep = "")
    new("CompartmentModel", markers = markers, pureCenters = pureCenters, 
        sigma0 = sigma0)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
