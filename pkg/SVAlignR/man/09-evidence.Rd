\name{Circular DNA Evidence}
\alias{getEvidence}
\alias{showCoverage}

\title{Support and Evidence for Circular DNA}
\description{
  Tools to assess the evidence supporting the idea that a long-read
  breakpoint sequence actually represents a piece of extrachromosomal
  circular DNA.
}
\usage{
getEvidence(cycles,  seqs, alfa, totalOnly = TRUE)
showCoverage(ID, evidence, alfa, colsams, locsams = NULL)
}
\arguments{
  \item{cycles}{A character vector containg the breakpoint strings of
    the potential circles to be evaluated.}
  \item{seqs}{A chaarcter vector continaing the breakpoint sequences of
    the long read seqeunces to be used for support of circles.}
  \item{alfa}{An object of the \code{\link{Cipher}} class.}
  \item{totalOnly}{A logical value; only return cycles that are totally
    supported?}
  \item{ID}{An integer index into the "coverage" list indicating which
    cycle to display.}
  \item{evidence}{A list as output by the \code{getEvidence} function.}
  \item{colsams}{A character vector indicating which colors to use for
    supporting long reads. If missing, all long reads are colored grey.}
  \item{locsams}{A character vector indicating which samples contained
    the long reads.}
}
\details{
TDA
}
\value{
  The \code{getEvidence} function returns a list with three entries:
  \itemize{
    \item{\code{coverage} is itself a list with one entry for each
    cycle. These entries contain \code{N} (an integer number of
    supporting reads), \code{RDS} (a character string with long-read
    identifiers seaparated by semicolons), \code{TOT} (an integer, the
    minimum support level for the nodes/brealpoints in the cycle),
    \code{WTOT} (an integer, the minimum weighted count of support
    across nodes), \code{COV} (an intege vector showing the number of
    reads supporting each individual node, and \code{WCOV} the weighted
    version of \code{COV}).}
    \item{\code{circles} is a data frame with one row for each cycles and
    columns \code{N}, \code{RDS}, \code{TOT}, \code{WTOT}, and
    \code{cycle}. If \code{totalOnly} was \code{TRUE}, only totally
    supported cycles (those with \code{TOT>0}) are returned.}
    \item{\code{vspan} a character vector containng the strings
    describing an expanded set of potential cycles.}
  }
}
\author{Kevin R. Coombes <krc@silicovore.com>}
\references{
  TDA
}
\seealso{
  TDA
}
\examples{
getEvidence
}
\keyword{hplot}
\keyword{models}
