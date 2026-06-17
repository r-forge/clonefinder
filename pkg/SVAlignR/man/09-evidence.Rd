\name{Circular DNA Evidence}
\alias{Vexpand}
\alias(bothWays)
\alias(getSupport)
\alias(getCoverage)
\alias(showCoverage)

\title{Support and Evidence for Circular DNA}
\description{
  Tools to assess the evidnce supporting the idea that a long-read
  breakpoint sequence actually represents a piece of extrachromosomal
  circular DNA.
}
\usage{
Vexpand(cycles, seqs, alfa)
bothWays(xc, seqs)
getSupport(comb, clues)
getCoverage(xc, rawSP, seqs, doubles, alfa)
showCoverage(ID, colsams)
}
\arguments{
  \item{cycles}{}
  \item{seqs}{}
  \item{alfa}{An object of the \code{\link{Cipher}} class.}
  \item{xc}{}
  \item{comb}{}
  \item{clues}{}
  \item{rawSP}{}
  \item{doubles}{}
  \item{ID}{}
  \item{colsams}{}
}
\details{

}
\value{
  \itemize{
    \item \code{Vexpand} returns .
    \item \code{bothWays} returns a list with two elements, \code{comb}
      and \code{doubles}.
    \item \code{getSupport} returns a data frame.
    \item \code{getCoverage} returns a list.
  }
}
\author{Kevin R. Coombes <krc@silicovore.com>}
\references{
}
\seealso{

}
\examples{

\keyword{color}
