\name{StringGraph-class}
\docType{class}
\alias{StringGraph}
\alias{StringGraph-class}
\alias{plot,StringGraph,ANY-method}
\alias{MotifGraph}
\alias{DecompositionGraph}
\alias{exportSG}

\title{Class \code{"StringGraph"}}
\description{
The \code{StringGraph} class is used to represent graphs that arise from
strings reprsenting long-read breakpoint sequences. The basic examples
are: (1) "Motif Graphs" where the edges are subtring relations, and (2)
"Decomposition Graphs" where the edges are restricted subtring relations
that decompose a long read.
}
\section{Objects from the Class}{
  Objects should be defined using the \code{MotifGraph} or
  \code{DecompositionGraph} constructor. You typically pass in a
  "motifNodes" object, which is a list of sequence-strings separated by
  length, along with some auxiliary information.
}
\usage{
MotifGraph(motifNodes, alfa, name = "motif")
DecompositionGraph(decomp, alfa, motifNodes, name = "decomp")
exportSG(sg)
\S4method{plot}{StringGraph,ANY}(x, y, \dots)
}
\arguments{
  \item{motifNodes}{A list of node names and counts, separated by
    length. In particular, \code{motifNodes[[L]]} should contain the
    nodes of length L.}
  \item{alfa}{A \code{Cipher} object.}
  \item{name}{A character vector of length one.}
  \item{decomp}{A decomposition object; see details.}
  \item{sg}{An object of the \code{StringGraph} class.}
  \item{x}{An object of the \code{StringGraph} class.}
  \item{y}{Anything; it is ignored.}
  \item{\dots}{Extra graphical parameters.}
}
\section{Slots}{
  \describe{
      \item{\code{name}:}{A character vector of length one.}
      \item{\code{edgelist}:}{A matrix representing a graph as a list of edges.}
      \item{\code{nodelist}:}{A matrix representing the nodes of the
	graph, along with their properties.}
      \item{\code{graph}:}{An \code{igraph} object.}
      \item{\code{layout}:}{A matrix containing x-y locations for the nodes.}
  }
}
\value{
The \code{MotifGraph} and \code{DecompositionGraph} functions return a
new object of the \code{StringGraph} class. The \code{plot} method and
\code{exportSG} functions return nothing and are called for their side
effects.
}
\author{Kevin R. Coombes <krc@silicovore.com>}
\note{
Attempting to manipulate a \code{StringGraph} object using text containing
NAs, missing values, or previously unknown symbols will result in an error.
}
\examples{
}
\keyword{ math }
