\name{DeBruijn-class}
\docType{class}
\alias{deBruijn}
\alias{DeBruijn-class}

\title{Class \code{"DeBruijn"}}
\description{
  Classes for contructing de Bruijn graphs from collections of long read
  sequences mappe over brakpoints.
}
\section{Objects from the Class}{
  Objects should be defined using the \code{deBruijn} constructor.
}
\usage{
deBruijn(rawseq, M)
}
\arguments{
  \item{rawseq}{A character vecvtor of the long read sequences,
    expressed as hyphen-separated breakpoint ids.}
  \item{M}{An integer; the length of the motifs/words to be used in
    constructing the graph.}
}
\section{Slots}{
  \describe{
      \item{\code{G}:}{An object of hhe \code{\link[igraph]{igraph}} class}
      \item{\code{adjmat}:}{An adjacency matrix.}
      \item{\code{motifs}:}{A table of motifs/words.}
  }
}
\value{
The \code{deBruijn} constructor returns a newly created object of the
\code{DeBruijn} class.
}
\author{Kevin R. Coombes <krc@silicovore.com>}
\keyword{ manip }
