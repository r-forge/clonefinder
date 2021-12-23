\name{SequenceCluster-class}
\docType{class}
\alias{SequenceCluster}
\alias{SequenceCluster-class}
\alias{plot,SequenceCluster,missing-method}

\title{Class \code{"SequenceCluster"}}
\description{
The \code{SequenceCluster} class is used to cluster seqeunces of "words"
from an arbitrarily long alphabet. The \code{SequenceCluster} function
returns a new object of the \code{SequenceCluster} class.
}
\section{Objects from the Class}{
  Objects should be defined using the \code{SequenceCluster} constructor. You
  typically pass in a character vector of "words"  to be clustered.
}
\usage{
SequenceCluster(encoded, method = c("needelman", "levenshtein"))
\S4method{plot}{SequenceCluster,missing}(x, NK = 5, type = "rooted", \dots)
}
\arguments{
  \item{encoded}{A character vector that contains all words or "sequences"
    to be clustered.}
  \item{method}{The algorithm to use to compute distances between
    sequences. The choices are "levenshstein", which uses the
    Levenshtein edit distance, or "needelman", which uses the
    Needelman-Wunsch global alignment algorithm.}
  \item{x}{An object of the \code{SequenceCluster} class.}
  \item{NK}{An integer; the number of branches/clusters to cut from the dendrogram..}
  \item{type}{A character strnig; the type of lpot to make. Val;id types
    are "rooted", "clippe", or "unrooted".}
  \item{\dots}{extra arguments for generic or plotting routines}
}
\section{Slots}{
  \describe{
      \item{\code{method}:}{A character vector describing which
	algorithm was used.}
      \item{\code{distance}:}{A \code{dist} object.}
      \item{\code{hc}:}{An \code{hcult} object.}
  }
}
\value{
The \code{SequenceCluster} function returns a new object of the \code{SequenceCluster}
class.
}
\author{Kevin R. Coombes <krc@silicovore.com>}
\examples{
data(longreads)
sequences <- longreads$connection[1:30]
sequences <- sequences[!duplicated(sequences)] # dedup
alfa <- Cipher(sequences)
enc <- encode(alfa, sequences)                 # encode
sc <- SequenceCluster(enc)                     # cluster
plot(sc@hc)                                    # visualize
}
\keyword{ cluster }
