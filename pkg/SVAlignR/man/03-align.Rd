\name{AlignedCluster-class}
\docType{class}
\alias{AlignedCluster-class}
\alias{alignCluster}
\alias{alignAllClusters}
\alias{makeSubsMatrix}
\alias{image,AlignedCluster-method}

\title{Class \code{"AlignedCluster"}}
\description{
The \code{AlignedCluster} class is used to align a set of clustered
sequences. The \code{alignClusters} function creates a new object of the
\code{AlignedCluster} class. The \code{alignAllClusters} function takes
a \code{SequenceCluster} object and returns a list of
\code{AlignedCluster} objects. Clustering is performed using the
ClustalW algorithm. The associated class and fuinctions take care of
encoding and decoding sequences into a form that can be used by the
implementation of ClustalW in the \code{msa} package.
}
\section{Objects from the Class}{
  Objects should be defined using the \code{alignCluster} or
  \code{alignAllCluster} functions. You typically pass in a character
  vector of sequences that have already been found to forma  cluster.
}
\usage{
alignCluster(sequences, mysub = NULL, gapO = 10, gapE = 0.2)
alignAllClusters(sc, mysub = NULL, gapO = 10, gapE = 0.2)
makeSubsMatrix(match = 5, mismatch = -2)
\S4method{image}{AlignedCluster}(x, col = "black", cex = 1, main = "", ...)
}
\arguments{
  \item{sequences}{A character vector that contains all sequences
    to be aligned.}
  \item{mysub}{A square (usually symmetric) substitution matrix.}
  \item{gapO}{A numeric value defining the penalty for opening a gap.}
  \item{gapE}{A numeric value defining the penalty for extending a gap.}
  \item{sc}{An object of the \code{SequenceCluster} class.}
  \item{match}{A numeric value defining the reward for matching symbols
    from two sequences.}
  \item{mismatch}{A numeric value defining the penalty for mismatching
    symbols from two sequences.}
  \item{x}{An object of the \code{AlignedCluster} class.}
  \item{col}{A character setting the color of annotations in the image.}
  \item{main}{Character; the plot title.}
  \item{cex}{Numeric; size of teh text inside the image of the alignment
    matrix.}
  \item{\dots}{Extra arguments for generic or plotting routines.}
}
\section{Slots}{
  \describe{
      \item{\code{alignment}:}{A matrix of aligned sequences; rows are
	sequences and columns are aligned positions..}
      \item{weights}{A numeric vactor; the numbof times each unique raw
	sequence occurs.}
      \item{\code{consensus}:}{A character vector; the consensus
	sequence of a successful alignment.}
  }
}
\value{
The \code{alignCluster} function returns a new object of the \code{AlignedCluster}
class. The \code{alignAllClusters} function returns a list of
\code{AlignedCluster} objects. The \code{makeSubMatrix} function returns
a symmetric substitution matrix.
}
\author{Kevin R. Coombes <krc@silicovore.com>}
\examples{
data(longreads)
seqs <- longreads$connection[1:15]
pad <- c(rep("0", 9), rep("", 6))
names(seqs) <- paste("LR", pad, 1:length(seqs), sep = "")
seqs <- seqs[!duplicated(seqs)]
mysub <- makeSubsMatrix(match = 2, mismatch = -6)
ab <- alignCluster(seqs, mysub)
image(ab)
}
\keyword{ cluster }
