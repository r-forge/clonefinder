\name{Words}
\docType{class}
\alias{makeWords}
\alias{countWords}
\alias{plotWords}

\title{Class \code{"Words"}}
\description{
Provides the ability to find, count, and plot words of specific length
in collections of strings in any sequence language.
}
\usage{
makeWords(opstrings, K, nb = 1)
countWords(opstrings, K, alpha = NULL)
plotWords(K, m)
}
\arguments{
  \item{opstrings}{A character vector containing a set of words that
    have been encoded into an alphabet where each character uses the same
    number of bytes in the encoding.}
  \item{K}{An integer; the length of the words of interest.}
  \item{nb}{An integer; the number of bytes used to encode each character.}
  \item{alpha}{A \code{Cipher} object, used to decode the word-strings.}
  \item{m}{A list of word-counts produced by the \code{makeWords} function.}
}
\details{
  For constructing motifs, or for producing De Bruijn graphs, we need to
  be able to decompose a set of input strings into "words" of a fixed
  length. In our application, the words are derived from long-read
  sequences that cross multiple breakpoints. Each breakpoint is given a
  unique name/label, thatwhich can be of arbirtrary length in order to be
  maningful to the researchers. Using the \code{\link{Cipher}} class, we
  encode the breakpoint names into character strings of the same
  size. (In the original version of this package, we used single
  characters. That approach eventually proved to be inadequate when we
  looked at long-read data from samples with a very large number of
  breakpoints. We then extended the package to work with two-byte
  codes. This solution may eventually be extended to even longer coding
  sequences.)

  The \code{makeWords} and \code{countWords} functions take as inputs a
  vector of character strings (typically describing long-read
  sequences) that have already been encoded into fixed-byte-length
  characters.  They then find all words in those strings of a given
  fixed length. They only differ in the form of their output. The former
  function returns the word counts in their encoded form; the latter
  decodes them back to the original names (as long as you provide the
  optional appropriate Cipher argument).

  The \code{plotWords} function gives a visible representaiton of words
  of length \code{K} sorted by their frequency. The x-axis contains the
  sorted word list; the y-axis is the frequency. The idea is that one
  can quickly figure out which words are most common in the input "text".
}
\value{
  The \code{makeWords} function returns a table of words (of length
  \code{K}) along with the counts of the number of times each one was
  seen in the input strings. The \code{countWords} function returns the
  same table, but with the words decoded back to the original language.
  The \code{plotWords} function returns a vector of the word counts for
  all words of length \code{K} in the list \code{m}.
}
\author{Kevin R. Coombes <krc@silicovore.com>}
\examples{
data(longreads)             # read sample data
raw <- longreads$connection # get the raw strings
alfa <- Cipher(raw)         # make a translation cipher
coded <- encode(alfa, raw)  # encode all the input strings
makeWords(coded, 3)
countWords(coded, 3, alfa)
m <- lapply(1:8, function(J) countWords(coded, J, alfa))
plotWords(3, m)
}
\keyword{ manip }
