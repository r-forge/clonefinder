\name{Cipher-class}
\docType{class}
\alias{Cipher}
\alias{Cipher-class}
\alias{encode}
\alias{decode}

\title{Class \code{"Cipher"}}
\description{
The \code{Cipher} class is used to change between different alphabets
(and so behaves as a simple substitution cipher). The \code{Cipher} function
creates a new object of the \code{Cipher} class.
}
\section{Objects from the Class}{
  Objects should be defined using the \code{Cipher} constructor. You
  typically pass in a character vector of "words" that contain all the
  symbols that are contained in the text to be translated (i.e., encoded
  and decoded) between languages. A standard target alphabet is created
  along with forward and reverse transliteration rules.
}
\usage{
Cipher(sampleText, split = "-", extras = c("-" = ":", "?" = "?"))
encode(cipher, text)
decode(cipher, text)
}
\arguments{
  \item{sampleText}{A character vector that contains all symbols you
    want to be able to transliterate. Duplicate symbols are
    automatically removed.}
  \item{split}{A single character used to split words into
    symbols. Defaults to a hyphen for our applications.}
  \item{extras}{Additional characters to be added for reverse
    tranlsiteration, since they may appear as the results of
    alignments in consensus sequences.}
  \item{cipher}{An object of the \code{Cipher} class.}
  \item{text}{A character vector of words to be transliterated.}
}
\section{Slots}{
  \describe{
      \item{\code{forward}:}{A named character vector.}
      \item{\code{reverse}:}{A named character vector.}
  }
}
\value{
The \code{Cipher} function returns a new object of the \code{Cipher}
class. The \code{encode} and \code{decode} functions return character
vectors that are the same size as their input \code{text} parameters.
}
\author{Kevin R. Coombes <krc@silicovore.com>}
\note{
Attempting to manipulate a \code{Cipher} object using text containing
NAs, missing values, or previously unknown symbols will result in an error.
}
\examples{
motif <- "0-50-74-0-50-74-25-26-35"
alfa <- Cipher(motif)
alfa
en <-encode(alfa, motif)
en
de <- decode(alfa, en)
de
}
\keyword{ math }
