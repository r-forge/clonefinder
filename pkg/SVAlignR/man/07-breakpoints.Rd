\name{Breakpoints-class}
\docType{class}
\alias{Breakpoints}
\alias{Breakpoints-class}
\alias{plot,Breakpoints,missing-method}

\title{Class \code{"Breakpoints"}}
\description{
  Classes for working with collections of breakpoints.
}
\section{Objects from the Class}{
  Objects should be defined using the \code{Breakpoints} constructor. You
  typically pass in a data frame containing columns with the name/id of
  the breakpoint, and their chromosome name, start, and stop positions for
  each side of the break.
}
\usage{
Breakpoints(working)
\S4method{plot}{Breakpoints,missing}(x, y, \dots)
}
\arguments{
  \item{working}{A data frame containing the locations of break
    points. These should be seven consecutive columns, starting with the
    break point id followed by three columns each (chromosome, start,
    stop) for each side of the break point.}
  \item{x}{An object of the \code{Breakpoints} class.}
  \item{y}{Anything; it is ignored.}
  \item{\dots}{Extra graphical parameters.}
}
\section{Slots}{
  \describe{
      \item{\code{relLocation}:}{A numeric vector giving relative
	coordinates (in the unit interval) of the breakpoints along a
	chromosome, with first and last break points mapped to 0 and 1.}
      \item{\code{labels}:}{A character vector containing the names
	of the chromosomes.}
      \item{\code{ypos}:}{A numeric vector indicating the chromosomes
	involved in the full set of break points.}
      \item{\code{spread}:}{How far apart the display of different
	chromosomes should be spread apart on the y-axis.}
      \item{\code{id}:}{The character vector of break point names.}
  }
}
\value{
The \code{Breakpoints} constructor returna a newly created object of the
\code{Breakpoints} class. The \code{plot} method invisible returns its
first argument.
}
\author{Kevin R. Coombes <krc@silicovore.com>}
\keyword{ manip }
