\name{SVAkignR-data}
\alias{SVAlignR-data}
\alias{longreads}
\alias{breakpoints}
\docType{data}
\title{SVAlignR Sample Data}
\description{
  These data sets contain binary versions of data describing breakpoints
  and long read sequences from an HPV-positive head-and-neck cancer.
}
\usage{
data("longreads")
}
\format{
  \describe{
    \item{\code{longreads}}{A data frame with 197 rows and 5 columns.
      Each row represents a single Oxford Nanopore long read from a study
      of a cell line from an HPV-positive head-and-neck squamous cell
      tumor. The five columns contain (i) a uniqie identifier of each
      long read, (ii) the length of the read, in bytes, (iii) the
      ordered sequence of break points, represented as a hyphen
      separated list of numeric identifiers, (iv) manually estimated
      natural groups of reads, and (v) a manual curated indication of
      whether certain long reads should be omitted from the analysis.
    }
    \item{\code{breakpoints}}{A data frame with 82 rows and 11 columns.
      Each row represents a single breakpoint from a study of a cell line
      from an HPV-positive head-and-neck squamous cell tumor. The columns
      contain (1) a unique identifier that is used in the long read
      connections, (2-4) a description of the chromosomal segment to the
      left of the breakpoint, (5-7) a description of the chromosomal
      segment to the right of the breakpoint, (8-9) the orientation of the
      two chromosomal segments, (10) a shorthand description of the
      breakpoint with the segment names separated by a vertical bar and
      negative strands contained in parentheses, and (11) a shorthand
      representation of the reverse orientation of the breakpoint.
    }
  }
}
\source{
  Long read (Oxford Nanopore) sequencing was performed on samples
  prepared at the laboratory of Maura Gillison and David   Symer.
  Characterization of long reads as a sequence of well-defined
  break points was performed by Keiko Akagi.
}
\author{
  Kevin R. Coombes <krc@silicovore.com>
}
\examples{
data(longreads)
head(longreads)

alphabet <- Cipher(longreads$connection)
en <- encode(alphabet, "0-50-74-0-50-74-35")
en
decode(alphabet, en)
}
\keyword{datasets}

