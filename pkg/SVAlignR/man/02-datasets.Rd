\name{SVAkignR-data}
\alias{SVAlignR-data}
\alias{longreads}
\docType{data}
\title{SVAlignR Sample Data}
\description{
  These data sets contain binary versions of subsets of cytogenetic
  karyotype data from patients with chronic myelogenous leukemia (CML). 
}
\usage{
data("longreads")
}
\format{
  \describe{
    \item{\code{longreads}}{A data frame with 197 rows and 5 columns.
      Each row rpresents a single Oxford Nanopore long read frema study
      of a cell line from an HPV-positive head-and-neck squamous cell
      tumor. The five columns contain (i) a uniqie identifier of each
      long read, (ii) the length of the read, in bytes, (iii) the
      ordered sequence of break point junctions, represented as a hyphen
      separated list of numeric identifiers, (iv) manually estimated
      natural groups of reads, and (v) a manul indication of wheher
      certain long reds shouold be omitted from the analysis.
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

