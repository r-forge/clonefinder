setClass("DeBruijn",
         slots = c("G" = "igraph",
                   "adjmat" = "matrix",
                   "motifs" = "table"))

## rawseq are the raw sequences, typically separated by hyphens since
## they are represented as multicharacter symbols (often chromosome number
## plus an arbitrary alphabetic character).
##
## M = length of motifs to be used in the DeBruijn graph.
##
deBruijn <- function(rawseq, M) {
  ## Start by re-encoding in a single character alphabet.
  bethe <- Cipher(rawseq)
  rs <- encode(bethe, rawseq)
  ## Get a list of motifs (character strings) of length M
  ## present in bhe data. Note that the `countWordss` function
  ## wants to convert the element names back to the orginal alphabnet
  motifs <- countWords(rs, M, bethe)
  ## decompose each long-read breakpoint sequences (LRBPS) as a
 ## sequence of overlapping M-mers. Uses this function:
  dekomp <- function (opstrings, K)  {
    temp <- sapply(opstrings, function(use, K) {
      pop <- strsplit(use, "")[[1]]
      L <- length(pop)
      if (L < K)
        return(NULL)
      if (L == K) {
        val <- 1
        names(val) <- use
        return(use)
      }
      sapply(0:(L - K), function(S) paste0(pop[S + (1:K)], 
                                           collapse = ""))
    }, K = K)
    temp
  }
  ## Don't mess with really small sequences
  rs6 <- rs[nchar(rs) > M]
  ## actually decompose everything that ios decomposable
  layers <- sapply(rs6, function(sqn) {
    dek <- decode(bethe, dekomp(sqn, M))
    sapply(dek, function(src) {
      which(names(motifs) == src)
    })
  })
  ## Create a (directed) adjacency matrix from the decompositions
  adjmat <- matrix(0, nrow = length(motifs), ncol = length(motifs))
  rownames(adjmat) <- colnames(adjmat) <- names(motifs)
  for (LRval in names(layers)) {
    LYR <- layers[[LRval]]
    for (J in 1:(length(LYR) - 1)) {
      a <- names(LYR)[J]
      b <- names(LYR)[J + 1]
      DBG <- adjmat[a, b]
      adjmat[a, b] <- DBG + 1
#      cat(a, "\t", b, "\t", DBG, "\t", adjmat[a, b],  "\n")
      if (all(adjmat == 0)) stop("fuck")
    }
  }
  ## Create a graph from the adjacency matrix
  G <- graph_from_adjacency_matrix(adjmat, mode = "directed", weighted = TRUE)
  G <- set_edge_attr(G, "width", value = 1 + log(edge_attr(G, "weight")))
  G <- set_vertex_attr(G, "size", value = 3)

  new("DeBruijn",
      G = G,
      adjmat = adjmat,
      motifs = motifs)
}

