

makeWords <- function(opstrings, K) {
  temp <- sapply(opstrings, function(use, K) {
    pop <- strsplit(use, "")[[1]]
    L <- length(pop)
    if (L < K) return(NULL)
    if (L == K) {
      val <- 1
      names(val) <- use
      return(use)
    }
    sapply(0:(L-K), function(S) paste0(pop[S + (1:K)], collapse=""))
  }, K = K)
  table(unlist(temp)) # this is the slow step. look here if you need to speed it up
}
countWords <- function(opstrings, K, alpha = NULL) {
  m <- makeWords(opstrings, K)
  if (!is.null(alpha)) {
    names(m) <- decode(alpha, names(m))
  }
  m
}

plotWords <- function(k, m) {
  V <- as.vector(m[[k]])
  N <- rownames(as.matrix(m[[k]]))
  oo <- order(V)
  V <- V[oo]
  N <- N[oo]
  L <- length(V)
  if (L > 50) {
    keep <- (L-50):L
    V <- V[keep]
    N <- N[keep]
    L <- length(V)
  }
  plot(1:L, V, type = "n", ylim = c(-5, 5 + max(V)),
       xlab = "Index", ylab = "Count", main = paste(k, "mers", sep = "-"))
  text(1:L, V, N, srt = 90)
  invisible(V)
}
