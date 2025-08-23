
## opstrings is a character vector containing already encoded sequences
## K is the length of the words to be listed
## nb is the number of bytes used to encode each character
dekomp <- function(opstrings, K, nb = 1) {
  temp <- sapply(opstrings, function(use, K) {
    ## cat("Use:", class(use), " ", length(use), "\n", file = stderr())
    ### use is always a single character string
    pop <- strsplit(use, "")[[1]]
    ## cat("Pop:", class(pop), " ", length(use), "\n", file = stderr())
    ## cat("Sent:", pop, "\n")
    ### pop is also a single character string
    if (nb > 1) {
      tmat <- matrix(pop, ncol = nb, byrow = TRUE)
      pop <- apply(tmat, 1, paste, collapse = "")
      ## cat("After Pop:", class(pop), " ", length(use), "\n", file = stderr())
      ## cat("Got: ", pop, "\n", file = stderr())
    }
    L <- length(pop)
    if (L < K) return(NULL)
    if (L == K) {
      val <- 1
      names(val) <- use
      return(use)
    }
    sapply(0:(L-K), function(S) paste0(pop[S + (1:K)], collapse=""))
  }, K = K)
  temp
}

makeWords <- function(opstrings, K, nb = 1) {
  temp <- dekomp(opstrings, K, nb)
  table(unlist(temp)) # this is the slow step. look here if you need to speed it up
}

## opstrings is a character vector containing already encoded sequences
countWords <- function(opstrings, K, alpha = NULL) {
  nb <- ifelse(is.null(alpha), 1, alpha@bytes)
  m <- makeWords(opstrings, K, nb)
  if (!is.null(alpha)) {
    names(m) <- decode(alpha, names(m))
  }
  m
}

## 'm' is a list consisting of the ouput of countWords for different integers
## 'k' is the word length
plotWords <- function(K, m) {
  V <- as.vector(m[[K]])
  N <- rownames(as.matrix(m[[K]]))
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
       xlab = "Index", ylab = "Count", main = paste(K, "mers", sep = "-"))
  text(1:L, V, N, srt = 90)
  invisible(V)
}
