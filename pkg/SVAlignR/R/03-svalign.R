makeSubsMatrix <- function(match = 5, mismatch = -2) {
  if (!is.numeric(match) | !is.numeric(mismatch)) {
    stop("Both parameters must be numeric.\n")
  }
  if (length(match) != 1 | length(mismatch) != 1) {
    stop("Both parametersm must be of length one.\n")
  }
  lab <- c(LETTERS[c(-15, -21)], "*" )
  MYSUB <- matrix(mismatch, 25, 25, dimnames = list(lab, lab))
  diag(MYSUB) <- match
  MYSUB
}

### Inputs:
### 'sequences' is a character vector representing sequences
### in an arbitrary alphabet with at most 25 distinct characters
###
### Values = a list with three components
###    babel = original sequences translated into an AAStringSet
###    aligned = results of an alignment algorithm
###    alignedOriginal = results translated back to starting alphabet
align <- function(sequences, mysub = NULL, gapO = 10, gapE = 0.2) {
  if (is.null(mysub)) mysub <- makeSubsMatrix
  U2 <- sort(unique(unlist(strsplit(sequences, ""))))
  if (length(U2) > 25) {
    stop("Cannot handle an alphabet with more than 25 letter!")
  }
  ## Translate from the starting alphabet to amino acids
  protein <- Cipher(sequences, split = "", extras = c("-" = "-", "?" = "?"))
  tempaa <- .xlate(sequences, protein@forward, "", "")
  babel <- AAStringSet(tempaa) # pretend we have proteins
  ## Align the sequences with a custom matrix
  aligned <- msa(babel, method = "ClustalW",
                 substitutionMatrix = mysub, gapOpening = gapO, gapExtension = gapE)
  ## Translate back to the original alphabet
  reverse <- .xlate(as.character(aligned), protein@reverse, "", "")
  xcons  <- .xlate(msaConsensusSequence(aligned), protein@reverse, "", "")

  list(babel = babel,
       aligned = aligned,
       alignedOriginal = reverse, 
       cons = xcons)
}

setClass("AlignedCluster",
         slots = c(
           alignment = "matrix",
           weights = "integer",
           consensus = "character")
         )

alignCluster <- function(sequences, mysub = NULL, gapO = 10, gapE = 0.2) {
  ## assume we have names on the individual sequences
  if (is.null(mysub)) mysub <- makeSubsMatrix
  seqs <- sequences[!duplicated(sequences)]  # dedup
  alfa <- Cipher(seqs)
  enc  <- encode(alfa, seqs)                 # encode as though amino acids
  sva  <- align(enc, mysub, gapO, gapE)      # align using ClustalW
  cons <- decode(alfa, sva$cons)             # decode the consensus sequence
  back <- decode(alfa, sva$alignedOriginal)  # decode the aligned sequences
  rack <- strsplit(back, "-")
  new("AlignedCluster",
      alignment = as.matrix(as.data.frame( rack )),
      weights = integer(0),
      consensus = cons)
}


alignAllClusters <- function(sc, mysub = NULL, gapO = 10, gapE = 0.2) {
  if (!inherits(sc, "SequenceCluster")) {
    stop("Bad input, sc = ", class(sc), "\n")
  }
  if (is.null(mysub)) mysub <- makeSubsMatrix
  NC <- sc@NC
  result <- lapply(1:NC, function(K) {
    ab <- alignCluster(sc@rawSequences[sc@clusters == K],
                      mysub = mysub, gapO = gapO, gapE = gapE)
    ab@weights <- sc@weights[colnames(ab@alignment)]
    ab
  })
  nzero <- trunc(log10(1:NC))
  pad <- sapply(max(nzero) - nzero, function(n) paste(rep("0", n), collapse = ""))
  names(result) <- paste("B", pad, 1:NC, sep = "")
  result
}


setMethod("image", "AlignedCluster",  function(x, col = "black", cex = 1, main = "", ...) {
  meet <- x@alignment
  wt <- x@weights
  cons <- unlist(strsplit(x@consensus, "-"))
  matr <- matrix(1, nrow = nrow(meet), ncol = ncol(meet))
  matr[cons == "?", ] <- 2
  matr[!(cons %in% c("?", ":"))] <- 3
  image(1:nrow(meet), 1:ncol(meet),
        matr, col = c("gray89", "navajowhite", "lightsteelblue1"),
        yaxt = "n", xlab = "Position", ylab = "", zlim = c(1,3))
  for (i in 1:nrow(meet)) {
    text(i, 1:ncol(meet), meet[i,], cex = cex)
  }
  tag <- apply(meet, 1, function(x) names(which.max(table(x))))
  tagcex <- ifelse(nrow(meet) > 20, 0.8, 1)
  mtext(tag, side = 3, line = 1, at = 1:nrow(matr), col = col, cex = tagcex)
  txtcex <- ifelse(ncol(meet) > 15, 
                   ifelse(ncol(meet) > 22, 0.5, 0.7), 1)
  mtext(colnames(meet), side = 2, las = 2, at = 1:ncol(meet), line = 0.5,
        col = col, cex = txtcex)
  if (length(wt) == ncol(meet)) {
    mtext(paste("n=", wt, sep = ""), side = 4, line = 1/2, at = 1:length(wt), las=2)
  }
  invisible(x)
})
