

computeSequenceDistance <- function() {
  method <- "nn" # options are "nn" or "leven"

mypar <- defaultNeedleParams
mypar$MATCH <- 2
mypar$MISMATCH <- -6
mypar$GAP <- -2

mat <- matrix(NA, nrow = length(pseudo), ncol = length(pseudo))
rownames(mat) <- colnames(mat) <- names(pseudo)
for (i in 1:length(pseudo)) {
  if (method == "leven") {
    scores <- NA
    mat[i,] <-  adist(pseudo[i], pseudo)
  } else {
    scores <- needleScores(pseudo[i], pseudo, mypar)
    mat[i, ] <- (scores - min(scores))/(max(scores) - min(scores))
    }
}
if (method == "leven") {
  symm <- as.dist(mat)
} else {
  symm <- as.dist(1 - (mat + t(mat))/2)
}
S <- as.matrix(symm)

}


### Inputs:
### 'sequences' is a character vector representing sequences
### in an arbitrary alphabet with at most 26 distinct characters
###
### Values = a list with three components
###    babel = original sequences translated into an AAStringSet
###    aligned = results of an alignment algorithm
###    alignedOriginal = results translated back to starting alphabet
align <- function(sequences, mysub = MYSUB, gapO = 10, gapE = 0.2) {
  U2 <- sort(unique(unlist(strsplit(sequences, ""))))
  if (length(U2) > 26) {
    stop("Cannot handle an alphabet with more than 26 letter!")
  }
  ## Translate from the starting alphabet to amino acids
  rewrite <- AA_ALPHABET[1:length(U2)]
  names(rewrite) <- U2
  tempaa <- sapply(strsplit(sequences, ""), function(X) {
    paste(rewrite[X], collapse = "")
  })
  babel <- AAStringSet(tempaa) # pretend we have proteins
  ## Align the sequences with a custom matrix
  aligned <- msa(babel, method = "ClustalW",
                 substitutionMatrix = mysub, gapOpening = gapO, gapExtension = gapE)
  ## Translate back to the original alphabet
  backme <- names(rewrite)
  names(backme) <- rewrite
  backme <- c(backme, "-" = "-", "?" = "?")
  reverse <- sapply(strsplit(as.character(aligned), ""), function(X) {
    paste(backme[X], collapse = "")
  })
  
  cons <- msaConsensusSequence(aligned)
  rcons <- strsplit(msaConsensusSequence(aligned), "")
  xcons <- xlate(rcons, backme, "")

  list(babel = babel,
       aligned = aligned,
       alignedOriginal = reverse, 
       cons = xcons)
}

plotCluster <- function() {

hcd %>% set("labels_col", value = myColorSet[1:NK], k = NK) %>% 
  set("branches_k_color", value = myColorSet[1:NK], k = NK) %>%
  set("branches_lwd", 2) %>% 
  set("labels_cex", 0.7) %>%
  plot(main = "Colored clusters")

hcd %>% set("labels_col", value = myColorSet[1:NK], k = NK) %>% 
  set("branches_k_color", value = myColorSet[1:NK], k = NK) %>%
  set("branches_lwd", 2) %>% 
  set("labels_cex", 0.7) %>%
  plot(main = "Colored clusters")

 mycut <- mean(rev(hc$height)[NK + -1:0])
dend2 <- cut(hcd, h = mycut)
dend2$upper %>% set("labels_col", value = myColorSet[1:NK], k = NK) %>% 
  set("branches_k_color", value = myColorSet[1:NK], k = NK) %>%
  set("branches_lwd", 2) %>% 
  set("labels_cex", 1.3) %>%
  plot(main = "Colored clusters")

plot(as.phylo(hc), type = "unrooted", 
     edge.color = "steelblue", edge.width = 2,
     tip.color = myColorSet[typer], cex = 0.7)

}

heatmap <- function() {
  heatmap(S, ColSideColors = myColorSet[typer],
        Rowv = as.dendrogram(hc),
        Colv = as.dendrogram(hc), col = viridis(64), scale = "none")

}

aliognBranch <- function() {
}

alignBranches <- function()  {
}

computeConsensus <- function() {
}

revert <- function(stuff) {
  A <- stuff$alignedOriginal
  xcons <- xlate(strsplit(stuff$cons, ""), omega, "-")
  temp <- as.character(A)
  names(temp) <- names(A)
  back <- xlate(strsplit(temp, ""), omega, "-")
  rack <- strsplit(back, "-")
  list(meet = as.matrix(as.data.frame( rack )), cons = xcons)
}

showme <- function(twoparts, col = "black", cex = 1) {
  meet <- twoparts$meet
  cons <- unlist(strsplit(twoparts$cons, "-"))
  matr <- matrix(1, nrow = nrow(meet), ncol = ncol(meet))
  matr[cons == "?", ] <- 2
  matr[!(cons %in% c("?", ":"))] <- 3
  image(1:nrow(meet), 1:ncol(meet),
        matr, col = palette36[c(2, 24, 28)],
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
  invisible(meet)
}

analyzeMotif <- function() {
}
