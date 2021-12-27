### Initial Pass: Clustering Sequences

setOldClass("dist")
setOldClass("hclust")

setClass("SequenceCluster",
         slots = c(
           method = "character",
           rawSequences = "character",
           distance = "dist",
           hc = "hclust",
           NC = "numeric",
           clusters = "numeric")
         )

SequenceCluster <- function(rawseq, method = c("needelman", "levenshtein"), NC = 5) {
  method <- match.arg(method)

  doNeedelman <- function(pseudo) {
    mypar <- defaultNeedleParams
    mypar$MATCH <- 2
    mypar$MISMATCH <- -6
    mypar$GAP <- -2
    mat <- matrix(NA, nrow = length(pseudo), ncol = length(pseudo))
    rownames(mat) <- colnames(mat) <- names(pseudo)
    for (i in 1:length(pseudo)) {
      scores <- needleScores(pseudo[i], pseudo, mypar)
      mat[i, ] <- (scores - min(scores))/(max(scores) - min(scores))
    }
    as.dist(1 - (mat + t(mat))/2)
  }

  doLevenshtein <- function(pseudo) {
    mat <- matrix(NA, nrow = length(pseudo), ncol = length(pseudo))
    rownames(mat) <- colnames(mat) <- names(pseudo)
    for (i in 1:length(pseudo)) {
      mat[i,] <-  adist(pseudo[i], pseudo)
    }
    as.dist(mat)
  }
  ## Make sure everything is named
  if (is.null(names(rawseq))) {
    warning("Adding names to input sequences.")
    L <- length(rawseq)
    nzero <- trunc(log10(1:L))
    pad <- sapply(max(nzero) - nzero, function(n) paste(rep("0", n), collapse = ""))
    names(rawseq) <- paste("SEQ", pad, 1:L, sep = "")
  }

  ## dedup
  if (any(dup <- duplicated(rawseq))) {
    warning("Removing ", sum(dup), " duplicated sequences.\n")
    rawseq <- rawseq[!dup]
  }

  ## encode the input sequences
  alphabet <- Cipher(rawseq)
  encoded <- encode(alphabet, rawseq)
  ## compute distance matrix as 'dist' object
  symm <- switch(method,
                 needelman = doNeedelman(encoded),
                 levenshtein = doLevenshtein(encoded),
                 stop("Method |", method, "| does not match any known methods." )
                 )
  hc <- hclust(symm, "ward.D2")
  clust <- cutree(hc, k = NC)
  new("SequenceCluster",
      method = method,
      rawSequences = rawseq,
      distance = symm,
      hc = hc,
      NC = NC,
      clusters = clust)
}

updateClusters <- function(sc, NC) {
  sc@NC <- NC
  sc@clusters <- cutree(sc@hc, k = NC)
  sc
}


myColorSet <- c("cornflowerblue", "hotpink2", "green4", "red",
                "darkorchid1", "gray16", "darkgoldenrod", "mediumorchid4",
                "coral1", "blue", "lightseagreen", "maroon1",
                "deeppink", "rosybrown", "darkolivegreen4", "skyblue4",
                "tomato4", "mediumpurple2", "hotpink4", "blue3",
                "orchid", "tan4", "midnightblue", "firebrick")

setMethod("plot", signature("SequenceCluster", "missing"),
function(x, type = "rooted", ...) {
  plotRooted <- function(x, NC = x@NC, ...) {
##    cat("rooted\n", file = stderr())
    hcd <- as.dendrogram(x@hc)
    hcd %>% set("labels_col", value = myColorSet[1:NC], k = NC) %>% 
      set("branches_k_color", value = myColorSet[1:NC], k = NC) %>%
      set("branches_lwd", 2) %>% 
      set("labels_cex", 0.7) %>%
      plot(main = "Colored clusters")
  }
  plotClipped <- function(x, NC, ...) {
 ##   cat("clipped\n", file = stderr())
    mycut <- mean(rev(x@hc$height)[NC + -1:0])
    hcd <- as.dendrogram(x@hc)
    dend2 <- cut(hcd, h = mycut)
    dend2$upper %>% set("labels_col", value = myColorSet[1:NC], k = NC) %>% 
      set("branches_k_color", value = myColorSet[1:NC], k = NC) %>%
      set("branches_lwd", 2) %>% 
      set("labels_cex", 1.3) %>%
      plot(main = "Colored clusters")
  }
  plotUnrooted <- function(x, NC, ...) {
 ##   cat("unrooted\n", file = stderr())
    typer <- cutree(x@hc, k = NC)
    plot(as.phylo(x@hc), type = "unrooted",
         edge.color = "steelblue", edge.width = 2,
         tip.color = myColorSet[typer], cex = 0.7)
  }

  NC <- x@NC
  types <- c("rooted", "clipped", "unrooted")
  type <- match.arg(type, types)
##  cat("Calling", type, "\n")
  val <- switch(type,
                rooted = plotRooted(x, NC, ...),
                clipped = plotClipped(x, NC, ...),
                unrooted = plotUnrooted(x, NC, ...),
                print("Invalid type.")
                )
  invisible(x)
})

heat <- function(x, ...) {
  S <- as.matrix(x@distance)
  hc <- x@hc
  typer <- cutree(hc, k = x@NC)
  if(requireNamespace("viridisLite")) {
    colscheme <- viridisLite::viridis(64)
  } else {
    colscheme <- topo.colors(64)
  }
  heatmap(S, ColSideColors = myColorSet[typer],
        Rowv = as.dendrogram(hc),
        Colv = as.dendrogram(hc), col = colscheme, scale = "none")
}

