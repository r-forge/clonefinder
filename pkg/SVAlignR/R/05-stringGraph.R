setOldClass("igraph")
setClass("StringGraph",
         slots = c(name = "character",
                   edgelist = "matrix",
                   nodelist = "matrix",
                   graph = "igraph",
                   layout = "matrix"
                   )
)

setMethod("plot", "StringGraph", function(x,y, ...) {
  plot(x@graph, layout = x@layout, ...)
})


exportSG <- function(sg) {
  if (!inherits(sg, "StringGraph")) {
    stop("Wrong object type.\n")
  }
  fileedge <- paste(sg@name, "edgelist.csv", sep = "-")
  filenode <- paste(sg@name, "nodelist.csv", sep = "-")
  write.csv(sg@edgelist, file = fileedge, row.names = FALSE)
  temp <- data.frame(sg@nodelist, sg@layout)
  write.csv(temp, file = filenode, row.names = FALSE)
}

makeEdges <- function(A, B, debug = FALSE) {
  temp <- sapply(names(A), function(nm) {
    grepl(nm, names(B))
  })
  if (debug) cat(class(temp), "\n",  file = stderr())
  if (inherits(temp, "matrix")) {
    rownames(temp) <- names(B)
  } else {
    names(temp) <- names(B)
  }
  R <- expand.grid(names(B), names(A))
  R$Edge <- as.vector(temp)
  R[R$Edge,]
}

#################### MOTIFS ##########################

.makeEdgeListFromMotifs <- function(motifNodes) {
  motifCounts <- sapply(motifNodes, length)
  W <- which(motifCounts > 0)
  top <- max(W)
  EL <- list()
  for (i in 1:(top-1)) {
    A <- motifNodes[[i]]
    if (length(A) == 0) next # no nodes of this length
    for (j in (i+1):(top)) {
      B <- motifNodes[[j]]
      if (length(B) == 0) next # no nodes of this length
      EL <- rbind(EL, makeEdges(A, B))
    }
  }
  as.matrix(EL)
}
.makeNodeListFromMotifs <- function(motifNodes, crazyColors) {
  NL <- unique(unlist(sapply(motifNodes, names)))
  crazyColors <- createPalette(length(NL),
                               seedcolors = c("#ff0000", "#00ff00"), 
                               range = c(20, 80), M = 80000)
  names(crazyColors) <- NL
  stdcol <- Polychrome:::xform(crazyColors)
  luvmat <- as(hex2RGB(stdcol), "LUV")
  x <- luvmat@coords
  labelcols <- c("white", "black")[1 + 1*(x[,1] > 50)]
  names(labelcols) <- names(stdcol)
  nodeList <- data.frame(Names = NL, 
                         Type = "Motif",
                         Color = stdcol,
                         txtColor = labelcols)
  as.matrix(nodeList)
}

.makeLayoutFromMotifs <- function(NL, alfa, motifNodes) {
  XY <- data.frame(X = NA, Y = nchar(encode(alfa, NL))) # set Y values based on length.
  nMotifs <- sapply(motifNodes, length)
  for (I in which(nMotifs > 0)) {
    W <- which(XY$Y == I)
    top <- nMotifs[I]
    if (is.na(top)) stop("abruptly")
    XY$X[W] <- ((1:top)-top/2)*123
  }
  as.matrix(XY)
  rownames(XY) <- NL
  XY
}

.createGraph <- function(edgelist, nodelist) {
  X <- as.matrix(edgelist[, 1:2])
  M <- graph_from_edgelist(X, directed = TRUE)
  N <- vertex_attr(M, "name")
  if (length(N) != nrow(nodelist)) {
    stop("Disagreement in number of nodes")
  }
  LL <- 1 + str_count(N, "\\-")
  M <- set_vertex_attr(M, "shape", value = "vrectangle")
  M <- set_vertex_attr(M, "color", value = nodelist[N, "Color"]) # get the colors
  M <- set_vertex_attr(M, "size", value = 4*LL + 3)
  M <- set_vertex_attr(M, "size2", value = 4)
  M <- set_edge_attr(M, "arrow.size", value = 0.3)
  M <- set_edge_attr(M, "weight", value = 10)
  M <- set_edge_attr(M, "color", value = "black")
  M
}

MotifGraph <- function(motifNodes, alfa, name = "motif") {
  edgelist <- .makeEdgeListFromMotifs(motifNodes)
  nodelist <- .makeNodeListFromMotifs(motifNodes, alfa)
  layout <- .makeLayoutFromMotifs(nodelist[, "Names"], alfa, motifNodes)
  graph <- .createGraph(edgelist, nodelist)
  gname <- vertex_attr(graph, "name")
  layout <- as.matrix(layout[gname,])
  new("StringGraph",
      name = name,
      edgelist = edgelist,
      nodelist = nodelist,
      layout = layout,
      graph = graph)
}


################## DECOMP ############################

.makeEdgesFromDecomposition <- function(decomp) {
  SAP <- sapply(decomp, function(L) {
    a <- paste(L, collapse = "-")
    t(t(data.frame(motif = L, longread = a)))
  })
  EL2 <- do.call(rbind, SAP)  # edgelist
  EL2
}

.makeNodesFromDecomp <- function(edgelist, motifNodes) {
  temp <- unique(c(edgelist[,1], edgelist[,2]))
  nodeList2 <- data.frame(Names = temp,
                          Type = factor(rep(c("LRBP", "Motif"), times = length(temp)/2)),
                          Color = "#C0C0C0",
                          txtColor = "black")
  rownames(nodeList2) <- nodeList2[,"Names"]
  nodelist <- .makeNodeListFromMotifs(motifNodes)
  rownames(nodelist) <- nodelist[,"Names"]
  nodeList2[rownames(nodelist), ] <- nodelist
  nodeList2
}

.makeDecompLayout <- function(NL, alfa) {
  nuke <- nchar(encode(alfa, NL))
  XY <- data.frame(X = NA, Y = (-1)*nuke) # set Y based on length.
  nMotifs <- table(nuke)
  for (I in which(nMotifs > 0)) {
    W <- which(XY$Y == -I)
    top <- nMotifs[I]
    if (is.na(top)) stop("abruptly")
    XY$X[W] <- ((1:top)-top/2)*123
  }
  XY <- as.matrix(XY)
  XY[is.na(XY)] <- 0
  XY
}

DecompositionGraph <- function(decomp, alfa, motifNodes, name = "decomp") {
  edgelist <- .makeEdgesFromDecomposition(decomp)
  nodelist <- .makeNodesFromDecomp(edgelist, motifNodes)
  layout <- .makeDecompLayout(nodelist$Names, alfa)
  rownames(layout) <- nodelist$Names
  graph <- .createGraph(edgelist, nodelist)
  VA <- vertex_attr(graph, "name")
  layout <- layout[VA,]
  new("StringGraph",
      name = name,
      edgelist = as.matrix(edgelist),
      nodelist = as.matrix(nodelist),
      layout = layout,
      graph = graph)
}

