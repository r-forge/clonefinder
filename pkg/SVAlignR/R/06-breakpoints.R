
setClass("Breakpoints",
         slots = c(relLocation = "numeric",
                   labels = "character",
                   ypos = "numeric",
                   spread = "numeric",
                   id = "character"))

Breakpoints <- function(working) {
  A <- data.frame(working[, 1:3], side = "left")
  B <- data.frame(working[, c(1, 5:6)], side = "right")
  colnames(A) <- colnames(B) <- c("id", "chrom", "base", "side")
  R <- rbind(A, B)
  R <- R[order(R$chrom, R$base),]
  U <- unique(R$chrom)
  pick <- sapply(U, function(u) R$chrom == u) 
  mini <- aggregate(R$base, list(R$chrom), min)$x
  maxi <- aggregate(R$base, list(R$chrom), max)$x
  delta <- maxi - mini
  names(maxi) <- names(mini) <- names(delta) <- U
  unitcoords <- apply(R, 1, function(arow) {
    ch <- as.character(arow[2])
    (as.numeric(arow[3]) - mini[ch])/delta[ch]
  })

  ct <- as.numeric(factor(R$chrom))
  spread <- c(0.4, 0, -0.4, 0.1, -0.3, 0.2, -0.2 ,0.3, -0.1)
  spread <- rep(spread, times = 1 + trunc(length(ct/9)))[1:length(ct)]
  new("Breakpoints",
      relLocation = unitcoords,
      labels = U,
      ypos = ct,
      spread = spread,
      id = R$id)
}

  setMethod("plot", c("Breakpoints", "missing"), function(x, y, ...) {
    N <- length(x@labels)
    plot(x@relLocation, x@ypos + x@spread,
         type = "n", ylim = c(0.5, N + 0.5), yaxt = "n",
         xlab = "Relative Base Posiiton", ylab = "Chromosome")
    mtext(x@labels, side = 2, line = 1, at = 1:N, las = 3)
    text(x@relLocation, x@ypos + x@spread, font = 2,
         x@id, srt=90, col = c("red", "forestgreen", "black")[x@ypos])
    abline(h = (1:N) +0.5, col = "gray")
    invisible(x)
})
