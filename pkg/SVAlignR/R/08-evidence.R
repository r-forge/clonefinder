Vexpand <- function(cycles, seqs, alfa) {
  countVs <- function(cycs) {
    ss <- strsplit(cycs, 'v')
    sapply(ss, length) - 1
  }
  if( any(countVs(cycles) > 1) ) {
    stop("Input cycle is not a simple path!")
  }
  repVs <- function(seqs, V) {
    K <- 0
    count <- 1
    while(count > 0) {
      K <- K + 1
      vrep <- paste(rep(V, K), collapse = "-")
      count <- sum(grepl(vrep, seqs))
    }
    K - 1
  }
  # find repeated breakpoints
  enc <- encode(alfa, seqs)
  matches <- gregexpr("(.)\\1", enc)
  repeated_chars <- unique(unlist(regmatches(enc, matches)))
  stutter <- decode(alfa, substring(repeated_chars, 1, 1))
  for (V in stutter) {
    N <- repVs(seqs, V)
    if (N > 1) {
      hasV <- seqs[grep(V, cycles)]
      sap <- sapply(2:N, function(J, cycs) {
        vrep <- paste(rep(V, J), collapse = "-")
        sub(V, vrep, cycs)
      }, cycs = cycles)
      cycles <- c(cycles, as.vector(sap))
    }
  }
  longcM1 <- cycles[!duplicated(cycles)]
  longcM1
}

bothWays <- function(xc, seqs) {
  revert <- function(cycle) {
    X <- strsplit(cycle, "-")[[1]]
    paste(rev(X), collapse = "-")
  }
  doubles <- paste0(xc, "-", xc)
  support <- sapply(seqs, grep, x = doubles)
  RE <- sapply(seqs, revert)
  troppus <- sapply(RE, grep, x = doubles)
  comb <- lapply(1:length(support), function(J) {
    unique(c(support[[J]], troppus[[J]]))
  })
  names(comb) <- names(support)
  list(comb=comb, doubles = doubles)
}

getSupport <- function(comb, clues) {
  L <- sapply(comb, length)
  SP <- data.frame(seq=NA, cyc = NA)[-1,]
  for (J in 1:length(comb)) {
    if (L[J] > 0) {
      tp <- data.frame(seq = names(comb)[J], cyc = comb[J])
      colnames(tp) <- c("seq", "cyc")
      SP <- rbind(SP, tp)
    }
  }
  SP$weight <- clues@weights[SP$seq]
  SP
}

getCoverage <- function(xc, rawSP, seqs, doubles, alfa) {
  cvg <- lapply(1:length(doubles), function(K) {
    fd <- seqs[rawSP$seq[rawSP$cyc == K]]
    gp <- encode(alfa, fd)
    DD <- encode(alfa, doubles[[K]])
    M <- rep(0, length(strsplit(xc[K], "-")[[1]]))
    strt <- sapply(gp, function(S) str_locate(DD, S))
    if (is.null(dim(strt))) return(list(N = 0, RDS = "", TOT = 0, COV = M, WCOV = M))
    M2 <- c(M, M)
    W2 <- M2
    rds <- NULL
    for (J in 1:ncol(strt)) {
      if (!is.na(strt[1,J])) {
        A <- strt[1,J]
        B <- strt[2,J]
        M2[A:B] <- M2[A:B] + 1
        W2[A:B] <- W2[A:B] + rawSP$weight[J]
        rds <- c(rds, colnames(strt)[J])
      }
    }
    half <- length(M)
    M <- M2[1:half] + M2[half + (1:half)]
    W <- W2[1:half] + W2[half + (1:half)]
    list(N = ncol(strt), RDS = paste(rds, collapse = ";"),
         TOT = min(M), COV=M, WCOV = W)
  })
  cvg
}

showCoverage <- function(ID, evidence, colsams) {
  CC <- evidence$coverage[[ID]]
  SQ <- evidence$vspan[ID]
  N <- length(sq <- strsplit(SQ, "-")[[1]])
  lrnames <- strsplit(CC$RDS, ";")[[1]]
  lrnames <- lrnames[order(whereSeen[lrnames])]
  L <- length(lrnames)
  if (missing(colsams)) {
    colsams <- rep("grey", L)
    names(colsams) <- lrnames
  }
  target <- encode(alfa, SQ)
  target <- paste0(target, target)
  plot(c(1/2,N+1/2), c(1,(1+L)),  type = "n", xaxt = "n", yaxt = "n",
       xlab = "breakpoints", ylab = "long reads",
       main = paste("Putative Cycle", ID))
  mtext(lrnames, side = 2, at = 1:L, las = 2, cex = 0.5)
  text(1:N, L+1, sq)
  for (J in 1:L) {
    lr <- lrnames[J]
    clr <- colsams[lr]
    S <- encode(alfa, sqn[lr])
    finder <- str_locate(target, S)
    strt <- finder[1, 1]
    end <- finder[1, 2]
    if (end <= N) {
      lines(c(strt-0.5, end+0.5),c(J, J),  col = clr, lwd = 5 - round(log10(L)))
    } else {
      lines(c(strt-0.5, N+0.5),c(J, J),  col = clr, lwd = 5)
      lines(c(0.5, end-N+0.5),c(J, J),  col = clr, lwd = 5)
    }
    lines(c(strt-0.5, strt-0.5), c(J-0.4, J+0.4), col = "black", lwd=4)
  }
}

getEvidence <- function(cycM1, noxsqn, alfa, totalOnly = TRUE) {
  xcM1 <- Vexpand(cycM1, noxsqn, alfa)
  BW <- bothWays(xcM1, noxsqn)
  SP <- getSupport(BW$comb, seqclust)
  coverage <- getCoverage(xcM1, SP, noxsqn, BW$doubles, alfa)
  wmin <- sapply(coverage, function(X) min(X$WCOV))
  drat <- data.frame(N = unlist(lapply(coverage, function(X) X["N"])),
                     RDS = unlist(lapply(coverage, function(X) X["RDS"])),
                     TOT = unlist(lapply(coverage, function(X) X["TOT"])),
                     WTOT = wmin,
                     cycle = xcM1)
  if (totalOnly) {
    drat <- drat[drat$TOT > 0,]
  }
  list(coverage = coverage, circles = drat, vspan = xcM1)
}
