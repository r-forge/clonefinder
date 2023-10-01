library(SVAlignR)
data(longreads)
seqs <- longreads$connection[1:15]
pad <- c(rep("0", 9), rep("", 6))
names(seqs) <- paste("LR", pad, 1:length(seqs), sep = "")
seqs <- seqs[!duplicated(seqs)]

alfa <- Cipher(seqs)
seqclust <- SequenceCluster(seqs, NC = 5)
plot(seqclust)
coded <- encode(alfa, seqclust@rawSequences)
kmers <- lapply(1:11, function(N) {countWords(coded, N)})
#motifNodes <- lapply(kmers, function(x) x[names(x) %in% names(motifs)])
#save(alfa, motifs, kmers, motifNodes, decomp,
#     file = file.path(paths$scratch, "motifNodes.Rda"))
