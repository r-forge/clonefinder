library(SVAlignR)
data(longreads)
seqs <- longreads$connection[1:15]
pad <- c(rep("0", 9), rep("", 6))
names(seqs) <- paste("LR", pad, 1:length(seqs), sep = "")
seqs <- seqs[!duplicated(seqs)]

try( makeSubsMatrix(match = "13") )
try( makeSubsMatrix(mismatch = "13") )

mysub <- makeSubsMatrix(match = 2, mismatch = -6)

## test the 'alignBranch' method
try( alignBranch(13) )

ab <- alignBranch(seqs, mysub)
showme(ab)
