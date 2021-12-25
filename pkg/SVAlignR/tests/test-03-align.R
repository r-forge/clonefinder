library(SVAlignR)
data(longreads)
seqs <- longreads$connection[1:15]
pad <- c(rep("0", 9), rep("", 6))
names(seqs) <- paste("LR", pad, 1:length(seqs), sep = "")
seqs <- seqs[!duplicated(seqs)]
alfa <- Cipher(seqs)
### Since the encoding here is identical to the re-encoding internally to 'align',
### we manually create an alternate one for debugging purposes
ff <- rev(LETTERS)[1:6]
ff <- LETTERS[c(13:14, 15:18)]
rr <- names(ff) <- names(alfa@forward)
names(rr) <- ff
rr <- c(rr, "-" = ":", "?" = "?")
beta <- new("Cipher", forward = ff, reverse = rr)

mysub <- SVAlignR:::makeMYSUB(match = 2, mismatch = -6)
### Use beta-encoding
enc <- encode(beta, seqs)
svb <- align(enc, mysub = mysub)
svb
db <- decode(beta, svb$alignedOriginal)
ob <- db[order(names(db))]
decode(beta, svb$cons)
### Use alfa-encoding
reup <- encode(alfa, seqs)
sva <- align(reup, mysub = mysub)
da <- decode(alfa, sva$alignedOriginal)
oa <- da[order(names(da))]
decode(alfa, sva$cons)
### compare
data.frame(oa, ob)
all(oa == ob)

## test the simple shortcut method
ab <- alignBranch(seqs, mysub)
showme(ab)
