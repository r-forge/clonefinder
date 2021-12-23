library(SVAlignR)
data(longreads)
sequences <- longreads$connection
sequences <- sequences[!duplicated(sequences)]
alfa <- Cipher(sequences)
enc <- encode(alfa, sequences)

head(data.frame(sequences, enc))

# specify Levenshtein
sc1 <- SequenceCluster(enc, method = "levenshtein")
plot(sc1@hc)

# default Needelman-Wunsch
sc <- SequenceCluster(enc)
plot(sc)
# change number of branches
plot(sc, NK = 20)
# Other plot forms
plot(sc, type = "clipped", NK=12)
plot(sc, type = "unrooted", NK=12)


# make sure unknown methods fail.
try( sc <- SequenceCluster(enc, method = "clustal") )
try( sc3 <- SequenceCluster(enc, method = "leven") )

# same with unkown lpot types
try( plot(sc, type = "doodle") )
