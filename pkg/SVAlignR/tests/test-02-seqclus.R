library(SVAlignR)
data(longreads)
sequences <- longreads$connection
sequences <- sequences[!duplicated(sequences)]
alfa <- Cipher(sequences)
enc <- encode(alfa, sequences)

head(data.frame(sequences, enc))

# specify Levenshtein
sc1 <- SequenceCluster(enc, method = "levenshtein", NC = 12)
plot(sc1@hc)

# default Needelman-Wunsch
sc <- SequenceCluster(enc, NC = 20)
plot(sc)
# Other plot forms
plot(sc, type = "clipped")
plot(sc, type = "unrooted")


# make sure unknown methods fail.
try( sc <- SequenceCluster(enc, method = "clustal") )
try( sc3 <- SequenceCluster(enc, method = "leven") ) # Now OK

# same with unkown plot types
try( plot(sc, type = "doodle") )

#### update cluster number
sc <- updateClusters(sc, NC = 14)
plot(sc, type = "unrooted")
