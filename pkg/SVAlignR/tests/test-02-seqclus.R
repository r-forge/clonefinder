library(SVAlignR)
data(longreads)
N <- 32
sequences <- longreads$connection[1:N]

# specify Levenshtein (should get two warnings)
sc1 <- SequenceCluster(sequences, method = "levenshtein", NC = 4)
plot(sc1@hc)

# after assigning names, should get one warning
names(sequences) <- rownames(longreads)[1:N]
sc1 <- SequenceCluster(sequences, method = "levenshtein", NC = 4)

# after removing duplicates, should get no warnings
sequences <- sequences[!duplicated(sequences)]
sc1 <- SequenceCluster(sequences, method = "levenshtein", NC = 4)

# default Needelman-Wunsch
sc <- SequenceCluster(sequences, NC = 4)
plot(sc)
# Other plot forms
plot(sc, type = "clipped")
plot(sc, type = "unrooted")


# make sure unknown methods fail.
try( sc <- SequenceCluster(sequences, method = "clustal") )
try( sc3 <- SequenceCluster(sequences, method = "leven") ) # Now OK

# same with unknown plot types
try( plot(sc, type = "doodle") )

#### update cluster number
sc <- updateClusters(sc, NC = 3)
plot(sc, type = "unrooted")
