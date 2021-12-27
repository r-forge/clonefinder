library(SVAlignR)
data(longreads)
sequences <- longreads$connection

# specify Levenshtein (should get two warnings)
sc1 <- SequenceCluster(sequences, method = "levenshtein", NC = 12)
plot(sc1@hc)

# after assigning names, should get one warning
names(sequences) <- rownames(longreads)
sc1 <- SequenceCluster(sequences, method = "levenshtein", NC = 12)

# after removing dupklcairtes, should get no warnings
sequences <- sequences[!duplicated(sequences)]
sc1 <- SequenceCluster(sequences, method = "levenshtein", NC = 12)

# default Needelman-Wunsch
sc <- SequenceCluster(sequences, NC = 20)
plot(sc)
# Other plot forms
plot(sc, type = "clipped")
plot(sc, type = "unrooted")


# make sure unknown methods fail.
try( sc <- SequenceCluster(sedquences, method = "clustal") )
try( sc3 <- SequenceCluster(enc, method = "leven") ) # Now OK

# same with unknown plot types
try( plot(sc, type = "doodle") )

#### update cluster number
sc <- updateClusters(sc, NC = 14)
plot(sc, type = "unrooted")
