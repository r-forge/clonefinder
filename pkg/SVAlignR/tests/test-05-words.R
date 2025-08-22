library(SVAlignR)

## single byte
data(longreads)
raw <- longreads$connection
alfa <- Cipher(raw)
coded <- encode(alfa, raw)
makeWords(coded, 3)
countWords(coded, 3)
countWords(coded, 3, alfa)

IT <- intToUtf8(L <- c(38:44, 46:57, 59:62, 64:126)) ## omit hyphen, colon, and question
JT <- paste(strsplit(IT, "|")[[1]], collapse = "-")
beta <- Cipher(JT)

twob <- encode(beta, JT)
double <- c(twob, twob)
countWords(double, 3)       # don't do this
countWords(double, 3, beta) # this is right

makeWords(double, 3)    # don't do this
makeWords(double, 3, 2) # this is right

