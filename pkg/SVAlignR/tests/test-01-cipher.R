library(SVAlignR)

motif <- "0-50-74-0-50-74-25-26-35"
alfa <- Cipher(motif)
alfa

en <-encode(alfa, motif)
en

de <- decode(alfa, en)
de
de == motif

## bad Cipher inputs
try( Cipher(13) )
try( Cipher(motif, extras = 13) )
try( Cipher(motif, extras = c(z = "zero", 13)) )
try( Cipher(motif, extras = c(A = "zero")) )
try( Cipher(motif, extras = c("x" = 26)) )

odd <- c(motif, "0-0-50-74-61")
try( encode(alfa, odd) )
try( encode(odd, alfa) )

## Longer alphabets
nchar(IT <- intToUtf8(L <- c(38:44, 46:57, 59:62, 64:126)))
## omit hyphen, colon, and question
JT <- paste(strsplit(IT, "|")[[1]], collapse = "-")
gamma <- Cipher(JT)

b2 <- encode(gamma, JT)
decode(gamma, b2) == JT
