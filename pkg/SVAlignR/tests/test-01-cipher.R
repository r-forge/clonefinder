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

toolong <- paste(c(LETTERS, letters, 0:9,
                   "!","@","#", "$", "%", "^", "&", "*", "(", ")"),
                 collapse = "-")
try( beta <- Cipher(toolong) )

odd <- c(motif, "0-0-50-74-61")
try( encode(alfa, odd) )
try( encode(odd, alfa) )


