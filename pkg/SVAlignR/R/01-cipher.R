### Working With Ciphers

setClass("Cipher",
         slots = c(
           forward = "character",
           reverse = "character")
         )

Cipher <- function(sampleText, split = "-", extras = c("-" = ":", "?" = "?")) {
  shatter <- strsplit(sampleText, split) # decompose into lists of symbols
  U <- sort(unique(unlist(shatter))) # find the unique symbols
  base <- c(LETTERS, letters, 1:9, 0, "!", "@", "#", "$", "%", "^", "&", "|")
  if (length(U) > length(base)) {
    stop("We can only handle codes with at most ", length(base), " letters.\n")
  }
  forward  <- base[1:length(U)]
  names(forward) <- U
  reverse <- names(forward)
  names(reverse) <- forward
  reverse <- c(reverse, extras)
  new("Cipher", forward = forward, reverse = reverse)
}

.xlate <- function(txt, language, split, collapse) {
  txt <- strsplit(txt, split)
  temp <- unique(unlist(txt))
  fail <- !(temp %in% names(language))
  if (any(fail)) {
    stop("Text includes unknown characters: ", temp[fail])
  }
  sapply(txt, function(X) {
    paste(language[X], collapse = collapse)
  })
}

encode <- function(cipher, text) {
  .xlate(text, cipher@forward, "\\-", "")
}

decode <- function(cipher, text) {
  .xlate(text, cipher@reverse, "", "-")
}

