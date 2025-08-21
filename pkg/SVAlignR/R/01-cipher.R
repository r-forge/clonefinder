### Working With Ciphers

setClass("Cipher",
         slots = c(
           forward = "character",
           reverse = "character",
           bytes = "numeric")
         )

Cipher <- function(sampleText, split = "-", extras = c("-" = ":", "?" = "?")) {
  nb <- 1
  shatter <- strsplit(sampleText, split) # decompose into lists of symbols
  U <- sort(unique(unlist(shatter))) # find the unique symbols
  base <- c(LETTERS[-c(10,24)], letters, 1:9, 0,
            "!", "@", "#", "$", "%", "^", "&", "|", "_", ";", "/", ",")
  if (length(U) > length(base)) {
    warning("Input has more than ", length(base), " letters. ",
            "Using two-byte codes.\n")
    nb <- 2
    multi <- 1 + trunc(length(U)/length(base))
    base <-paste0(rep(LETTERS[1:multi], each = length(base)),
                    rep(base, times = multi))
  }
  forward  <- base[1:length(U)]
  names(forward) <- U
  if (is.null(names(extras)) | any(names(extras) == "")) {
    stop("All elements of 'extras' must have names.\n")
  }
  if (any(duplicated(c(names(forward), extras)))) {
    stop("Values of 'extras' cannot match letter in the 'forward' cipher.\n")
  }
  reverse <- names(forward)
  names(reverse) <- forward
  reverse <- c(reverse, extras)
  if (any(duplicated(names(reverse)))) {
    stop("Names of 'extras' cannot include symbols in 'sampleText'.\n")
  }
  new("Cipher", forward = forward, reverse = reverse, bytes = nb)
}

.xlate <- function(txt, language, split, collapse, nb = 1) {
  txt <- strsplit(txt, split)
  if (nb > 1) {
    txt <- lapply(txt, function(TX) {
      tmat <- matrix(TX, ncol = nb, byrow = TRUE)
      apply(tmat, 1, paste, collapse = "")
    })
  }
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
  .xlate(text, cipher@forward, "\\-", "", 1)
}

decode <- function(cipher, text) {
  .xlate(text, cipher@reverse, "", "-", cipher@bytes)
}

