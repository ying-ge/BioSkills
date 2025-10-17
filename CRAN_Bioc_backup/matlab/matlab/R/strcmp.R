###
### $Id: strcmp.R 29 2022-05-30 23:02:22Z proebuck $
###
### Compare strings.
###


##-----------------------------------------------------------------------------
strcmp <- function(S, T) {
    if (!is.character(S)) {
        stop(sprintf("argument %s must be character", sQuote("S")))
    }
    if (!is.character(T)) {
        stop(sprintf("argument %s must be character", sQuote("T")))
    }

    if (length(S) == length(T)) {
        all(S == T)
    } else {
        FALSE
    }
}

