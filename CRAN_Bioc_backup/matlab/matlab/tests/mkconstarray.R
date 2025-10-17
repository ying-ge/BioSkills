###
### $Id: mkconstarray.R 22 2022-05-30 18:03:47Z proebuck $
###


##-----------------------------------------------------------------------------
test.mkconstarray <- function(input, expected) {
    output <- do.call(getFromNamespace("mkconstarray", "matlab"), input)
    identical(output, expected)
}

test.mkconstarray(list(class.type = "double", value = pi, size = 4), rep(pi, 4))

