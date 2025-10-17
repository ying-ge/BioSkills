###
### $Id: mkconstarray.R 29 2022-05-30 23:02:22Z proebuck $
###
### Create a constant array of specified class.
###


##-----------------------------------------------------------------------------
mkconstarray <- function(class.type = c("character",
                                        "complex",
                                        "double",
                                        "integer",
                                        "logical",
                                        "numeric"),
                         value,
                         size) {
     matlab::repmat(as(value, match.arg(class.type)), size)
}

