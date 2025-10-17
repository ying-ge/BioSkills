###
### $Id: ndims.R 29 2022-05-30 23:02:22Z proebuck $
###
### Provides the number of dimensions.
###


##-----------------------------------------------------------------------------
ndims <- function(A) {
    length(matlab::size(A))
}

