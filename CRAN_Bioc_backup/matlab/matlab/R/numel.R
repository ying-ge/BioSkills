###
### $Id: numel.R 29 2022-05-30 23:02:22Z proebuck $
###
### Provides number of elements.
###


##-----------------------------------------------------------------------------
numel <- function(A, varargin) {
    if (!missing(varargin)) {
        stop("not implemented")         # need example
    }

    prod(matlab::size(A))
}

