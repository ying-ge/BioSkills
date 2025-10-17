###
### $Id: isempty.R 29 2022-05-30 23:02:22Z proebuck $
###
### Determine if object is empty.
###


##-----------------------------------------------------------------------------
isempty <- function(A) {
    any(matlab::size(A) == 0)
}

