### -----------------------------------------------------------------
### Getters
###
setMethod("ID", "XMatrix", function(x) x@ID)
setMethod("name", "XMatrix", function(x) x@name)
setMethod("matrixClass", "XMatrix", function(x) x@matrixClass)
setMethod("Matrix", "XMatrix", function(x) x@profileMatrix)
setMethod("strand", "XMatrix", function(x) x@strand)
setMethod("bg", "XMatrix", function(x) x@bg)
setMethod("tags", "XMatrix", function(x) x@tags)
setMethod("matrixType", "PFMatrix", function(x) "PFM")
setMethod("matrixType", "ICMatrix", function(x) "ICM")
setMethod("matrixType", "PWMatrix", function(x) "PWM")
setMethod("pseudocounts", "PWMatrix", function(x) x@pseudocounts)
setMethod("schneider", "ICMatrix", function(x) x@schneider)
setMethod("colSums", "XMatrix", function(x) colSums(x@profileMatrix))
setMethod("rowSums", "XMatrix", function(x) rowSums(x@profileMatrix))
setMethod("dim", "XMatrix", function(x) dim(x@profileMatrix))

### -----------------------------------------------------------------
### Setters
###
setReplaceMethod("ID", "XMatrix",
                 function(x, value){
                   x@ID = value
                   return(x)
                 }
                 )

setReplaceMethod("name", "XMatrix",
                 function(x, value){
                   x@name = value
                   return(x)
                 }
                 )

setReplaceMethod("matrixClass", "XMatrix",
                 function(x, value){
                   x@matrixClass = value
                   return(x)
                 }
                 )

setReplaceMethod("strand", "XMatrix",
                 function(x, value){
                   x@strand = value
                   return(x)
                 }
                 )

setReplaceMethod("bg", "XMatrix",
                 function(x, value){
                   x@bg = value
                   return(x)
                 }
                 )

setReplaceMethod("Matrix", "XMatrix",
                 function(x, value){
                   x@profileMatrix= value
                   return(x)
                 }
                 )

setReplaceMethod("pseudocounts", "PWMatrix",
                 function(x, value){
                   x@pseudocounts = value
                   return(x)
                 }
                 )

setReplaceMethod("schneider", "ICMatrix",
                 function(x, value){
                   x@schneider = value
                   return(x)
                 }
                 )

### -----------------------------------------------------------------
### Updating and cloing
###
### An object is either 'update'd in place (usually with a replacement
### method) or 'clone'd (copied), with specified slots/fields overridden.
setMethod("update", "ANY",
          function(object, ..., check=TRUE){
            initialize(object, ...)
          }
          )
setMethod("clone", "ANY",
          function(x, ...){
            if(nargs() > 1L)
              initialize(x, ...)
            else
              x
          }
          )

### ---------------------------------------------------------
### Some utilities functions for XMatrix object
###
setMethod("length", "XMatrix",
# gets the pattern length in nucleotides (i.e. number of columns in the matrix)
          function(x){
            ncol(Matrix(x))
          }
          )

setMethod("ncol", "XMatrix",
          function(x){
            ncol(Matrix(x))
          }
          )

setMethod("reverseComplement", "XMatrix",
          function(x){
            ans = x
            Matrix(ans) = reverseComplement(Matrix(x))
            if(length(strand(x)) != 0L)
              strand(ans) = ifelse(strand(x) == "+", "-", "+")
            return(ans)
          }
          )


