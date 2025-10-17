### =========================================================================
### The as.data.frame() generic
### -------------------------------------------------------------------------
###
### Note that base::as.data.frame is an S3 generic.
###
### Need to explicitly define this generic otherwise the implicit generic in
### package "base" would dispatch on all its arguments. Here we set dispatch
### on the 1st arg (the 'x' arg) only!

setGeneric("as.data.frame", signature="x")

