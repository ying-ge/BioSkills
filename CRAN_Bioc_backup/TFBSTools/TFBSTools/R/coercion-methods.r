### -----------------------------------------------------------------
### XMatrix Coercion
###
setAs("XMatrix", "matrix", function(from) Matrix(from))
setAs("matrix", "PFMatrix", function(from) PFMatrix(profileMatrix=from))
setAs("matrix", "PWMatrix", function(from) PWMatrix(profileMatrix=from))
setAs("matrix", "ICMatrix", function(from) ICMatrix(profileMatrix=from))
setMethod("as.matrix", "XMatrix",
          function(x){Matrix(x)}
          )

### -----------------------------------------------------------------
### XMatrixList Coercion
setAs("XMatrixList", "data.frame", function(from){})

### -----------------------------------------------------------------
### SiteSet coercion
setAs("SiteSet", "data.frame", function(from){
      as(as(from, "DataFrame"), "data.frame")
          })

setAs("SiteSet", "DataFrame", function(from){
  if(length(from) == 0L)
    return(DataFrame())
  
  seqs <- DNAStringSet(views(from))
  seqs[strand(from) == "-"] <- 
    reverseComplement(seqs[strand(from) == "-"])
  absScore <- score(from) 
  relScore <- relScore(from)
  ans <- DataFrame(seqnames=Rle(from@seqname),
                   source=Rle(from@sitesource),
                   feature=Rle(from@sitesource),
                   start=start(views(from)),
                   end=end(views(from)),
                   absScore=absScore,
                   relScore=relScore,
                   strand=Rle(strand(from)),
                   ID=Rle(from@pattern@ID),
                   TF=Rle(from@pattern@name),
                   class=Rle(paste(from@pattern@matrixClass, collapse="; ")),
                   siteSeqs=seqs
                   )
  return(ans)
  }
)

setMethod("as.data.frame", "SiteSet",
          function(x){
            as(x, "data.frame")
          })

setAs("SiteSet", "GRanges", function(from){
      from.DataFrame <- as(from, "DataFrame")
      ans <- GRanges(seqnames=from.DataFrame[["seqnames"]],
                     ranges=IRanges(start=from.DataFrame[["start"]],
                                    end=from.DataFrame[["end"]]),
                     strand=from.DataFrame[["strand"]],
                     from.DataFrame[!colnames(from.DataFrame) %in% 
                                    c("seqnames", "start",
                                      "end", "strand")])
      return(ans)
          })

### -----------------------------------------------------------------
### SiteSetList coercion
setAs("SiteSetList", "data.frame", function(from){
      ans <- as(as(from, "DataFrame"), "data.frame")
      return(ans)
          })

setAs("SiteSetList", "DataFrame", function(from){
      ans <- do.call(rbind, lapply(from, as, "DataFrame"))
      return(ans)
          })

setMethod("as.data.frame", "SiteSetList",
          function(x){
            ans <- as(x, "data.frame")
            return(ans)
          })

setAs("SiteSetList", "GRanges", function(from){
      from.DataFrame <- as(from, "DataFrame")
      ans <- GRanges(seqnames=from.DataFrame[["seqnames"]],
                     ranges=IRanges(start=from.DataFrame[["start"]],
                                    end=from.DataFrame[["end"]]),
                     strand=from.DataFrame[["strand"]],
                     from.DataFrame[!colnames(from.DataFrame) %in%
                                    c("seqnames", "start",
                                      "end", "strand")])
      return(ans)
          })

### ----------------------------------------------------------------
### SitePairSet coercion
setAs("SitePairSet", "data.frame", function(from){
      as(as(from, "DataFrame"), "data.frame")
          })

setAs("SitePairSet", "DataFrame", function(from){
      ans <- cbind(as(from@siteset1, "DataFrame"),
                   as(from@siteset2, "DataFrame"))
      return(ans)
          })

setMethod("as.data.frame", "SitePairSet",
          function(x){
            ans <- as(x, "data.frame")
            return(ans)
          })

setAs("SitePairSet", "GRanges", function(from){
      from.DataFrame <- as(from, "DataFrame")
      columns <- ncol(from.DataFrame)
      from.DataFrame.target <- from.DataFrame[ ,1L:(columns/2L)]
      from.DataFrame.query <- setNames(
                                from.DataFrame[ ,(columns/2L+1L):columns],
                                names(from.DataFrame[ ,1L:(columns/2L)]))
      ansTarget <- GRanges(seqnames=from.DataFrame.target[["seqnames"]],
                           ranges=IRanges(start=
                                          from.DataFrame.target[["start"]],
                                          end=
                                          from.DataFrame.target[["end"]]),
                           strand=from.DataFrame.target[["strand"]],
                           from.DataFrame.target[
                                 !colnames(from.DataFrame.target) %in%
                                 c("seqnames", "start",
                                   "end", "strand")]
                           )
                           
      ansQuery <- GRanges(seqnames=from.DataFrame.query[["seqnames"]],
                          ranges=IRanges(start=
                                         from.DataFrame.query[["start"]],
                                         end=
                                         from.DataFrame.query[["end"]]),
                          strand=from.DataFrame.query[["strand"]],
                          from.DataFrame.query[
                                 !colnames(from.DataFrame.query) %in%
                                 c("seqnames", "start", "end", 
                                   "strand")]
                          )
      ans <- GRangesList(targetTFBS=ansTarget, queryTFBS=ansQuery)
      return(ans)
          })

### -----------------------------------------------------------------
### SitePairSetList coersion
setAs("SitePairSetList", "data.frame", function(from){
      as(as(from, "DataFrame"), "data.frame")
          })

setAs("SitePairSetList", "DataFrame", function(from){
      ans <- do.call(rbind, lapply(from, as, "DataFrame"))
      return(ans)
          })

setMethod("as.data.frame", "SitePairSetList",
          function(x){
            ans <- as(x, "data.frame")
            return(ans)
          })

setAs("SitePairSetList", "GRanges", function(from){
      from.DataFrame <- as(from, "DataFrame")
      columns <- ncol(from.DataFrame)
      from.DataFrame.target <- from.DataFrame[ ,1L:(columns/2L)]
      from.DataFrame.query <- setNames(
                                from.DataFrame[ ,(columns/2L+1L):columns],
                                names(from.DataFrame[ ,1L:(columns/2L)]))
      ansTarget <- GRanges(seqnames=from.DataFrame.target[["seqnames"]],
                           ranges=IRanges(start=
                                          from.DataFrame.target[["start"]],
                                          end=
                                          from.DataFrame.target[["end"]]),
                           strand=from.DataFrame.target[["strand"]],
                           from.DataFrame.target[
                                 !colnames(from.DataFrame.target) %in%
                                 c("seqnames", "start",
                                   "end", "strand")]
                           )

      ansQuery <- GRanges(seqnames=from.DataFrame.query[["seqnames"]],
                          ranges=IRanges(start=
                                         from.DataFrame.query[["start"]],
                                         end=
                                         from.DataFrame.query[["end"]]),
                          strand=from.DataFrame.query[["strand"]],
                          from.DataFrame.query[
                                 !colnames(from.DataFrame.query) %in%
                                 c("seqnames", "start", "end",
                                   "strand")]
                          )
      ans <- GRangesList(targetTFBS=ansTarget, queryTFBS=ansQuery)
      return(ans)
          })
