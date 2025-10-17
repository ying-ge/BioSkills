### ---------------------------------------------------------------
### The real wrapper function for MEME
###
run_MEME <- function(inputFastaFn, binary="meme", seqtype="DNA", 
                    arguments=list()){
  valueArguments <- c("-mod", "-nmotifs", "-evt", "-nsites", 
                        "-minsites", "-maxsites", "-wnsites",
                        "-w", "-minw", "-maxw", "-wg", "-ws",
                        "-bfile", "-maxiter", "-distance",
                        "-psp", "-prior", "-b", "-plib",
                        "-spfuzz", "-spmap", "-cons", "-heapsize",
                        "-maxsize", "-p", "-time", "-sf")
  booleanArguments <- c("-nomatrim", "-noendgaps", "-revcomp", "-pal",
                       "-x_branch", "-w_branch")
  arguments1 <- arguments[names(arguments) %in% valueArguments]
  arguments1 <- paste(names(arguments1), format(arguments1, scientific=FALSE),
                      collapse=" ")
  arguments2 <- arguments[names(arguments) %in% booleanArguments]
  arguments2 <- paste(names(arguments2), collapse=" ")
  arguments <- paste(arguments1, arguments2)
  options <- c(inputFastaFn, "-text", 
               ifelse(seqtype=="DNA", "-dna", "-protein"), 
               arguments)
  message(paste(binary, paste(options, collapse=" ")))
  memeOutput <- system2(command=binary, args=options, stderr=FALSE, stdout=TRUE)
  if(!is.null(attributes(memeOutput))){
    if(attributes(memeOutput)$status == 1L){
      testFn <- paste0("MEME-", basename(inputFastaFn))
      file.copy(inputFastaFn, testFn)
      stop("MEME run with error! Please test MEME with input file ", testFn)
    }else if(attributes(memeOutput)$status == 127L){
      stop("Can't run MEME!")
    }
  }
  ans <- try(parseMEMEOutput(memeOutput), silent=TRUE)
  if(class(ans) == "try-error"){
    ans <- parseMEMEOutput412(memeOutput)
  }
  return(ans)
}

### -------------------------------------------------------------
### The MEME method
### Exported!
setMethod("runMEME", "character",
          function(x, binary="meme", seqtype="DNA", arguments=list(), 
                   tmpdir=tempdir()){
            seqtype = match.arg(seqtype, c("DNA", "AA"))
            ans = run_MEME(x, binary=binary, seqtype=seqtype, 
                           arguments=arguments)
            subjectSeqs = switch(seqtype,
                                 "DNA"=readDNAStringSet(filepath=x, 
                                                        format="fasta"),
                                 "AA"=readAAStringSet(filepath=x, 
                                                      format="fasta")
                                 )
            ans = MotifSet(motifList=ans[["motifList"]], 
                           motifEvalues=ans[["motifEvalues"]], 
                           subjectSeqs=subjectSeqs)
            return(ans)
          }
          )

setMethod("runMEME", "DNAStringSet",
          function(x, binary="meme", seqtype="DNA", arguments=list(), 
                   tmpdir=tempdir()){
            tmpFile = tempfile(pattern="MEME_", tmpdir=tmpdir, 
                               fileext = ".fasta")
            if(is.null(names(x)))
              names(x) <- seq_along(x)
            writeXStringSet(x, filepath=tmpFile, format="fasta")
            on.exit(unlink(tmpFile))
            ans = run_MEME(tmpFile, binary=binary, 
                           seqtype=seqtype, arguments=arguments)
            ans = MotifSet(motifList=ans[["motifList"]], 
                           motifEvalues=ans[["motifEvalues"]], subjectSeqs=x)
            return(ans)
          }
          )

### -----------------------------------------------------------------
### Parse the output from MEME
### Exported!
parseMEMEOutput <- function(x){
  if(length(x) == 1L){
    ## The path of file
    memeOutput <- readLines(x)
  }else{
    memeOutput <- x
  }
  
  # get the version of meme used.
  oneLine <- memeOutput[grep("^MEME version", memeOutput)]
  version = strsplit(oneLine, " ")[[1]][3]
  
  # get the command
  oneLine = memeOutput[grep("^command:", memeOutput)]
  command = gsub("^command: ", "", oneLine)
  
  
  # get the motifs information
  revcomp = grepl("revcomp", command)
  oneLines = memeOutput[grep("^MOTIF +\\d+", memeOutput)]
  splittedLines = strsplit(oneLines, "[[:blank:]]+")
  motifWidths = as.integer(sapply(splittedLines, "[", 6))
  motifOccurrences = as.integer(sapply(splittedLines, "[", 9))
  motifEvalues = as.numeric(sapply(splittedLines, "[", 15))
  
  # get the motifs names
  indexNames = grep("sorted by position p-value", memeOutput)
  oneLines = memeOutput[indexNames]
  motifNames = lapply(strsplit(oneLines, "[[:blank:]]+"), "[", c(2,3))
  motifNames = sapply(motifNames, paste0, collapse=" ")
  
  # get the motifs ranges
  motifList = list()
  for(i in seq_len(length(indexNames))){
    oneLines = memeOutput[seq(from=indexNames[i]+4, 
                              to=indexNames[i]+4+motifOccurrences[i]-1)]
    splittedLines = strsplit(oneLines, "[[:blank:]]+")
    strands <- NULL
    starts <- NULL
    if(revcomp){
      strands <- sapply(splittedLines, "[", 2)
      starts <- as.integer(sapply(splittedLines, "[", 3))
    }else{
      strands <- Rle("+", length(splittedLines))
      starts <- as.integer(sapply(splittedLines, "[", 2))
    }
    oneRange <- GRanges(seqnames=sapply(splittedLines, "[", 1), 
                        ranges=IRanges(start=starts,
                                       width=motifWidths[i]),
                        strand=strands,
                        score=as.numeric(sapply(splittedLines, 
                                                "[", ifelse(revcomp, 4, 3)))
    )
    motifList = c(motifList, oneRange)
  }
  motifList = GRangesList(motifList)
  names(motifList) = motifNames
  ans = list(motifList=motifList, motifEvalues=motifEvalues)
  return(ans)
}

parseMEMEOutput412 <- function(x){
  if(length(x) == 1L){
    ## The path of file
    memeOutput <- readLines(x)
  }else{
    memeOutput <- x
  }
  
  # get the version of meme used.
  oneLine <- memeOutput[grep("^MEME version", memeOutput)]
  version = strsplit(oneLine, " ")[[1]][3]
  
  # get the command
  oneLine = memeOutput[grep("^command:", memeOutput)]
  command = gsub("^command: ", "", oneLine)
  
  
  # get the motifs information
  revcomp = grepl("revcomp", command)
  oneLines = memeOutput[grep("^MOTIF.*MEME-\\d+", memeOutput)]
  splittedLines = strsplit(oneLines, "[[:blank:]]+")
  motifWidths = as.integer(sapply(splittedLines, "[", 6))
  motifOccurrences = as.integer(sapply(splittedLines, "[", 9))
  motifEvalues = as.numeric(sapply(splittedLines, "[", 15))
  
  # get the motifs names
  indexNames = grep("sorted by position p-value", memeOutput)
  oneLines = memeOutput[indexNames]
  motifNames = lapply(strsplit(oneLines, "[[:blank:]]+"), "[", c(2,3))
  motifNames = sapply(motifNames, paste0, collapse=" ")
  
  # get the motifs ranges
  motifList = list()
  for(i in seq_len(length(indexNames))){
    oneLines = memeOutput[seq(from=indexNames[i]+4, 
                              to=indexNames[i]+4+motifOccurrences[i]-1)]
    splittedLines = strsplit(oneLines, "[[:blank:]]+")
    strands <- NULL
    starts <- NULL
    if(revcomp){
      strands <- sapply(splittedLines, "[", 2)
      starts <- as.integer(sapply(splittedLines, "[", 3))
    }else{
      strands <- Rle("+", length(splittedLines))
      starts <- as.integer(sapply(splittedLines, "[", 2))
    }
    oneRange <- GRanges(seqnames=sapply(splittedLines, "[", 1), 
                        ranges=IRanges(start=starts,
                                       width=motifWidths[i]),
                        strand=strands,
                        score=as.numeric(sapply(splittedLines, 
                                                "[", ifelse(revcomp, 4, 3)))
    )
    motifList = c(motifList, oneRange)
  }
  motifList = GRangesList(motifList)
  names(motifList) = motifNames
  ans = list(motifList=motifList, motifEvalues=motifEvalues)
  return(ans)
}
