## ----style, echo = FALSE, results = 'asis', message = FALSE-------------------
BiocStyle::markdown()
library(knitr)

## ----loadingPackage, warning=FALSE, message=FALSE-----------------------------
library(nucleoSim)

## ----demoMap, warning=FALSE, message=FALSE, collapse=TRUE---------------------
wp.num           <- 20         ### Number of well-positioned nucleosomes
wp.del           <- 5          ### Number of well-positioned nucleosomes to delete
wp.var           <- 30         ### variance associated with the starting 
                               ###   position of the sequences of the 
                               ###   well-positioned nucleosomes
fuz.num          <- 5          ### Number of fuzzy nucleosomes
fuz.var          <- 50         ### Variance associated with the starting 
                               ###   positions of the sequences for the 
                               ###   fuzzy nucleosomes
max.cover        <- 70         ### Maximum sequences associated with one 
                               ###   nucleosome (default: 100)
nuc.len          <- 147        ### Length of the nucleosome (default: 147)
len.var          <- 12         ### variance associated with the length of 
                               ###   the sequences (default: 10)
lin.len          <- 20         ### Length of the DNA linker (default: 20)
distr            <- "Normal"   ### Type of distribution to use

rnd.seed         <- 210001     ### Set seed when result needs to be reproducible
        
#### Create a synthetic nucleosome map
nucleosomeMap <- syntheticNucMapFromDist(wp.num=wp.num, wp.del=wp.del, 
                        wp.var=wp.var, fuz.num=fuz.num, fuz.var=fuz.var, 
                        max.cover=max.cover, nuc.len=nuc.len, len.var=len.var, 
                        lin.len=lin.len, rnd.seed=rnd.seed, distr=distr)

#### The start positions of all well-positioned nucleosomes
nucleosomeMap$wp.starts

#### The number of sequences associated with each well-positioned nucleosome
nucleosomeMap$wp.nreads

#### IRanges object containing all sequences for the well-positioned nucleosomes
head(nucleosomeMap$wp.reads, n = 2)

#### The start positions of all fuzzy nucleosomes
nucleosomeMap$fuz.starts

#### The number of sequences associated with each fuzzy nucleosome
nucleosomeMap$fuz.nreads

#### A IRanges object containing all sequences for the fuzzy nucleosomes
head(nucleosomeMap$fuz.reads, n = 2)

#### A IRanges object containing all sequences for all nucleosomes
head(nucleosomeMap$syn.reads, n = 2)

## ----showMap, fig.align='center', fig.height=5, fig.width=8-------------------
#### Create visual representation of the synthetic nucleosome map
plot(nucleosomeMap, xlab="Position", ylab="Coverage")

## ----demoMapTiling, warning=FALSE, message=FALSE, collapse=TRUE---------------
as.ratio         <- TRUE       ### Activate the simulation of hybridization data
rnd.seed         <- 212309     ### Set seed when result needs to be reproducible
        
#### Create a synthetic nucleosome map with hybridization data
nucleosomeMapTiling <- syntheticNucMapFromDist(wp.num=10, wp.del=2, wp.var=20, 
                                    fuz.num=1, fuz.var=32, max.cover=50,
                                    nuc.len=145, len.var=3, lin.len=40, 
                                    rnd.seed=rnd.seed, as.ratio=as.ratio,
                                    distr="Uniform")

#### Control sequences for hybridization data (only when as.ratio = TRUE)
head(nucleosomeMapTiling$ctr.reads, n=4)

#### Ratio for hybridization data (only when as.ratio = TRUE)
head(nucleosomeMapTiling$syn.ratio, n=4)

#### Create visual representation of the synthetic nucleosome map
plot(nucleosomeMapTiling)

## ----demoSample, warning=FALSE, message=FALSE, collapse=TRUE------------------
wp.num           <- 30            ### Number of well-positioned nucleosomes
wp.del           <- 10            ### Number of well-positioned nucleosomes 
                                  ###   to delete
wp.var           <- 30            ### variance associated with the starting 
                                  ###   positions of the sequences for the 
                                  ###   well-positioned nucleosomes
fuz.num          <- 10            ### Number of fuzzy nucleosomes
fuz.var          <- 50            ### Variance associated with the starting 
                                  ###   positions of the sequences for the 
                                  ###   fuzzy nucleosomes
max.cover        <- 90            ### Maximum paired-end reads associated with 
                                  ###   one nucleosome (default: 100)
nuc.len          <- 147           ### Length of the nucleosome (default: 147)
len.var          <- 12            ### variance associated with the distance 
                                  ###   between start positions of 
                                  ###   paired-end reads (default: 10)
lin.len          <- 20            ### Length of the DNA linker (default: 20)
read.len         <- 45            ### Length of the generated forward and 
                                  ###   reverse reads (default: 40)
distr            <- "Uniform"     ### Type of distribution to use
offset           <- 10000         ### The number of bases used to offset 
                                  ###   all nucleosomes and reads

rnd.seed         <- 202           ### Set seed when result needs to be 
                                  ###   reproducible
        
#### Create Uniform sample
nucleosomeSample <- syntheticNucReadsFromDist(wp.num=wp.num, wp.del=wp.del, 
                        wp.var=wp.var, fuz.num=fuz.num, fuz.var=fuz.var, 
                        max.cover=max.cover, nuc.len=nuc.len, len.var=len.var, 
                        read.len=read.len, lin.len=lin.len, rnd.seed=rnd.seed, 
                        distr=distr, offset=offset)

#### The central position of all well-positioned nucleosomes with the
#### number of paired-end reads each associated with each one
head(nucleosomeSample$wp, n = 2)

#### The central position of all fuzzy nucleosomes with the
#### number of paired-end reads each associated with each one
head(nucleosomeSample$fuz, n = 2)

#### A data.frame with the name of the synthetic chromosome, the starting 
#### position, the ending position and the direction of all forward and
#### reverse reads
head(nucleosomeSample$dataIP, n = 2)

## ----showSample, fig.align='center', fig.height=5, fig.width=8----------------
#### Create visual representation of the synthetic nucleosome sample
plot(nucleosomeSample, xlab="Position", ylab="Coverage (number of reads)")

## ----demoSampleFromMap, warning=FALSE, message=FALSE, collapse=TRUE-----------

#### A nucleosome map has already been created 
class(nucleosomeMap)

#### 
read.len    <- 45   ### The length of the reverse and forward reads
offset      <- 500  ### The number of bases used to offset all nucleosomes and reads

#### Create nucleosome sample
nucleosomeSampleFromMap <- syntheticNucReadsFromMap(nucleosomeMap, 
                                        read.len=read.len, offset=offset)

#### A data.frame with the name of the synthetic chromosome, the starting 
#### position, the ending position and the direction of all forward and
#### reverse reads
head(nucleosomeSampleFromMap$dataIP, n = 2)

## ----sessionInfo, echo=FALSE--------------------------------------------------
sessionInfo()

