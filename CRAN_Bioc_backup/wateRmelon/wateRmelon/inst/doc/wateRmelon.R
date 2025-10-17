## ----eval= FALSE, echo= TRUE--------------------------------------------------
#  install.packages('BiocManager')
#  BiocManager::install('wateRmelon')
#  library(wateRmelon)

## ----echo=FALSE, eval=TRUE, message=FALSE-------------------------------------
library(wateRmelon)
data(melon)

## ----eval= FALSE, echo= TRUE--------------------------------------------------
#  install.packages('devtools')
#  devtools::install_github('schalkwyk/wateRmelon')

## ----qs-----------------------------------------------------------------------
data(melon)
dim(melon)
# Quality filter using default thresholds
melon.pf <- pfilter(melon)

# Normalize using one of the many available methods
melon.dasen.pf <- dasen(melon.pf)

# Extract Betas for downstream analysis
norm_betas <- betas(melon.dasen.pf)

## ----read, eval=FALSE, echo=TRUE----------------------------------------------
#  mlumi <- readEPIC('path/to/idats')

## ----outlyx-------------------------------------------------------------------
outliers <- outlyx(melon, plot=TRUE)
print(outliers)
# remove outliers with melon[,!outliers$out]

## ----bscon--------------------------------------------------------------------
bsc <- bscon(melon)
hist(bsc, xlab = c(0, 100))

## ----pfilter------------------------------------------------------------------
melon.pf <- pfilter(melon)

## ----horv---------------------------------------------------------------------
agep(melon, method='all')

agep(melon, method='horvath')

## ----sex----------------------------------------------------------------------
estimateSex(betas(melon), do_plot=TRUE)

## ----cct, eval=FALSE----------------------------------------------------------
#  estimateCellCounts.wateRmelon(melon, referencePlatform = "IlluminaHumanMethylation450k")
#  estimateCellCounts.wateRmelon(melon, referencePlatform = "IlluminaHumanMethylationEPIC") # change reference

## ----norm---------------------------------------------------------------------
dasen.melon <- dasen(melon) # Use whichever method you would like to use. 

## ----qual---------------------------------------------------------------------
das <- dasen(melon)
qu <- qual(betas(melon), betas(das))
plot(qu[,1], qu[,2])

## ----dmrse--------------------------------------------------------------------
dmrse_row(melon.pf)
dmrse_row(melon.dasen.pf) # Slightly better standard errores

## ----genki--------------------------------------------------------------------
genki(melon.pf)
genki(melon.dasen.pf)

## ----XCI----------------------------------------------------------------------
seabi(melon.pf, sex=pData(melon.pf)$sex, X=fData(melon.pf)$CHR=='X')
seabi(melon.dasen.pf, sex=pData(melon.dasen.pf)$sex, X=fData(melon.dasen.pf)$CHR=='X')

## -----------------------------------------------------------------------------
bet <- betas(melon)
pwod_bet <- pwod(bet)

# Statistical Analysis using pwod_bet

## ----sessionInfo, echo=FALSE--------------------------------------------------
sessionInfo()

