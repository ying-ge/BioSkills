## ----init, echo=FALSE, results='hide'-----------------------------------------
## crop=NULL or FALSE =>  fix vignette rendering based on yihui/knitr#1796
## added also error=FALSE to include_graphics
knitr::opts_chunk$set(
  collapse = TRUE,
  tidy = FALSE,
  ## comment = "#>",
  error = FALSE,
  warning = FALSE,
  message = FALSE,
  crop = NULL                           
)

## check the output type
out_type <- knitr::opts_knit$get("rmarkdown.pandoc.to")
if (is.null(out_type))
    out_type <- "html"

## add styling
if (out_type == "html") {
    BiocStyle::markdown()
    ## BiocStyle::markdown(css.files = c('custom.css'))
} else if (out_type == "latex") {
    BiocStyle::latex()
}

## ----library------------------------------------------------------------------
library(seqLogo) 

## ----makePWM------------------------------------------------------------------
mFile <- system.file("extdata/pwm1", package="seqLogo")
m <- read.table(mFile)
m
p <- makePWM(m)

## ----slots--------------------------------------------------------------------
slotNames(p)
pwm(p)
ic(p)
consensus(p)

## ----seqLogo1, fig.height=4, fig.width=6, fig.cap="Sequence logo with column heights proportional to information content."----
seqLogo(p)

## ----seqLogo2, fig.height=4, fig.width=6, fig.cap="Sequence logo with uniform column heights."----
seqLogo(p, ic.scale=FALSE)

## ----seqLogo3, fig.height=4, fig.width=6, fig.cap="Sequence logo with user specified colors."----
seqLogo(p, fill=c(A="#4daf4a", C="#377eb8", G="#ffd92f", T="#e41a1c"),
        ic.scale=FALSE)

## ----seqLogo4, fig.height=4, fig.width=6, fig.cap="RNA Sequence logo."--------
r <- makePWM(m, alphabet="RNA")
seqLogo(r, ic.scale=FALSE)

## ----session-info-------------------------------------------------------------
sessionInfo()

