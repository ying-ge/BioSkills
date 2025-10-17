## ----style, echo=FALSE, results="asis", message=FALSE-------------------------
knitr::opts_chunk$set(tidy = FALSE,
                      warning = FALSE,
                      message = FALSE)

## ----echo=FALSE, results='hide', message=FALSE--------------------------------
CRANpkg <- function (pkg) {
    cran <- "https://CRAN.R-project.org/package"
    fmt <- "[%s](%s=%s)"
    sprintf(fmt, pkg, cran, pkg)
}

Biocpkg <- function (pkg) {
    sprintf("[%s](http://bioconductor.org/packages/%s)", pkg, pkg)
}

library(org.Hs.eg.db)
library(DOSE)
library(ReactomePA)

