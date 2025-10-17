## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(tidyestimate)

## -----------------------------------------------------------------------------
dim(ov)

## -----------------------------------------------------------------------------
head(ov[, 1:5])

## -----------------------------------------------------------------------------
filtered <- filter_common_genes(ov, 
                                id = "hgnc_symbol", 
                                tidy = FALSE, 
                                tell_missing = TRUE, 
                                find_alias = TRUE)

## -----------------------------------------------------------------------------
scored <- estimate_score(filtered,
                         is_affymetrix = TRUE)
scored

## ---- fig.width=6, fig.height=6-----------------------------------------------
plot_purity(scored, is_affymetrix = TRUE)

