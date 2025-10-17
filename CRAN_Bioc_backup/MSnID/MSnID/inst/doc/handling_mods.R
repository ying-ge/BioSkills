### R code from vignette source 'handling_mods.Rnw'

###################################################
### code chunk number 1: style
###################################################
BiocStyle::latex(use.unsrturl=FALSE)


###################################################
### code chunk number 2: libs
###################################################
library(MSnID)
library(xtable)
library(dplyr)
library(Biostrings)


###################################################
### code chunk number 3: handling_mods.Rnw:40-43
###################################################
m <- MSnID(".")
mzids <- system.file("extdata","phospho.mzid.gz",package="MSnID")
m <- read_mzIDs(m, mzids)


###################################################
### code chunk number 4: handling_mods.Rnw:51-53
###################################################
# to know the present mod masses
report_mods(m)


###################################################
### code chunk number 5: handling_mods.Rnw:66-69
###################################################
m <- add_mod_symbol(m, mod_mass="79.966330925", symbol="*")
x <- psms(m) %>% 
    distinct(modification, peptide, peptide_mod)


###################################################
### code chunk number 6: handling_mods.Rnw:74-85
###################################################
sel_idx <- c(1,6,10,13,19)
x <- x %>%
  `[`(sel_idx,) %>%
    xtable()
align(x) <- "lp{1.8in}p{1.8in}p{1.8in}"
print(x,
      include.rownames = FALSE,
      size="\\fontsize{8pt}{6pt}\\selectfont",
      add.to.row = list(pos = as.list(seq_along(sel_idx)), 
                        command = rep("\\\\[3pt]",length(sel_idx))),
      floating = FALSE)


###################################################
### code chunk number 7: handling_mods.Rnw:93-97
###################################################
m <- add_mod_symbol(m, mod_mass="229.1629", symbol="#")
m <- add_mod_symbol(m, mod_mass="57.021463735", symbol="^")
x <- psms(m) %>% 
    distinct(modification, peptide, peptide_mod)


###################################################
### code chunk number 8: handling_mods.Rnw:102-112
###################################################
x <- x %>% 
    `[`(sel_idx,) %>%
    xtable()
align(x) <- "lp{1.8in}p{1.8in}p{1.8in}"
print(x,
      include.rownames = FALSE,
      size="\\fontsize{8pt}{6pt}\\selectfont",
      add.to.row = list(pos = as.list(seq_along(sel_idx)), 
                        command = rep("\\\\[3pt]",length(sel_idx))),
      floating = FALSE)


###################################################
### code chunk number 9: handling_mods.Rnw:126-128
###################################################
fst_path <- system.file("extdata","for_phospho.fasta.gz",package="MSnID")
fst <- readAAStringSet(fst_path)


###################################################
### code chunk number 10: handling_mods.Rnw:134-135
###################################################
names(fst) <- sub("(^[^ ]*) .*$", "\1", names(fst))


###################################################
### code chunk number 11: handling_mods.Rnw:142-150
###################################################
m <- map_mod_sites(m, fst, 
                   accession_col = "accession", 
                   peptide_mod_col = "peptide_mod", 
                   mod_char = "*",
                   site_delimiter = "lower")

x <- psms(m) %>% 
  distinct(peptide_mod, SiteID)


###################################################
### code chunk number 12: handling_mods.Rnw:153-163
###################################################
x <- x %>% 
    `[`(sel_idx,) %>%
    xtable()
align(x) <- "lp{1.8in}p{3.6in}"
print(x,
      include.rownames = FALSE,
      size="\\fontsize{8pt}{6pt}\\selectfont",
      add.to.row = list(pos = as.list(seq_along(sel_idx)), 
                        command = rep("\\\\[3pt]",length(sel_idx))),
      floating = FALSE)


###################################################
### code chunk number 13: handling_mods.Rnw:183-185
###################################################
conv_tab <- fetch_conversion_table("Homo sapiens", "UNIPROT", "SYMBOL")
head(conv_tab)


###################################################
### code chunk number 14: handling_mods.Rnw:197-200
###################################################
head(accessions(m))
m <- remap_accessions(m, conv_tab, extraction_pttrn = "\\|([^|-]+)(-\\d+)?\\|")
head(accessions(m))


###################################################
### code chunk number 15: handling_mods.Rnw:208-214
###################################################
fst_path <- system.file("extdata","for_phospho.fasta.gz",package="MSnID")
fst_path_2 <- remap_fasta_entry_names(fst_path, conv_tab, "\\|([^|-]+)(-\\d+)?\\|")

library(Biostrings)
readAAStringSet(fst_path)
readAAStringSet(fst_path_2)


###################################################
### code chunk number 16: handling_mods.Rnw:220-229
###################################################
fst <- readAAStringSet(fst_path_2)
m <- map_mod_sites(m, fst, 
                   accession_col = "accession", 
                   peptide_mod_col = "peptide_mod", 
                   mod_char = "*",
                   site_delimiter = "lower")

x <- psms(m) %>% 
  distinct(peptide_mod, SiteID)


###################################################
### code chunk number 17: handling_mods.Rnw:232-242
###################################################
x <- x %>% 
    `[`(sel_idx,) %>%
    xtable()
align(x) <- "lp{1.8in}p{3.6in}"
print(x,
      include.rownames = FALSE,
      size="\\fontsize{8pt}{6pt}\\selectfont",
      add.to.row = list(pos = as.list(seq_along(sel_idx)), 
                        command = rep("\\\\[3pt]",length(sel_idx))),
      floating = FALSE)


###################################################
### code chunk number 18: handling_mods.Rnw:248-249
###################################################
unlink(".Rcache", recursive=TRUE)


