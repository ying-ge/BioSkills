### R code from vignette source 'R4RNA.Rnw'

###################################################
### code chunk number 1: R4RNA.Rnw:40-54
###################################################
library(R4RNA)

message("TRANSAT prediction in helix format")
transat_file <- system.file("extdata", "helix.txt", package = "R4RNA")
transat <- readHelix(transat_file)

message("RFAM structure in dot bracket format")
known_file <- system.file("extdata", "vienna.txt", package = "R4RNA")
known <- readVienna(known_file)

message("Work with basepairs instead of helices for more flexibility")
message("Breaks all helices into helices of length 1")
transat <- expandHelix(transat)
known <- expandHelix(known)


###################################################
### code chunk number 2: R4RNA.Rnw:63-65
###################################################
plotHelix(known, line = TRUE, arrow = TRUE)
mtext("Known Structure", side = 3, line = -2, adj = 0)


###################################################
### code chunk number 3: R4RNA.Rnw:73-76
###################################################
plotDoubleHelix(transat, known, line = TRUE, arrow = TRUE)
mtext("TRANSAT\nPredicted\nStructure", side = 3, line = -5, adj = 0)
mtext("Known Structure", side = 1, line = -2, adj = 0)


###################################################
### code chunk number 4: R4RNA.Rnw:84-86
###################################################
message("Filter out helices above a certain p-value")
transat <- transat[which(transat$value <= 1e-3), ]


###################################################
### code chunk number 5: R4RNA.Rnw:92-100
###################################################
message("Assign colour to basepairs according to p-value")
transat$col <- col <- colourByValue(transat, log = TRUE)

message("Coloured encoded in 'col' column of transat structure")
plotDoubleHelix(transat, known, line = TRUE, arrow = TRUE)

legend("topright", legend = attr(col, "legend"), fill = attr(col, "fill"),
    inset = 0.05, bty = "n", border = NA, cex = 0.75, title = "TRANSAT P-values")


###################################################
### code chunk number 6: R4RNA.Rnw:112-113
###################################################
plotOverlapHelix(transat, known, line = TRUE, arrow = TRUE, scale = FALSE)


###################################################
### code chunk number 7: R4RNA.Rnw:130-137
###################################################
message("Multiple sequence alignment of interest")
library(Biostrings)
fasta_file <- system.file("extdata", "fasta.txt", package = "R4RNA")
fasta <- as.character(readBStringSet(fasta_file))

message("Plot covariance in alignment")
plotCovariance(fasta, known, cex = 0.5)


###################################################
### code chunk number 8: R4RNA.Rnw:147-148
###################################################
plotCovariance(fasta, transat, cex = 0.5, conflict.col = "grey")


###################################################
### code chunk number 9: R4RNA.Rnw:158-162
###################################################
col <- colourByCovariation(known, fasta, get = TRUE)
plotCovariance(fasta, col, grid = TRUE, legend = FALSE)
legend("topright", legend = attr(col, "legend"), fill = attr(col, "fill"),
    inset = 0.1, bty = "n", border = NA, cex = 0.37, title = "Covariation")


###################################################
### code chunk number 10: R4RNA.Rnw:168-173
###################################################
custom_colours <- c("green", "blue", "cyan", "red", "black", "grey")
plotCovariance(fasta, col <- colourByConservation(known, fasta, get = TRUE),
    palette = custom_colours, cex = 0.5)
legend("topright", legend = attr(col, "legend"), fill = attr(col, "fill"),
    inset = 0.15, bty = "n", border = NA, cex = 0.75, title = "Conservation")


###################################################
### code chunk number 11: R4RNA.Rnw:178-182
###################################################
col <- colourByCanonical(known, fasta, custom_colours, get = TRUE)
plotCovariance(fasta, col, base.colour = TRUE, cex = 0.5)
legend("topright", legend = attr(col, "legend"), fill = attr(col, "fill"),
    inset = 0.15, bty = "n", border = NA, cex = 0.75, title = "% Canonical")


###################################################
### code chunk number 12: R4RNA.Rnw:187-189
###################################################
col <- colourByUnknottedGroups(known, c("red", "blue"), get = TRUE)
plotCovariance(fasta, col, base.colour = TRUE, legend = FALSE, species = 23, grid = TRUE, text = TRUE, text.cex = 0.2, cex = 0.5)


###################################################
### code chunk number 13: R4RNA.Rnw:197-199
###################################################
toLatex(sessionInfo())



