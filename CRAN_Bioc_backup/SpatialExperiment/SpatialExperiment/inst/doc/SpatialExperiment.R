## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(cache = TRUE, autodep = TRUE, cache.lazy = FALSE)

## ----echo=FALSE, out.width = "100%", fig.cap="Overview of the SpatialExperiment class structure."----
knitr::include_graphics("SPE.png")

## ----message = FALSE----------------------------------------------------------
library(SpatialExperiment)
example(read10xVisium, echo = FALSE)
spe

## -----------------------------------------------------------------------------
head(spatialCoords(spe))

## -----------------------------------------------------------------------------
spatialCoordsNames(spe)

## -----------------------------------------------------------------------------
imgData(spe)

## -----------------------------------------------------------------------------
(spi <- getImg(spe))
identical(spi, imgData(spe)$data[[1]])

## ----fig.small = TRUE---------------------------------------------------------
plot(imgRaster(spe))

## ----fig.small = TRUE, eval = TRUE--------------------------------------------
url <- "https://i.redd.it/3pw5uah7xo041.jpg"
spe <- addImg(spe, 
    sample_id = "section1", 
    image_id = "pomeranian",
    imageSource = url, 
    scaleFactor = NA_real_, 
    load = TRUE)
img <- imgRaster(spe, 
    sample_id = "section1", 
    image_id = "pomeranian")
plot(img)

## -----------------------------------------------------------------------------
imgData(spe <- rmvImg(spe, "section1", "pomeranian"))

## -----------------------------------------------------------------------------
n <- length(z <- letters)
y <- matrix(nrow = n, ncol = n)
cd <- DataFrame(x = seq(n), y = seq(n), z)

spe1 <- SpatialExperiment(
    assay = y, 
    colData = cd, 
    spatialCoordsNames = c("x", "y"))

## -----------------------------------------------------------------------------
xy <- as.matrix(cd[, c("x", "y")])

spe2 <- SpatialExperiment(
    assay = y, 
    colData = cd["z"], 
    spatialCoords = xy)

## -----------------------------------------------------------------------------
identical(spe1, spe2)

## -----------------------------------------------------------------------------
spe <- SpatialExperiment(
    assays = y)
isEmpty(spatialCoords(spe))

## ----results = "hide"---------------------------------------------------------
n <- 10; m <- 20
y <- matrix(nrow = n, ncol = m)
cd <- DataFrame(x = seq(m), y = seq(m))
xy <- matrix(nrow = m, ncol = 2)
colnames(xy) <- c("x", "y")

SpatialExperiment(
    assay = y, 
    colData = cd,
    spatialCoordsNames = c("x", "y"),
    spatialCoords = xy)

## -----------------------------------------------------------------------------
dir <- system.file(
   file.path("extdata", "10xVisium", "section1", "outs"),
   package = "SpatialExperiment")

# read in counts
fnm <- file.path(dir, "raw_feature_bc_matrix")
sce <- DropletUtils::read10xCounts(fnm)

# read in image data
img <- readImgData(
    path = file.path(dir, "spatial"),
    sample_id = "foo")

# read in spatial coordinates
fnm <- file.path(dir, "spatial", "tissue_positions_list.csv")
xyz <- read.csv(fnm, header = FALSE,
    col.names = c(
        "barcode", "in_tissue", "array_row", "array_col",
        "pxl_row_in_fullres", "pxl_col_in_fullres"))

# construct observation & feature metadata
rd <- S4Vectors::DataFrame(
    symbol = rowData(sce)$Symbol)

# construct 'SpatialExperiment'
(spe <- SpatialExperiment(
    assays = list(counts = assay(sce)),
    rowData = rd, 
    colData = DataFrame(xyz), 
    spatialCoordsNames = c("pxl_col_in_fullres", "pxl_row_in_fullres"),
    imgData = img,
    sample_id = "foo"))

## -----------------------------------------------------------------------------
dir <- system.file(
    file.path("extdata", "10xVisium"),
    package = "SpatialExperiment")

sample_ids <- c("section1", "section2")
samples <- file.path(dir, sample_ids, "outs")

(spe10x <- read10xVisium(samples, sample_ids,
    type = "sparse", data = "raw",
    images = "lowres", load = FALSE))

## -----------------------------------------------------------------------------
samples2 <- file.path(dir, sample_ids)

(spe10x2 <- read10xVisium(samples2, sample_ids,
    type = "sparse", data = "raw",
    images = "lowres", load = FALSE))

## ----message = FALSE, warning = FALSE-----------------------------------------
n <- 1e3  # number of molecules
ng <- 50  # number of genes
nc <- 20  # number of cells
# sample xy-coordinates in [0, 1]
x <- runif(n)
y <- runif(n)
# assign each molecule to some gene-cell pair
gs <- paste0("gene", seq(ng))
cs <- paste0("cell", seq(nc))
gene <- sample(gs, n, TRUE)
cell <- sample(cs, n, TRUE)
# assure gene & cell are factors so that
# missing observations aren't dropped
gene <- factor(gene, gs)
cell <- factor(cell, cs)
# construct data.frame of molecule coordinates
df <- data.frame(gene, cell, x, y)
head(df)

## ----message = FALSE, warning = FALSE-----------------------------------------
# construct 'BumpyMatrix'
library(BumpyMatrix)
mol <- splitAsBumpyMatrix(
    df[, c("x", "y")], 
    row = gene, col = cell)

## ----message = FALSE, warning = FALSE-----------------------------------------
# get count matrix
y <- with(df, table(gene, cell))
y <- as.matrix(unclass(y))
y[1:5, 1:5]
# construct SpatialExperiment
spe <- SpatialExperiment(
    assays = list(
        counts = y, 
        molecules = mol))
spe

## ----message = FALSE, warning = FALSE-----------------------------------------
molecules(spe)

## -----------------------------------------------------------------------------
sub <- spe10x[, spe10x$sample_id == "section1"]

## -----------------------------------------------------------------------------
sub <- spe10x[, colData(spe10x)$in_tissue]
sum(colData(spe10x)$in_tissue) == ncol(sub)

## -----------------------------------------------------------------------------
spe1 <- spe2 <- spe
spe3 <- cbind(spe1, spe2)
unique(spe3$sample_id)

## -----------------------------------------------------------------------------
# make sample identifiers unique
spe1 <- spe2 <- spe
spe1$sample_id <- paste(spe1$sample_id, "A", sep = ".")
spe2$sample_id <- paste(spe2$sample_id, "B", sep = ".")

# combine into single object
spe3 <- cbind(spe1, spe2)

## ----error=TRUE---------------------------------------------------------------
new <- spe3$sample_id
new[1] <- "section2.A"
spe3$sample_id <- new
new[1] <- "third.one.of.two"
spe3$sample_id <- new

## -----------------------------------------------------------------------------
# backup original sample IDs
tmp <- spe$sample_id
# try to remove sample IDs
spe$sample_id <- NULL
# sample IDs remain unchanged
identical(tmp, spe$sample_id)

## -----------------------------------------------------------------------------
# extract first image
spi <- getImg(spe10x)
# apply counter-/clockwise rotation
spi1 <- rotateImg(spi, -90)
spi2 <- rotateImg(spi, +90)
# visual comparison
par(mfrow = c(1, 3))
plot(as.raster(spi))
plot(as.raster(spi1))
plot(as.raster(spi2))

## -----------------------------------------------------------------------------
# specify sample & image identifier
sid <- "section1"
iid <- "lowres"
# counter-clockwise rotation
tmp <- rotateImg(spe10x, 
    sample_id = sid, 
    image_id = iid,
    degrees = -90)
# visual comparison
par(mfrow = c(1, 2))
plot(imgRaster(spe10x, sid, iid))
plot(imgRaster(tmp, sid, iid))

## -----------------------------------------------------------------------------
# extract first image
spi <- getImg(spe10x)
# mirror horizontally/vertically
spi1 <- mirrorImg(spi, "h")
spi2 <- mirrorImg(spi, "v")
# visual comparison
par(mfrow = c(1, 3))
plot(as.raster(spi))
plot(as.raster(spi1))
plot(as.raster(spi2))

## -----------------------------------------------------------------------------
# specify sample & image identifier
sid <- "section2"
iid <- "lowres"
# mirror horizontally
tmp <- mirrorImg(spe10x, 
    sample_id = sid, 
    image_id = iid,
    axis = "h")
# visual comparison
par(mfrow = c(1, 2))
plot(imgRaster(spe10x, sid, iid))
plot(imgRaster(tmp, sid, iid))

## ----tidy = TRUE--------------------------------------------------------------
sessionInfo()

