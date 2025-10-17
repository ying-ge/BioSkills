## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width=6, fig.height=4
)

## -----------------------------------------------------------------------------
library(biglasso)

data(colon)
X <- colon$X
y <- colon$y
dim(X)

## -----------------------------------------------------------------------------
X[1:5, 1:5]

## -----------------------------------------------------------------------------
## convert X to a big.matrix object
## X.bm is a pointer to the data matrix
X.bm <- as.big.matrix(X)
str(X.bm)

## -----------------------------------------------------------------------------
dim(X.bm)

## -----------------------------------------------------------------------------
X.bm[1:5, 1:5]
## same results as X[1:5, 1:5]

## -----------------------------------------------------------------------------
## fit entire solution path, using our newly proposed "Adaptive" screening rule (default)
fit <- biglasso(X.bm, y)
plot(fit)

## -----------------------------------------------------------------------------
## 10-fold cross-valiation in parallel
cvfit <- tryCatch(
         {
                cv.biglasso(X.bm, y, seed = 1234, nfolds = 10, ncores = 4)
         },
         error = function(cond) {
                cv.biglasso(X.bm, y, seed = 1234, nfolds = 10, ncores = 2)
         }
)

## -----------------------------------------------------------------------------
par(mfrow = c(2, 2), mar = c(3.5, 3.5, 3, 1) ,mgp = c(2.5, 0.5, 0))
plot(cvfit, type = "all")

## -----------------------------------------------------------------------------
summary(cvfit)

## -----------------------------------------------------------------------------
coef(cvfit)[which(coef(cvfit) != 0),]

## -----------------------------------------------------------------------------
data(Heart)
X <- Heart$X
y <- Heart$y
X.bm <- as.big.matrix(X)
fit <- biglasso(X.bm, y, family = "binomial")
plot(fit)

## -----------------------------------------------------------------------------
library(survival)
X <- heart[,4:7]
y <- Surv(heart$stop - heart$start, heart$event)
X.bm <- as.big.matrix(X)
fit <- biglasso(X.bm, y, family = "cox")
plot(fit)

## -----------------------------------------------------------------------------
set.seed(10101)
n=300; p=300; m=5; s=10; b=1
x = matrix(rnorm(n * p), n, p)
beta = matrix(seq(from=-b,to=b,length.out=s*m),s,m)
y = x[,1:s] %*% beta + matrix(rnorm(n*m,0,1),n,m)
x.bm = as.big.matrix(x)
fit = biglasso(x.bm, y, family = "mgaussian")
plot(fit)

## -----------------------------------------------------------------------------
## The data has 200 observations, 600 features, and 10 non-zero coefficients.
## This is not actually very big, but vignettes in R are supposed to render
## quickly. Much larger data can be handled in the same way.
if(!file.exists('BigX.bin')) {
  X <- matrix(rnorm(1000 * 5000), 1000, 5000)
  beta <- c(-5:5)
  y <- as.numeric(X[,1:11] %*% beta)
  write.csv(X, "BigX.csv", row.names = F)
  write.csv(y, "y.csv", row.names = F)
  ## Pretend that the data in "BigX.csv" is too large to fit into memory
  X.bm <- setupX("BigX.csv", header = T)
}

## ----warning=F----------------------------------------------------------------
rm(list = c("X", "X.bm", "y")) # Pretend starting a new session
X.bm <- attach.big.matrix("BigX.desc")
y <- read.csv("y.csv")[,1]

## -----------------------------------------------------------------------------
system.time({fit <- biglasso(X.bm, y)})

## -----------------------------------------------------------------------------
plot(fit)

## -----------------------------------------------------------------------------
# 10-fold cross validation in parallel
tryCatch(
	{
		system.time({cvfit <- cv.biglasso(X.bm, y, seed = 1234, ncores = 4, nfolds = 10)})
	},
	error = function(cond) {
  	  	system.time({cvfit <- cv.biglasso(X.bm, y, seed = 1234, ncores = 2, nfolds = 10)})
	}
)

## -----------------------------------------------------------------------------
par(mfrow = c(2, 2), mar = c(3.5, 3.5, 3, 1), mgp = c(2.5, 0.5, 0))
plot(cvfit, type = "all")

