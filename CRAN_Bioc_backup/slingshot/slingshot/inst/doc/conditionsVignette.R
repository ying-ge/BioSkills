## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(SingleCellExperiment, quietly = TRUE)
library(slingshot, quietly = TRUE)
library(RColorBrewer)
graphics:::par(pch = 16, las = 1)
set.seed(1)

## ----dataSetup----------------------------------------------------------------
data('slingshotExample')
rd <- slingshotExample$rd
cl <- slingshotExample$cl
condition <- factor(rep(c('A','B'), length.out = nrow(rd)))
condition[110:140] <- 'A'
ls()

## ----echo=FALSE---------------------------------------------------------------
plot(rd, asp = 1, pch = 16, col = brewer.pal(3,'Set1')[condition], las=1)
legend('topleft','(x,y)',legend = c('A','B'), title = 'Condition', pch=16, col = brewer.pal(3,'Set1')[1:2])

## -----------------------------------------------------------------------------
pto <- slingshot(rd, cl)

## ----echo=FALSE---------------------------------------------------------------
plot(rd, asp = 1, pch = 16, col = brewer.pal(3,'Set1')[condition], las=1)
lines(SlingshotDataSet(pto), lwd=3)
legend('topleft','(x,y)',legend = c('A','B'), title = 'Condition', pch=16, col = brewer.pal(3,'Set1')[1:2])

## ----fig.height=4, fig.width=8, echo = FALSE----------------------------------
n <- nrow(rd); L <- ncol(slingPseudotime(pto))
noise <- runif(n, -.1,.1)
plot(as.numeric(slingPseudotime(pto)), rep(1:L, each=n)+noise,pch=16, col = brewer.pal(9,'Set1')[condition], axes=FALSE, xlab='Pseudotime', ylab='Lineage', las=1)
axis(1); axis(2, at=1:L, las=1)

## ----echo=FALSE---------------------------------------------------------------
# tA1 <- slingPseudotime(pto, na=FALSE)[condition=='A',1]
# wA1 <- slingCurveWeights(pto)[condition=='A',1]; wA1 <- wA1/sum(wA1)
# dA1 <- density(tA1, weights = wA1)
# tB1 <- slingPseudotime(pto, na=FALSE)[condition=='B',1]
# wB1 <- slingCurveWeights(pto)[condition=='B',1]; wB1 <- wB1/sum(wB1)
# dB1 <- density(tB1, weights = wB1)
# tA2 <- slingPseudotime(pto, na=FALSE)[condition=='A',2]
# wA2 <- slingCurveWeights(pto)[condition=='A',2]; wA2 <- wA2/sum(wA2)
# dA2 <- density(tA2, weights = wA2)
# tB2 <- slingPseudotime(pto, na=FALSE)[condition=='B',2]
# wB2 <- slingCurveWeights(pto)[condition=='B',2]; wB2 <- wB2/sum(wB2)
# dB2 <- density(tB2, weights = wB2)
# 
# plot(range(slingPseudotime(pto),na.rm=TRUE), c(1,2+7*max(c(dA2$y,dB2$y))), col='white', xlab='Pseudotime', ylab='Lineage', axes = FALSE, las=1)
# axis(1); axis(2, at=1:2)
# polygon(c(min(dA1$x),dA1$x,max(dA1$x)), 1+7*c(0,dA1$y,0), col=rgb(1,0,0,.5))
# polygon(c(min(dB1$x),dB1$x,max(dB1$x)), 1+7*c(0,dB1$y,0), col=rgb(0,0,1,.5))
# polygon(c(min(dA2$x),dA2$x,max(dA2$x)), 2+7*c(0,dA2$y,0), col=rgb(1,0,0,.5))
# polygon(c(min(dB2$x),dB2$x,max(dB2$x)), 2+7*c(0,dB2$y,0), col=rgb(0,0,1,.5))

layout(matrix(1:2,nrow=1))
boxplot(slingPseudotime(pto)[,1] ~ condition, col = brewer.pal(3,'Set1')[1:2], main = 'Lineage 1', xlab='Condition', ylab='Pseudotime', las=1, pch = 16)
boxplot(slingPseudotime(pto)[,2] ~ condition, col = brewer.pal(3,'Set1')[1:2], main = 'Lineage 2', xlab='Condition', ylab='Pseudotime', las=1, pch = 16)
layout(1)

## ----eval=FALSE---------------------------------------------------------------
#  # Permutation test
#  t1 <- slingPseudotime(pto, na=FALSE)[,1]
#  w1 <- slingCurveWeights(pto)[,1]
#  d1 <- weighted.mean(t1[condition=='A'], w1[condition=='A']) -
#      weighted.mean(t1[condition=='B'], w1[condition=='B'])
#  dist1 <- replicate(1e4, {
#      condition.i <- sample(condition)
#      weighted.mean(t1[condition.i=='A'], w1[condition.i=='A']) -
#          weighted.mean(t1[condition.i=='B'], w1[condition.i=='B'])
#  })

## ----echo=FALSE, fig.height=4, fig.width=9------------------------------------
t1 <- slingPseudotime(pto, na=FALSE)[,1]
w1 <- slingCurveWeights(pto)[,1]
d1 <- weighted.mean(t1[condition=='A'], w1[condition=='A']) - 
    weighted.mean(t1[condition=='B'], w1[condition=='B'])
dist1 <- replicate(1e4, {
    condition.i <- sample(condition)
    weighted.mean(t1[condition.i=='A'], w1[condition.i=='A']) - 
        weighted.mean(t1[condition.i=='B'], w1[condition.i=='B'])
})
t2 <- slingPseudotime(pto, na=FALSE)[,2]
w2 <- slingCurveWeights(pto)[,2]
d2 <- weighted.mean(t2[condition=='A'], w2[condition=='A']) - 
    weighted.mean(t2[condition=='B'], w2[condition=='B'])
dist2 <- replicate(1e4, {
    condition.i <- sample(condition)
    weighted.mean(t2[condition.i=='A'], w2[condition.i=='A']) - 
        weighted.mean(t2[condition.i=='B'], w2[condition.i=='B'])
})

layout(matrix(1:2,nrow = 1))
hist(dist1, breaks=50, col='grey50', xlim = range(c(d1,dist1)), probability = TRUE, xlab = 'Difference of Weighted Means', main = 'Lineage 1 Permutation Test', las=1)
abline(v=d1, col=2, lwd=2)
legend('topright','(x,y)',legend = c('Null Distn.','Observed'), fill=c('grey50',NA), border=c(1,NA), lty=c(NA,1), lwd=c(NA,2), col=c(NA,2), merge = TRUE, seg.len = .6)

hist(dist2, breaks=50, col='grey50', xlim = range(c(d2,dist2)), probability = TRUE, xlab = 'Difference of Weighted Means', main = 'Lineage 2 Permutation Test', las=1)
abline(v=d2, col=2, lwd=2)
legend('topright','(x,y)',legend = c('Null Distn.','Observed'), fill=c('grey50',NA), border=c(1,NA), lty=c(NA,1), lwd=c(NA,2), col=c(NA,2), merge = TRUE, seg.len = .6)
layout(1)

## -----------------------------------------------------------------------------
paste0('Lineage 1 p-value: ', mean(abs(dist1) > abs(d1)))
paste0('Lineage 2 p-value: ', mean(abs(dist2) > abs(d2)))

## -----------------------------------------------------------------------------
# Kolmogorov-Smirnov test
ks.test(slingPseudotime(pto)[condition=='A',1], slingPseudotime(pto)[condition=='B',1])
ks.test(slingPseudotime(pto)[condition=='A',2], slingPseudotime(pto)[condition=='B',2])

## -----------------------------------------------------------------------------
pt <- slingPseudotime(pto, na=FALSE)
cw <- slingCurveWeights(pto)
assignment <- apply(cw, 1, which.max)
ptAs <- c() #assigned pseudotime
for(ii in 1:nrow(pt)) ptAs[ii] <- pt[ii,assignment[ii]] 

## ----echo=FALSE,fig.height=3.5------------------------------------------------
layout(matrix(1:2,nrow=1))
noise <- runif(n, -.05,.05)
plot(ptAs[assignment == 1], (as.numeric(condition)+noise)[assignment == 1], 
     col = brewer.pal(9,'Set1')[condition[assignment == 1]],
     xlab="Pseudotime", ylab="Condition", main="Lineage 1", pch=16, axes = FALSE)
axis(1); axis(2, at=seq_along(levels(condition)), labels = levels(condition), las = 1)

plot(ptAs[assignment == 2], (as.numeric(condition)+noise)[assignment == 2], 
     col = brewer.pal(9,'Set1')[condition[assignment == 2]],
     xlab="Pseudotime", ylab="Condition", main="Lineage 2", pch=16, axes = FALSE)
axis(1); axis(2, at=seq_along(levels(condition)), labels = levels(condition), las = 1)
layout(1)

## -----------------------------------------------------------------------------
# model for lineage 1: not significant
cond1 <- factor(condition[assignment == 1])
t1 <- ptAs[assignment == 1]
m1 <- glm(cond1 ~ t1, family = quasibinomial(link = "logit"))
summary(m1)

## -----------------------------------------------------------------------------
# model for lineage 2: significant
cond2 <- factor(condition[assignment == 2])
t2 <- ptAs[assignment == 2]
m2 <- glm(cond2 ~ t2, family = quasibinomial(link = "logit"))
summary(m2)

## -----------------------------------------------------------------------------
### note that logistic regression is monotone hence only allows for increasing or decreasing proportion of cells for a particular condition.
### hence, for complex trajectories, you could consider smoothing the pseudotime, i.e.
require(mgcv, quietly = TRUE)
m1s <- mgcv::gam(cond1 ~ s(t1), family="quasibinomial")
summary(m1s)

## -----------------------------------------------------------------------------
m2s <- mgcv::gam(cond2 ~ s(t2), family="quasibinomial")
summary(m2s)

## ----echo=FALSE---------------------------------------------------------------
layout(matrix(1:2,nrow=1))
plot(m1s, shade = TRUE, main = "Lineage 1", 
     xlab="Pseudotime", ylab="Logit Prob(B)", scheme=0)
plot(m2s, shade = TRUE, main = "Lineage 2",
     xlab="Pseudotime", ylab="Logit Prob(B)")
layout(1)

## -----------------------------------------------------------------------------
t1 <- pt[,1]
t2 <- pt[,2]
l1 <- as.numeric(assignment == 1)
l2 <- as.numeric(assignment == 2)
m <- gam(condition ~ s(t1, by=l1, id=1) + s(t2, by=l2, id=1),
         family = quasibinomial(link = "logit"))
summary(m)
### and then we're back to tradeSeq-like inference ...

## ----echo=FALSE---------------------------------------------------------------
layout(matrix(1:2,nrow=1))
plot(m, select=1, shade=TRUE, main='Lineage 1',
     xlab="Pseudotime", ylab="Logit Prob(B)")
plot(m, select=2, shade=TRUE, main='Lineage 2',
     xlab="Pseudotime", ylab="Logit Prob(B)")
layout(1)

## ----session------------------------------------------------------------------
sessionInfo()

