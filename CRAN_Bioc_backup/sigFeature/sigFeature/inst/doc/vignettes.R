## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----cars---------------------------------------------------------------------
library(sigFeature)
library(SummarizedExperiment)
data(ExampleRawData, package="sigFeature")
ExampleRawData   

## -----------------------------------------------------------------------------
x  <- t(assays(ExampleRawData)$counts)
y  <- colData(ExampleRawData)$sampleLabels

## -----------------------------------------------------------------------------
pvals <- sigFeaturePvalue(x,y)
hist(unlist(pvals),breaks=seq(0,0.08,0.0015),col="skyblue",
    xlab="p value",ylab="Frequency",main="")

## -----------------------------------------------------------------------------
#system.time(sigfeatureRankedList <- sigFeature(x, y))

## -----------------------------------------------------------------------------
data(sigfeatureRankedList)
print(sigfeatureRankedList[1:10])

## -----------------------------------------------------------------------------
library(e1071)
sigFeature.model=svm(x[ ,sigfeatureRankedList[1:1000]], y, 
                    type="C-classification", kernel="linear")
summary(sigFeature.model)

## -----------------------------------------------------------------------------
pred <- predict(sigFeature.model, x[ ,sigfeatureRankedList[1:1000]])
table(pred,y)

## -----------------------------------------------------------------------------
#system.time(featureRankedList <- svmrfeFeatureRanking(x, y))
data(featureRankedList)
print("Top 10 features are printed below:")
print(featureRankedList[1:10])

## -----------------------------------------------------------------------------
RFE.model=svm(x[ ,featureRankedList[1:1000]], y, 
            type="C-classification", kernel="linear")
summary(RFE.model)

## -----------------------------------------------------------------------------
pred <- predict(RFE.model, x[ ,featureRankedList[1:1000]])
table(pred,y)

## -----------------------------------------------------------------------------
pvalsigFe <- sigFeaturePvalue(x, y, 100, sigfeatureRankedList)
pvalRFE   <- sigFeaturePvalue(x, y, 100, featureRankedList)
par(mfrow=c(1,2))
hist(unlist(pvalsigFe),breaks=50, col="skyblue", main=paste("sigFeature"), 
    xlab="p value")
hist(unlist(pvalRFE),breaks=50, col="skyblue", 
    main=paste("SVM-RFE"), xlab="p value")

## -----------------------------------------------------------------------------
mytitle<-'Box Plot'
boxplot(unlist(pvalsigFe), unlist(pvalRFE), main=mytitle, 
        names=c("sigFeature", "SVM-RFE"),
        ylab="p value", ylim=c(min(unlist(pvalsigFe)), max(unlist(pvalRFE))))
stripchart(unlist(pvalsigFe), vertical=TRUE, method="jitter", add=TRUE, pch=16, 
        col=c('green'))
stripchart(unlist(pvalRFE), vertical=TRUE, at=2, method="jitter", add=TRUE, 
        pch=16, col=c('blue'))
grid(nx=NULL, ny=NULL, col="black", lty="dotted")

## -----------------------------------------------------------------------------
library("pheatmap")
library("RColorBrewer")
pheatmap(x[ ,sigfeatureRankedList[1:20]], scale="row", 
        clustering_distance_rows="correlation")

## -----------------------------------------------------------------------------
pheatmap(x[ ,featureRankedList[1:20]], scale="row", 
        clustering_distance_rows="correlation")

## -----------------------------------------------------------------------------
#set.seed(1234)
#results = sigFeature.enfold(x, y, "kfold", 10)
data("results")
str(results[1])

## -----------------------------------------------------------------------------
FeatureBasedonFrequency <- sigFeatureFrequency(x, results, 400, 400, pf=FALSE)
str(FeatureBasedonFrequency[1])

## -----------------------------------------------------------------------------
#inputdata <- data.frame(y=as.factor(y) ,x=x)
#To run the code given bellow will take huge time. Thus the process 
#data is given below. 
#featsweepSigFe = lapply(1:400, sigCVError, FeatureBasedonFrequency, inputdata)
data("featsweepSigFe")
str(featsweepSigFe[1])

## -----------------------------------------------------------------------------
PlotErrors(featsweepSigFe, 0, 0.4)

## -----------------------------------------------------------------------------
#WritesigFeature(results, x)

## -----------------------------------------------------------------------------
sessionInfo()

