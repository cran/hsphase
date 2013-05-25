### R code from vignette source 'hsphase.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: imageplot
###################################################

set.seed(718)
library(hsphase)
halfsibs <- .simulateHalfsib(numInd = 20)
imageplot(bmh(halfsibs),title = "Imageplot of simulated half-sib family")


###################################################
### code chunk number 2: rplot
###################################################
rplot(halfsibs,sort(sample(c(1:(1000*ncol(halfsibs))),size=ncol(halfsibs))))
title("Recombination of simulated half-sib family")


###################################################
### code chunk number 3: heatmap
###################################################
a <- .simulateHalfsib()
b <- .simulateHalfsib()
d <- rbind(a,b)
library(hsphase)
oh <- ohg(d)
heatmap(oh,symm=T,col=gray.colors(16,start=0,end=1),RowSideColors=as.character(c(rep(1,40),rep(2,40))),ColSideColors=as.character(c(rep(1,5),rep(2,40),rep(1,35))))


