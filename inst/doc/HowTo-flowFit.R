### R code from vignette source 'HowTo-flowFit.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: CompanionPkg (eval = FALSE)
###################################################
## 
## source("http://www.bioconductor.org/biocLite.R")
## 
## biocLite("flowFitExampleData")


###################################################
### code chunk number 2: HowTo-flowFit.Rnw:259-261
###################################################
library(flowFitExampleData)
library(flowFit)


###################################################
### code chunk number 3: HowTo-flowFit.Rnw:266-267
###################################################
library(flowCore)


###################################################
### code chunk number 4: HowTo-flowFit.Rnw:273-274
###################################################
data(QuahAndParish)


###################################################
### code chunk number 5: HowTo-flowFit.Rnw:291-292
###################################################
QuahAndParish[[1]]


###################################################
### code chunk number 6: HowTo-flowFit.Rnw:309-312
###################################################
parent.fitting.cfse <- parentFitting(QuahAndParish[[1]], "<FITC-A>")
parent.fitting.cpd <- parentFitting(QuahAndParish[[1]], "<APC-A>")
parent.fitting.ctv <- parentFitting(QuahAndParish[[1]], "<Alexa Fluor 405-A>")


###################################################
### code chunk number 7: HowTo-flowFit.Rnw:318-320 (eval = FALSE)
###################################################
## parent.fitting.cfse <- parentFitting(QuahAndParish[[1]], "<FITC-A>")
## plot(parent.fitting.cfse)


###################################################
### code chunk number 8: fig1
###################################################
plot(parent.fitting.cfse)


###################################################
### code chunk number 9: HowTo-flowFit.Rnw:338-340 (eval = FALSE)
###################################################
## parent.fitting.cpd <- parentFitting(QuahAndParish[[1]], "<APC-A>")
## plot(parent.fitting.cpd)


###################################################
### code chunk number 10: fig2
###################################################
plot(parent.fitting.cpd)


###################################################
### code chunk number 11: HowTo-flowFit.Rnw:357-359 (eval = FALSE)
###################################################
## parent.fitting.ctv <- parentFitting(QuahAndParish[[1]], "<Alexa Fluor 405-A>")
## plot(parent.fitting.ctv)


###################################################
### code chunk number 12: fig3
###################################################
plot(parent.fitting.ctv)


###################################################
### code chunk number 13: HowTo-flowFit.Rnw:379-380
###################################################
summary(parent.fitting.cfse)


###################################################
### code chunk number 14: HowTo-flowFit.Rnw:385-386
###################################################
confint(parent.fitting.cfse)


###################################################
### code chunk number 15: HowTo-flowFit.Rnw:396-397
###################################################
QuahAndParish[[1]]@parameters@data[7,]


###################################################
### code chunk number 16: HowTo-flowFit.Rnw:402-404 (eval = FALSE)
###################################################
## logTrans <- logTransform(transformationId="log10-transformation",
##  logbase=10, r=1, d=1)


###################################################
### code chunk number 17: HowTo-flowFit.Rnw:414-416
###################################################
data(PKH26data)
keyword(PKH26data[[1]])$`$P4M`


###################################################
### code chunk number 18: HowTo-flowFit.Rnw:421-423
###################################################
acquisiton.resolution <- QuahAndParish[[1]]@parameters@data$range[7]
log10(acquisiton.resolution)


###################################################
### code chunk number 19: HowTo-flowFit.Rnw:440-442
###################################################
parent.fitting.cfse@parentPeakPosition
parent.fitting.cfse@parentPeakSize


###################################################
### code chunk number 20: HowTo-flowFit.Rnw:447-451
###################################################
fitting.cfse <- proliferationFitting(QuahAndParish[[2]],
  "<FITC-A>", parent.fitting.cfse@parentPeakPosition,
  parent.fitting.cfse@parentPeakSize)
fitting.cfse


###################################################
### code chunk number 21: HowTo-flowFit.Rnw:455-456
###################################################
coef(fitting.cfse)


###################################################
### code chunk number 22: HowTo-flowFit.Rnw:459-460
###################################################
summary(fitting.cfse)


###################################################
### code chunk number 23: HowTo-flowFit.Rnw:488-490 (eval = FALSE)
###################################################
## library(flowFit)
## ?plot


###################################################
### code chunk number 24: HowTo-flowFit.Rnw:501-502 (eval = FALSE)
###################################################
## plot(fitting.cfse, which=1)


###################################################
### code chunk number 25: fig4
###################################################
plot(fitting.cfse, which=1)


###################################################
### code chunk number 26: HowTo-flowFit.Rnw:525-526 (eval = FALSE)
###################################################
## plot(fitting.cfse, which=2)


###################################################
### code chunk number 27: fig5
###################################################
plot(fitting.cfse, which=2)


###################################################
### code chunk number 28: HowTo-flowFit.Rnw:549-550 (eval = FALSE)
###################################################
## plot(fitting.cfse, which=3)


###################################################
### code chunk number 29: fig6
###################################################
plot(fitting.cfse, which=3)


###################################################
### code chunk number 30: HowTo-flowFit.Rnw:572-573 (eval = FALSE)
###################################################
## plot(fitting.cfse, which=4)


###################################################
### code chunk number 31: fig7
###################################################
plot(fitting.cfse, which=4)


###################################################
### code chunk number 32: HowTo-flowFit.Rnw:602-603 (eval = FALSE)
###################################################
## plot(fitting.cfse, which=5)


###################################################
### code chunk number 33: fig8
###################################################
plot(fitting.cfse, which=5)


###################################################
### code chunk number 34: HowTo-flowFit.Rnw:625-627 (eval = FALSE)
###################################################
## par(mfrow=c(2,2))
## plot(fitting.cfse, which=1:4, legend=FALSE)


###################################################
### code chunk number 35: fig9
###################################################
par(mfrow=c(2,2))
plot(fitting.cfse, which=1:4, legend=FALSE)


###################################################
### code chunk number 36: HowTo-flowFit.Rnw:641-643 (eval = FALSE)
###################################################
## par(mfrow=c(2,1))
## plot(fitting.cfse, which=c(3,5), legend=FALSE)


###################################################
### code chunk number 37: fig10
###################################################
par(mfrow=c(2,1))
plot(fitting.cfse, which=c(3,5), legend=FALSE)


###################################################
### code chunk number 38: GenFittingOut
###################################################
gen <- getGenerations(fitting.cfse)
class(gen)


###################################################
### code chunk number 39: GenFittingOutVector
###################################################
fitting.cfse@generations


###################################################
### code chunk number 40: ProliferationIndex
###################################################
proliferationIndex(fitting.cfse)


###################################################
### code chunk number 41: FittingCTVQuah
###################################################
fitting.ctv <- proliferationFitting(QuahAndParish[[4]],
  "<Alexa Fluor 405-A>", parent.fitting.ctv@parentPeakPosition,
  parent.fitting.ctv@parentPeakSize)


###################################################
### code chunk number 42: HowTo-flowFit.Rnw:739-740 (eval = FALSE)
###################################################
## plot(fitting.ctv, which=3)


###################################################
### code chunk number 43: fig11
###################################################
plot(fitting.ctv, which=3)


###################################################
### code chunk number 44: HowTo-flowFit.Rnw:760-761 (eval = FALSE)
###################################################
## plot(QuahAndParish[[3]],"<APC-A>", breaks=1024, main="CPD sample")


###################################################
### code chunk number 45: fig12
###################################################
plot(QuahAndParish[[3]],"<APC-A>", breaks=1024, main="CPD sample")


###################################################
### code chunk number 46: FittingCPDQuah
###################################################
fitting.cpd <- proliferationFitting(QuahAndParish[[3]],
  "<APC-A>",
  parent.fitting.cpd@parentPeakPosition,
  parent.fitting.cpd@parentPeakSize,
  fixedModel=TRUE,
  fixedPars=list(M=parent.fitting.cpd@parentPeakPosition))


###################################################
### code chunk number 47: HowTo-flowFit.Rnw:788-789 (eval = FALSE)
###################################################
## plot(fitting.cpd, which=3)


###################################################
### code chunk number 48: fig13
###################################################
plot(fitting.cpd, which=3)


###################################################
### code chunk number 49: CHISQUARETEST
###################################################
perc.cfse <- fitting.cfse@generations
perc.cpd <- fitting.cpd@generations
perc.ctv <- fitting.ctv@generations


###################################################
### code chunk number 50: CHISQUARETESTTRIM
###################################################
perc.cfse <- c(perc.cfse, rep(0,6))


###################################################
### code chunk number 51: CHISQUARE_TEST
###################################################
M <- rbind(perc.cfse, perc.cpd, perc.ctv)
colnames(M) <- 1:16
(Xsq <- chisq.test(M, B=100000, simulate.p.value=TRUE))


###################################################
### code chunk number 52: HowTo-flowFit.Rnw:827-834 (eval = FALSE)
###################################################
## plot(perc.cfse, type="b", axes=F, ylim=c(0,50), xlab="generations", ylab="Percentage of cells", main="")
## lines(perc.cpd, type="b", col="red")
## lines(perc.ctv, type="b", col="blue")
## legend("topleft", c("CFSE","CPD","CTV"), pch=1, col=c("black","red","blue"), bg = 'gray90',text.col = "green4")
## axis(2, at=seq(0,50,10), labels=paste(seq(0,50,10),"%"))
## axis(1, at=1:16,labels=1:16)
## text(8,40,paste("Chi-squared Test p=", round(Xsq$p.value, digits=4), sep=""))


###################################################
### code chunk number 53: fig14
###################################################
plot(perc.cfse, type="b", axes=F, ylim=c(0,50),
  xlab="generations", ylab="Percentage of cells", main="")
lines(perc.cpd, type="b", col="red")
lines(perc.ctv, type="b", col="blue")
legend("topleft", c("CFSE","CPD","CTV"), pch=1, col=c("black","red","blue"), bg = 'gray90',text.col = "green4")
axis(2, at=seq(0,50,10), labels=paste(seq(0,50,10),"%"))
axis(1, at=1:16,labels=1:16)
text(8,40,paste("Chi-squared Test p=", round(Xsq$p.value, digits=4), sep=""))


