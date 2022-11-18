if(require(flowFitExampleData)){

  data(QuahAndParish)
  parent.fitting.cfse <- parentFitting(QuahAndParish[[1]], "<FITC-A>")
  fitting.cfse <- proliferationFitting(QuahAndParish[[2]], "<FITC-A>", parent.fitting.cfse@parentPeakPosition, parent.fitting.cfse@parentPeakSize)
  summary(fitting.cfse)
  confint(fitting.cfse)
  coef(fitting.cfse)
  Data(fitting.cfse)
  plot(parent.fitting.cfse)
  plot(fitting.cfse)

  # for CPD sample we use a Fixed Model: we keep fixed in the model the Parent Peak Position
  parent.fitting.cpd <- parentFitting(QuahAndParish[[1]], "<APC-A>")
  fitting.cpd <- proliferationFitting(QuahAndParish[[3]], "<APC-A>", parent.fitting.cpd@parentPeakPosition, parent.fitting.cpd@parentPeakSize, fixedModel=TRUE, fixedPars=list(M=parent.fitting.cpd@parentPeakPosition))

  parent.fitting.ctv <- parentFitting(QuahAndParish[[1]], "<Alexa Fluor 405-A>")
  fitting.ctv <- proliferationFitting(QuahAndParish[[4]], "<Alexa Fluor 405-A>", parent.fitting.ctv@parentPeakPosition, parent.fitting.ctv@parentPeakSize)

  # let's compare the generations across the 3 samples:
  plot(parent.fitting.cfse, main="CFSE Non Stimulated")
  plot(fitting.cfse, which=3, main="CFSE")
  plot(fitting.cfse, which=4, main="CFSE")
  plot(fitting.cfse, which=5, main="CFSE")
  plot(parent.fitting.cpd, main="CPD Non Stimulated")
  plot(fitting.cpd, which=3, main="CPD")
  plot(fitting.cpd, which=4, main="CPD")
  plot(fitting.cpd, which=5, main="CPD")
  plot(parent.fitting.ctv, main="CTV Non Stimulated")
  plot(fitting.ctv, which=3, main="CTV")
  plot(fitting.ctv, which=4, main="CTV")
  plot(fitting.ctv, which=5, main="CTV")

  # ESTIMATE GOODNESS of FITTING with KS TEST
  perc.cfse <- fitting.cfse@generations
  perc.cpd <- fitting.cpd@generations
  perc.ctv <- fitting.ctv@generations
  perc.cfse <- c(perc.cfse, rep(0,6))

  # EXPLORATIVE PLOT
  plot(perc.cfse, type="b", axes=F, ylim=c(0,50), xlab="generations", ylab="Percentage of cells", main="")
  lines(perc.cpd, type="b", col="red")
  lines(perc.ctv, type="b", col="blue")
  legend("topleft", c("CFSE","CPD","CTV"), pch=1, col=c("black","red","blue"), bg = 'gray90',text.col = "green4")
  axis(2, at=seq(0,50,10), labels=paste(seq(0,50,10),"%"))
  axis(1, at=1:16,labels=1:16)

  # Pearson's Chi-squared Test for Count Data
  M <- rbind(perc.cfse, perc.cpd, perc.ctv)
  colnames(M) <- 1:16
  (Xsq <- chisq.test(M, B=100000, simulate.p.value=TRUE))
  text(8,40,paste("Chi-squared Test p=", round(Xsq$p.value, digits=4), sep=""))

  # COMPARING The sum of the squared residual vector (as log).
  cfse.lab <- paste("CFSE log residuals:",round(log(fitting.cfse@fittingDeviance), digits=2))
  ctv.lab <- paste("CTV log residuals:",round(log(fitting.ctv@fittingDeviance), digits=2))
  cpd.lab <- paste("CPD log residuals:",round(log(fitting.cpd@fittingDeviance), digits=2))
  plot(1:(fitting.cfse@lmOutput$niter+1), log(fitting.cfse@lmOutput$rsstrace), type="b", main="log residual sum of squares\n vs iteration number", xlab="iteration", ylab="log residual sum of squares", pch=21, bg="black", col="black", ylim=c(8,15))
  lines(1:(fitting.cpd@lmOutput$niter+1), log(fitting.cpd@lmOutput$rsstrace), type="b", pch=21, bg="red", col="red")
  lines(1:(fitting.ctv@lmOutput$niter+1), log(fitting.ctv@lmOutput$rsstrace), type="b", pch=21, bg="blue", col="blue")
  legend("topright", c(cfse.lab, ctv.lab, cpd.lab), pch=21, col=c("black","red","blue"), bg = 'gray90',text.col = "green4")

}