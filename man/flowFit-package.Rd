\name{flowFit-package}
\alias{flowFit-package}
\alias{flowFit}
\docType{package}

\title{ Estimate proliferation in cell-tracking dye studies }

\description{
  This package estimate the proliferation of a cell population in cell-tracking dye studies.

  In cells proliferation tracking experiments, cells are stained with a tracking dye before culture.
  During cell division, the tracking dye is partitioned between daughter cells, so that each division brings about a halving of fluorescence intensity; the intensity of a cell, by comparison with the intensity of resting cells, provides an indication of how many divisions the cell has undergone since stimulation

  This package uses an R implementation of the Levenberg-Marquardt algorithm (\code{\link{nls.lm}}) to fit a set of peaks (corresponding to different generations of cells) over the proliferation-tracking dye distribution in a FACS experiment.

  The package define two data structure (S4 classes): \code{\link[flowFit:proliferationFittingData-class]{proliferationFittingData}}, \code{\link[flowFit:parentFittingData-class]{parentFittingData}} and their methods and accessors.

  The package is integrated with other \url{www.bioconductor.org} libraries for analysis of flow cytometry data: \code{\link{flowCore}} and \code{\link{flowViz}}.
}

\details{
  \tabular{ll}{
    Package: \tab flowFit\cr
    Type: \tab Package\cr
    Version: \tab 0.2\cr
    Date: \tab 2012-11-29\cr
    License: \tab Artistic-2.0\cr
  }
}


\author{
  Maintainer: Davide Rambaldi <davide.rambaldi@gmail.com>
  Author: Davide Rambaldi
}

\references{
  \enumerate{
    \item Timur V. Elzhov, Katharine M. Mullen and Ben Bolker (2012). \bold{minpack.lm: R interface to the Levenberg-Marquardt nonlinear least-squares algorithm found in MINPACK.}
    \item \bold{Tracking antigen-driven responses by flow cytometry: Monitoring proliferation by dye dilution.} Paul K Wallace, Joseph D Tario, Jan L Fisher, Stephen S Wallace, Marc S Ernstoff, Katharine A Muirhead. Cytometry (2008) vol. 73A (11) pp. 1019-1034
    \item \bold{MEASURING MOLECULES OF EQUIVALENT FLUOROCHROME (MEF) USING SPHEROTM RAINBOW AND ULTRA RAINBOW CALIBRATION PARTICLES}, Spherotech, http://www.spherotech.com/tech_SpheroTech_Note_9.html
    \item Benjamin J.C. Quah, Christopher R. Parish, \bold{New and improved methods for measuring lymphocyte proliferation in vitro and in vivo using CFSE-like fluorescent dyes}, Journal of Immunological Methods, Volume 379, Issues 1-2, 31 May 2012, Pages 1-14, ISSN 0022-1759, 10.1016/j.jim.2012.02.012.
  }
}

\keyword{ package }

\seealso{
  \enumerate{
    \item \code{\link{proliferationFitting}} generations fitting function.
    \item \code{\link{parentFitting}} parent population fitting function.
    \item \code{\link{proliferationIndex}} proliferation index calculator.
    \item \code{\link{getGenerations}} get percentage of cells for generation.
    \item \code{\link{logTicks}} draw a log scale on your FACS plots.
    \item \code{\link{generationsDistance}} calculate the distance between 2 generations of cells on the FACS scale.
  }
}

\examples{
  if(require(flowFitExampleData)){
    data(QuahAndParish)
    parent.fitting.cfse <- parentFitting(QuahAndParish[[1]], "<FITC-A>")
    fitting.cfse <- proliferationFitting(QuahAndParish[[2]], "<FITC-A>", 
                                         parent.fitting.cfse@parentPeakPosition, 
                                         parent.fitting.cfse@parentPeakSize)

    summary(fitting.cfse)
    confint(fitting.cfse)
    coef(fitting.cfse)
    Data(fitting.cfse)

    plot(parent.fitting.cfse)
    plot(fitting.cfse)

    # for this sample we use a Fixed Model: we keep fixed in the model the Parent Peak Position
    parent.fitting.cpd <- parentFitting(QuahAndParish[[1]], "<APC-A>")
    fitting.cpd <- proliferationFitting(QuahAndParish[[3]], "<APC-A>", 
                                        parent.fitting.cpd@parentPeakPosition, 
                                        parent.fitting.cpd@parentPeakSize, 
                                        fixedModel=TRUE, 
                                        fixedPars=list(M=parent.fitting.cpd@parentPeakPosition))

    parent.fitting.ctv <- parentFitting(QuahAndParish[[1]], "<Alexa Fluor 405-A>")
    fitting.ctv <- proliferationFitting(QuahAndParish[[4]], "<Alexa Fluor 405-A>", 
                                        parent.fitting.ctv@parentPeakPosition, 
                                        parent.fitting.ctv@parentPeakSize)

    # let's compare the generations across the 3 samples:
    par(mfrow=c(3,4))
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
    par(mfrow=c(1,1), ask=FALSE)
    plot(perc.cfse, type="b", axes=FALSE, ylim=c(0,50), 
         xlab="generations", ylab="Percentage of cells", main="")
    lines(perc.cpd, type="b", col="red")
    lines(perc.ctv, type="b", col="blue")
    legend("topleft", c("CFSE","CPD","CTV"), pch=1, 
         col=c("black","red","blue"), bg = 'gray90',text.col = "green4")
    axis(2, at=seq(0,50,10), labels=paste(seq(0,50,10),"\%"))
    axis(1, at=1:16,labels=1:16)

    # Pearson's Chi-squared Test for Count Data
    M <- rbind(perc.cfse, perc.cpd, perc.ctv)
    colnames(M) <- 1:16
    (Xsq <- chisq.test(M, B=100000, simulate.p.value=TRUE))
    text(8,40,paste("Chi-squared Test p=", round(Xsq$p.value, digits=4), sep=""))

    # PKH26
    # load data
    data(PKH26data)
    parent.fitting <-  parentFitting(PKH26data[[1]], "FL2-Height LOG")
    my.fit <- proliferationFitting(PKH26data[[2]], "FL2-Height LOG", 
                                   parent.fitting@parentPeakPosition, 
                                   parent.fitting@parentPeakSize)
    my.fit
    summary(my.fit)
    confint(my.fit)
    coef(my.fit)
    Data(my.fit)
    # plot results
    plot(my.fit)

    # modeling with locked Peak Size
    my.fitb <- proliferationFitting(PKH26data[[2]], "FL2-Height LOG", 
                                  parent.fitting@parentPeakPosition, 
                                  parent.fitting@parentPeakSize, 
                                  fixedModel=TRUE, 
                                  fixedPars=list(S=16))
    plot(my.fitb)
    # modeling with locked Peak Size and Position
    my.fitc <- proliferationFitting(PKH26data[[2]], 
                                  "FL2-Height LOG", 
                                  parent.fitting@parentPeakPosition, 
                                  parent.fitting@parentPeakSize, 
                                  fixedModel=TRUE, 
                                  fixedPars=list(S=16, M=810))
    plot(my.fitc)
    # modeling with locked Peak Size, Position and Distance
    my.fitd <- proliferationFitting(PKH26data[[2]], "FL2-Height LOG", 
                                  parent.fitting@parentPeakPosition, 
                                  parent.fitting@parentPeakSize, 
                                  fixedModel=TRUE, 
                                  fixedPars=list(S=16, M=810, D=76))
    plot(my.fitd)
  }
}