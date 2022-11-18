\name{proliferationFitting}
\alias{proliferationFitting}
\alias{proliferationFitting-methods}
\alias{proliferationFitting,flowFrame-method}

\title{ Estimate proliferation in cell-tracking dye studies }

\description{
  The algorithm fit a set of N peaks on the \code{flowFrame} data using the \code{\link{nls.lm}} function.
  The number of peaks to be fitted is automatically estimated using \code{\link{generationsDistance}}.

  The algorithm take the position (\eqn{\mu}) and size (\eqn{\sigma}) of the \emph{Parent Population} as estimates and fit a set of peaks on a \code{flowFrame} data.

  The first peak correspond to the parent population:

  \deqn{ a^2\exp\frac{(x - \mu)^2}{2\sigma^2} }

  The next peak (corresponding to the next generation of cells) will be:

  \deqn{ b^2\exp\frac{(x - (\mu-D))^2}{2\sigma^2} }

  Where D is the estimated distance between 2 generations of cells.

  The complete formula for the fitting of the 15 peaks is the following:

  \deqn{ a^2\exp\frac{(x - M)^2}{2s^2} + b^2\exp\frac{(x - (M-D))^2}{2s^2} + ... + p^2\exp\frac{(x - (M-14 \cdot D))^2}{2s^2} }

  Where the parameters [a-q] represent an estimate of the number of cells for a given generation.

  In the Levenberg-Marquadt algorithm implementation we use this formula to estimate the error between the model and the real data:

  \deqn{residFun = (Observed - Model)^2}

  The ration between the intergral of a single peak and the integral of all model formula is an estimate of the percentage of cells in a given generation.
}

\usage{
proliferationFitting(
  flowframe,
  channel,
  estimatedParentPosition,
  estimatedParentSize,
  dataRange = NA,
  logDecades = NA,
  estimatedDistance = NA,
  binning = TRUE,
  breaks = 1024,
  dataSmooth = TRUE,
  smoothWindow = 2,
  fixedModel = FALSE,
  fixedPars = NA,
  verbose = FALSE
  )
}

\arguments{
  \item{flowframe}{An object of class \code{\link{flowFrame}} from \code{\link{flowCore}} }
  \item{channel}{FACS column/channel (\code{\link{flowFrame}} column) }
  \item{estimatedParentPosition}{Estimated parent peak position.}
  \item{estimatedParentSize}{Estimated parent peak size.}
  \item{dataRange}{Number of digital data points on the machine. If not provided will be extracted from \code{\link{flowFrame}} using \code{\link{keyword}} }
  \item{logDecades}{FACS dynamic range (log decades). If not provided will be extracted from \code{\link{flowFrame}} using \code{\link{keyword}}}
  \item{estimatedDistance}{Estimated distance between generations. If not provided will be estimated with \code{\link{generationsDistance}}}
  \item{binning}{Should I bin data? Some FACS have a large data range (Es: FACSCanto have 65536 data points, may be is convenient in this case to group data in bins to avoid acquiring too many cells). If you have you data log tranformed in range 0-5 it is mandatory to bin data}
  \item{breaks}{How many breaks if I bin data?}
  \item{dataSmooth}{ Should I smooth data with a Kolmogorov-Zurbenko low-pass linear filter (\code{\link{kz}})? }
  \item{smoothWindow}{ Window used to smooth data with the Kolmogorov-Zurbenko low-pass linear filter (\code{\link{kz}}). }
  \item{fixedModel}{ Should I use a model with fixed parameters? (Peak Position or Size). }
  \item{fixedPars}{ A list of fixed parameters. If you give me a value, I use that value, otherwise I use estimates (check examples) }
  \item{verbose}{ Verbose mode. }
}

\details{
  See the vignette for more details on this function.
}

\value{
  return a \code{\link{proliferationFittingData}} object
}

\references{
  \enumerate{
    \item Timur V. Elzhov, Katharine M. Mullen and Ben Bolker (2012). \bold{minpack.lm: R interface to the Levenberg-Marquardt nonlinear least-squares algorithm found in MINPACK.}
  }
}


\author{ Davide Rambaldi }

\examples{

  if(require(flowFitExampleData)){
    # PKH26
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
    my.fit <- proliferationFitting(PKH26data[[2]], "FL2-Height LOG", 
                                   parent.fitting@parentPeakPosition, 
                                   parent.fitting@parentPeakSize, 
                                   fixedModel=TRUE, 
                                   fixedPars=list(S=16))
                                   
    # modeling with locked Peak Size and Position
    my.fit <- proliferationFitting(PKH26data[[2]], "FL2-Height LOG", 
                                   parent.fitting@parentPeakPosition, 
                                   parent.fitting@parentPeakSize, 
                                   fixedModel=TRUE, 
                                   fixedPars=list(S=16, M=810))
                                   
    # modeling with locked Peak Size, Position and Distance
    my.fit <- proliferationFitting(PKH26data[[2]], "FL2-Height LOG", 
                                   parent.fitting@parentPeakPosition, 
                                   parent.fitting@parentPeakSize, 
                                   fixedModel=TRUE, fixedPars=list(S=16, M=810, D=76))

    # generations as vector
    my.fit@generations
    # generations as list
    getGenerations(my.fit)

    # CFSE, CPD and CTV data
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

    # for CPD samples we use a Fixed Model: we keep fixed in the model the Parent Peak Position
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
  }
}