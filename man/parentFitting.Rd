\name{parentFitting}
\alias{parentFitting}
\alias{parentFitting-methods}
\alias{parentFitting,flowFrame-method}

\title{ Fitting a parent population }

\description{ Estimate the proliferation of a cell population in cell-tracking dye studies. parentFitting: fit the parent population }

\usage{
parentFitting(flowframe, channel, 
              estimatedPeakPosition = NA, 
              estimatedPeakSize = NA, 
              dataRange = NA,  
              logDecades = NA, 
              binning = TRUE, 
              breaks = 1024, 
              dataSmooth = TRUE, 
              smoothWindow = 2, 
              fixedModel = FALSE, 
              fixedPars = NA, 
              verbose = FALSE )
}

\arguments{
  \item{flowframe}{An object of class \code{\link{flowFrame}} from \code{\link{flowCore}} }
  \item{channel}{FACS column/channel (\code{\link{flowFrame}} column) }
  \item{estimatedPeakPosition}{Estimated peak position. If not provided will be used the \code{\link{exprs}} mean}
  \item{estimatedPeakSize}{Estimated peak size. If not provided will be used the \code{\link{exprs}} standard deviation }
  \item{dataRange}{Number of digital data points on the machine. If not provided will be extracted from \code{\link{flowFrame}} using \code{\link{keyword}} }
  \item{logDecades}{FACS dynamic range (log decades). If not provided will be extracted from \code{\link{flowFrame}} using \code{\link{keyword}}}
  \item{binning}{Should I bin data? Some FACS have a large data range (Es: FACSCanto have 65536 data points, may be is convenient in this case to group data in bins to avoid acquiring too many cells). If you have you data log tranformed in range 0-5 it is mandatory to bin data}
  \item{breaks}{How many breaks if I bin data?}
  \item{dataSmooth}{ Should I smooth data with a Kolmogorov-Zurbenko low-pass linear filter? }
  \item{smoothWindow}{ Window used to smooth data with the Kolmogorov-Zurbenko low-pass linear filter. }
  \item{fixedModel}{ Should I use a model with fixed parameters? (Peak Position or Size). }
  \item{fixedPars}{ A list of fixed parameters. If you give me a value, I use that value, otherwise I use estimates (check examples) }
  \item{verbose}{ Verbose mode. }
}

\value{
return a \code{\link[flowFit:parentFittingData-class]{parentFittingData}} object
}

\details{
  The formula used to fit the parent population:

  \deqn{ a^2\exp\frac{(x - \mu)^2}{2\sigma^2} }

  The algorithm estimate the position (\eqn{\mu}) and size (\eqn{\sigma}) of the \emph{Parent Population}.
}

\author{
Davide Rambaldi
}

\seealso{
  \code{\link{proliferationFitting}}
}

\examples{

  if(require(flowFitExampleData)){
    # CFSE
    data(QuahAndParish)
    parent.fitting.cfse <- parentFitting(QuahAndParish[[1]], "<FITC-A>")

    parent.fitting.cfse
    summary(parent.fitting.cfse)
    confint(parent.fitting.cfse)
    coef(parent.fitting.cfse)
    plot(parent.fitting.cfse)
    Data(parent.fitting.cfse)

    # PKH26
    data(PKH26data)
    parent.fitting <-  parentFitting(PKH26data[[1]], "FL2-Height LOG")
    parent.fitting
    summary(parent.fitting)
    confint(parent.fitting)
    coef(parent.fitting)
    plot(parent.fitting)
    Data(parent.fitting)

    # fixedModel with estimates
    parent.fitting <-  parentFitting(PKH26data[[1]], "FL2-Height LOG", 
                                     fixedModel=TRUE, fixedPars=list(M=NA,S=NA))

    # fixedModel with user values
    parent.fitting <-  parentFitting(PKH26data[[1]], "FL2-Height LOG", 
                                     fixedModel=TRUE, fixedPars=list(M=810,S=16))

    # fixedModel with locked Peak Size (one fixed parameter)
    parent.fitting <-  parentFitting(PKH26data[[1]], "FL2-Height LOG", 
                                     fixedModel=TRUE, fixedPars=list(S=17))
  }
}

\references{
  \enumerate{
    \item Timur V. Elzhov, Katharine M. Mullen and Ben Bolker (2012). \bold{minpack.lm: R interface to the Levenberg-Marquardt nonlinear least-squares algorithm found in MINPACK.}
  }
}