\name{parentFittingData-class}
\Rdversion{1.0}
\docType{class}

\alias{parentFittingData}
\alias{parentFittingData-class}
\alias{initialize,parentFittingData-method}
\alias{show,parentFittingData-method}
\alias{summary,parentFittingData-method}
\alias{Data,parentFittingData-method}
\alias{coef,parentFittingData-method}
\alias{confint,parentFittingData-method}

\title{Class "parentFittingData"}

\description{
Provides S4 data structure and basic infrastructure and functions to store proliferation tracking data of the Parent Population.
}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("parentFittingData", flowframe, channel, ...)}.
This class is for internal use.
}
\section{Slots}{

  \describe{
    \item{\code{data}:}{Object of class \code{"flowFrame"} ~~ }
    \item{\code{channel}:}{Object of class \code{"character"} ~~ }
    \item{\code{estimatedPeakPosition}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{estimatedPeakSize}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{dataRange}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{logDecades}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{parentPeakPosition}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{parentPeakSize}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{fittingDeviance}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{fixedModel}:}{Object of class \code{"logical"} ~~ }
    \item{\code{fixedPars}:}{Object of class \code{"list"} ~~ }
    \item{\code{binning}:}{Object of class \code{"logical"} ~~ }
    \item{\code{breaks}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{dataSmooth}:}{Object of class \code{"logical"} ~~ }
    \item{\code{smoothWindow}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{parStart}:}{Object of class \code{"list"} ~~ }
    \item{\code{dataMatrix}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{dataPoints}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{modelPoints}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{residFun}:}{Object of class \code{"function"} ~~ }
    \item{\code{lmOutput}:}{Object of class \code{"nls.lm"} ~~ }
  }
}

\section{Methods}{
  \describe{
      \item{\link{plot}}{
      Basic plots for \code{parentFittingData} objects.
      \emph{Usage:}
      \code{plot(parentFittingData, main="Original data and Parent fitting", xlab="FACS data range", ylab="# of Events", showLegend=TRUE, logScale = TRUE, drawGrid = TRUE, ...)}
    }

    \item{\link{show}}{
       Display details about the \code{parentFittingData} object.
    }

    \item{\link{summary}}{
      Return a descriptive summary about the \code{parentFittingData} object.
    }

    \item{\link{Data}}{
      Return the \code{\link{flowFrame}} object.
    }

    \item{\link{coef}}{
      Return coefficients of the model.
    }

    \item{\link{confint}}{
      Return confidence intervals of the model.
    }

   }
}

\author{ Davide Rambaldi }

\seealso{
\code{\link{parentFitting}}
}

\examples{
  showClass("parentFittingData")
  if(require(flowFitExampleData)){
    data(PKH26data)
    parent.fitting <-  parentFitting(PKH26data[[1]], "FL2-Height LOG")
    parent.fitting
    summary(parent.fitting)
    confint(parent.fitting)
    coef(parent.fitting)
    plot(parent.fitting)
    Data(parent.fitting)
  }
}

\keyword{classes}