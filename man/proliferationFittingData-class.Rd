\name{proliferationFittingData-class}
\Rdversion{1.0}
\docType{class}

\alias{proliferationFittingData}
\alias{proliferationFittingData-class}
\alias{initialize,proliferationFittingData-method}
\alias{show,proliferationFittingData-method}
\alias{summary,proliferationFittingData-method}
\alias{Data,proliferationFittingData-method}
\alias{coef,proliferationFittingData-method}
\alias{confint,proliferationFittingData-method}

\title{Class "proliferationFittingData"}

\description{
Provides S4 data structure and basic infrastructure and functions to store proliferation tracking data of the Parent Population.
}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("proliferationFittingData", flowframe, channel, ...)}.
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
    \item{\code{estimatedDistance}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{parentPeakPosition}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{parentPeakSize}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{fittingDeviance}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{fixedModel}:}{Object of class \code{"logical"} ~~ }
    \item{\code{fixedPars}:}{Object of class \code{"list"} ~~ }
    \item{\code{generationsDistance}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{heights}:}{Object of class \code{"list"} ~~ }
    \item{\code{generations}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{binning}:}{Object of class \code{"logical"} ~~ }
    \item{\code{breaks}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{dataSmooth}:}{Object of class \code{"logical"} ~~ }
    \item{\code{smoothWindow}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{parStart}:}{Object of class \code{"list"} ~~ }
    \item{\code{dataMatrix}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{dataPoints}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{modelPoints}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{model}:}{Object of class \code{"function"} ~~ }
    \item{\code{residFun}:}{Object of class \code{"function"} ~~ }
    \item{\code{lmOutput}:}{Object of class \code{"nls.lm"} ~~ }
    \item{\code{numberOfPeaks}:}{Object of class \code{"numeric"} ~~ }
  }
}

\section{Methods}{
  \describe{
      \item{\link{plot}}{
      Basic plots for \code{proliferationFittingData} objects.
      \emph{Usage:}
      \code{plot(proliferationFittingData, main="Original data and Parent fitting", xlab="FACS data range", ylab="# of Events", showLegend=TRUE, logScale = TRUE, drawGrid = TRUE, ...)}
    }

    \item{\link{show}}{
       Display details about the \code{proliferationFittingData} object.
    }

    \item{\link{summary}}{
      Return a descriptive summary about the \code{proliferationFittingData} object.
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
  \code{\link{proliferationFitting}}
}

\examples{
  showClass("proliferationFittingData")
  if(require(flowFitExampleData)){
    data(PKH26data)
    parent.fitting <-  parentFitting(PKH26data[[1]], "FL2-Height LOG")
    my.fit <- proliferationFitting(PKH26data[[2]], 
                                   "FL2-Height LOG", 
                                   parent.fitting@parentPeakPosition, 
                                   parent.fitting@parentPeakSize)
    my.fit
    summary(my.fit)
    confint(my.fit)
    coef(my.fit)
    plot(my.fit)
    Data(my.fit)
  }
}

\keyword{classes}