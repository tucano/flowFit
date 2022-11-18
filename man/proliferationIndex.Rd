\name{proliferationIndex}
\alias{proliferationIndex}
\alias{proliferationIndex-methods}
\alias{proliferationIndex,proliferationFittingData-method}

\title{ proliferation index calculator }
\description{
Proliferation index calculator. Proliferation index is calculated as the sum of the cells in all generations including the parental divided by the computed number of original parent cells theoretically present at the start of the experiment. It is a measure of the fold increase in cell number in the culture over the course of the experiment.
}

\usage{
proliferationIndex(object)
}

\details{
  The formula is:
  \eqn{\frac{\sum_0^iN_i}{\sum_0^iN^i/2^i}}
  Where \eqn{i} is the generation number (parent generation = 0).
  In the absence of proliferation, that is, when all cells are in the parent generation, the formula gives:
  \eqn{\frac{N_0}{N_0/2^0} = 1}
  defining the lower limit of the PI.
}

\arguments{
  \item{object}{
    An object of class \code{\link{proliferationFittingData}}
  }
}

\value{
  return a \code{numeric}
}

\references{
  \enumerate{
    \item{ Munson ME. An improved technique for calculating relative response in cellular proliferation experiments. Cytometry A. 2010 Oct;77(10):909-10. doi: 10.1002/cyto.a.20935. Erratum in: Cytometry A. 2010 Dec;77(12):1177. PubMed PMID: 21290464.}
  }
}

\author{ Davide Rambaldi }

\seealso{
  \code{\link{proliferationFitting}}
  \code{\link{proliferationFittingData-class}}
}

\examples{
  # load data
  if(require(flowFitExampleData)){
    data(PKH26data)
    parent.fitting <-  parentFitting(PKH26data[[1]], "FL2-Height LOG")
    my.fit <- proliferationFitting(PKH26data[[2]], "FL2-Height LOG", 
                                   parent.fitting@parentPeakPosition, 
                                   parent.fitting@parentPeakSize)
    my.index <- proliferationIndex(my.fit)
  }
}