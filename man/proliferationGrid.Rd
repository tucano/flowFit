\name{proliferationGrid}
\alias{proliferationGrid}
\alias{proliferationGrid-methods}
\alias{proliferationGrid,ANY-method}

\title{ proliferationGrid function for plotting }

\description{
  This function draw a proliferation grid. The grid marks the distance between cell generations calculated with the function \code{\link{generationsDistance}}
}

\usage{
  proliferationGrid(parentPosition,
    fittedDistance = NA, dataRange = 1024, logDecades = 4,
    lwd=1, lty=3, col=rgb(0,0,0,0.5))
}

\arguments{
  \item{parentPosition}{ Position of the parent Peak from \code{\link{parentFitting}}}
  \item{fittedDistance}{ You can provide the distance estimated from the \code{\link{proliferationFitting}} function}
  \item{dataRange}{ Range of your data (number of data points in the FACS)}
  \item{logDecades}{ Number of log decades in the FACS}
  \item{lwd}{ Grid line size. See \code{\link{par}} }
  \item{lty}{ Grid line type. See \code{\link{par}} }
  \item{col}{ Grid color. See \code{\link{par}} and \code{\link{rgb}} }
}

\author{ Davide Rambaldi }

\examples{
  plot(c(0,1023),c(0,1000), 
       xlim=c(0,1023),
       ylim=c(0,1000),  
       xlab="FACS CHANNEL", 
       ylab="# OF EVENTS", 
       main="A flowFit Empty Plot")
       
  # create a grid with parent at 800
  proliferationGrid(1000, dataRange=1024, logDecades=4)
}