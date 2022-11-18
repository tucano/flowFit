\name{logTicks}
\alias{logTicks}
\alias{logTicks-methods}
\alias{logTicks,ANY-method}

\title{ Generate logTicks for FACS plotting }
\description{
  This function return a log scale to be used on your FACS plots
}

\usage{
  logTicks(dataRange, logDecades, doScale = TRUE)
}

\arguments{
  \item{dataRange}{ Range of your data (number of data points in the FACS)}
  \item{logDecades}{ Number of log decades in the FACS}
  \item{doScale}{ Scale according to dataRange and logDecades: scale.factor = dataRange / logDecades }
}

\value{
  Return a \code{list} with 3 elements:
  \item{major}{Position of the Major ticks}
  \item{all}{Position of the All ticks}
  \item{label}{Labels for Major Ticks}
}

\author{ Davide Rambaldi }

\examples{

if(require(flowFitExampleData)){
  # using flowViz

  # load data
  data(PKH26data)

  # plot data
  plot(PKH26data[[1]], "FL2-Height LOG", axes=FALSE, breaks=1024)

  # create ticks
  my.ticks <- logTicks(1024,4)

  # plot your ticks
  axis(1,my.ticks$all,label=FALSE)
  axis(1,at=my.ticks$major,labels=my.ticks$labels)
  axis(2)
}
}