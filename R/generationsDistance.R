## =========================================================================##
## =========================================================================##
##                          generationsDistance                             ##
## =========================================================================##
## =========================================================================##
generationsDistance <- function(dataRange,  logDecades)
{
  # Relative fluorescent units for the t0 peak (at position maxRange)
  .pmean = dataRange
  rfi <- 10^((.pmean * logDecades) / dataRange)
  # next peaks is at rfi/2
  rfi.next <- rfi/2
  # location of the next peak:
  firstgen.mean <- (dataRange*log(rfi.next))/(logDecades*log(10))
  # distance between generations
  return(.pmean - firstgen.mean)
}