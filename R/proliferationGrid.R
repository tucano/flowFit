## =========================================================================##
## =========================================================================##
##                          proliferationGrid                               ##
## =========================================================================##
## =========================================================================##
proliferationGrid <- function(parentPosition,
  fittedDistance = NA, dataRange = 1024, logDecades = 4,
  lwd=1, lty=3, col=rgb(0,0,0,0.5))
{
  .dist <- ifelse(is.na(fittedDistance),generationsDistance(dataRange, logDecades), fittedDistance)
  if (dev.cur() == 1)
    stop("Plot has not be called yet.")
  .this.step = parentPosition
  .steps = round(parentPosition / .dist)
  for (i in 1:.steps)
  {
    abline(v = .this.step, lwd=lwd, lty=lty, col=col)
    .this.step = .this.step - .dist
  }
}