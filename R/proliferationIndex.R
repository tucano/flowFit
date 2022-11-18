## =========================================================================##
## =========================================================================##
##                          proliferationIndex                              ##
## =========================================================================##
## =========================================================================##
proliferationIndex <- function(object)
{
  tot.events <- sum(object@modelPoints$y)
  # absolute values of cells is %/100*tot.events
  my.count <- round((object@generations/100) * tot.events)
  temp <- vector()
  for(i in 1:length(object@generations)) {
    temp[i] <- (my.count[i] / 2^(i-1))
  }
  prolif.index <- (sum(my.count)/sum(temp))
  return(prolif.index)
}