## =========================================================================##
## =========================================================================##
##                              logTicks                                    ##
## This function draw a log scale on your FACS plots                        ##
## =========================================================================##
## =========================================================================##
logTicks <- function(dataRange, logDecades, doScale = TRUE)
{
  .scale.factor = 1
  if (doScale) {
    .scale.factor = dataRange / logDecades
  }
  .ticks = NULL
  major<-(10^(0:logDecades))
  for(i in 1:(length(major)-1)) { .ticks <- c(.ticks, seq(major[i],major[i+1],l=10)) }
  .ticks<-unique(log(.ticks,10)); #Log base 10 and remove duplicates
  major<-log(major,10)
  my.labels <- formatC(10^major, format="e", digits=0) # labels
  my.ticks <- list(major = (major*.scale.factor), all = (.ticks*.scale.factor), labels = my.labels)
  return(my.ticks)
}