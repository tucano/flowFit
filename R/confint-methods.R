## ---------------------------------------------------------------------------
## Computes confidence intervals for one or more parameters in a fitted model.
## There is a default and a method for objects inheriting from class "lm".
## ---------------------------------------------------------------------------

###############################
# parentFittingData: confint
###############################
setMethod(
  f="confint",
  signature="parentFittingData",
  def=function(object, ...)
  {
    return(tryCatch(confint(object@lmOutput), error = function(e) paste("Can't calculate confidence intervals: ", e)))
  }
)


####################################
# proliferationFittingData: confint
####################################
setMethod(
  f="confint",
  signature="proliferationFittingData",
  def=function(object, ...)
  {
    return(tryCatch(confint(object@lmOutput), error = function(e) paste("Can't calculate confidence intervals: ", e)))
  }
)
