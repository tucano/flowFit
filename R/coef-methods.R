## ---------------------------------------------------------------------------
## coef is a generic function which extracts model coefficients from objects returned by
## modeling functions. coefficients is an alias for it.
## ---------------------------------------------------------------------------

###############################
# parentFittingData: coef
###############################
setMethod(
  f="coef",
  signature="parentFittingData",
  def=function(object)
  {
    return(coef(object@lmOutput))
  }
)

#################################
# proliferationFittingData: coef
#################################
setMethod(
  f="coef",
  signature="proliferationFittingData",
  def=function(object)
  {
    return(coef(object@lmOutput))
  }
)