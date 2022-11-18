## ==========================================================================
## Accessor to data slot. This returns the flow data object (flowFrame )
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

##########################################
# parentFittingData: return flowFrame
##########################################
setMethod(
  f="Data",
  signature=c("parentFittingData"),
  def=function(object) {
    return(object@data)
  }
)

#############################################
# proliferationFittingData: return flowFrame
#############################################
setMethod(
  f="Data",
  signature=c("proliferationFittingData"),
  def=function(object) {
    return(object@data)
  }
)
