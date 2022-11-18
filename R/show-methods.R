## ==========================================================================
## show and print methods display details about an object, and print usually
## allows for a bit mor fine-grained control.
## ==========================================================================

###############################
# parentFittingData: show
###############################

setMethod(
  f="show",
  signature="parentFittingData",
  def=function(object)
  {
    cat("********* flowFit: Parent Population Data Object ********* \n")
    cat("* Data Range = "); cat(object@dataRange); cat("\n")
    cat("* Log Decades = "); cat(object@logDecades); cat("\n")
    cat("* Channel used for fitting = "); cat(object@channel); cat("\n")
    cat("* Number of events in flowFrame = "); cat(as.numeric(dim(object@data)[1])); cat("\n")

    if (length(object@lmOutput) > 0) {
      cat("********************************************************** \n")
      cat("* FITTING with:\n")
      cat(.getModelType(object));cat("\n")
      cat("* Parent Peak Position = "); cat(object@parentPeakPosition); cat("\n")
      cat("* Parent Peak Size = "); cat(object@parentPeakSize); cat("\n")
      cat("* Fitting Deviance = "); cat(object@fittingDeviance); cat("\n")
    }
    cat("********************************************************** \n")
  }
)

###############################
# proliferationFittingData: show
###############################
setMethod(
  f="show",
  signature="proliferationFittingData",
  def=function(object)
  {
    cat("********* flowFit: Final Population Data Object ********* \n")
    cat("* Data Range = "); cat(object@dataRange); cat("\n")
    cat("* Log Decades = "); cat(object@logDecades); cat("\n")
    cat("* Channel used for fitting = "); cat(object@channel); cat("\n")
    cat("* Number of events in flowFrame = "); cat(as.numeric(dim(object@data)[1])); cat("\n")

    if (length(object@lmOutput) > 0)
    {
      cat("********************************************************** \n")
      cat("* FITTING with:\n")
      cat(.getModelType(object));cat("\n")
      cat("* Number of Peaks in the model = "); cat(object@numberOfPeaks); cat("\n")
      cat("* Parent Peak Position = "); cat(object@parentPeakPosition); cat("\n")
      cat("* Parent Peak Size = "); cat(object@parentPeakSize); cat("\n")
      cat("* Estimated generations distance:"); cat(object@generationsDistance); cat("\n")
      cat("* Fitting Deviance = "); cat(object@fittingDeviance); cat("\n")
      cat("*********************** Generations ********************** \n")
      for (i in 1:length(object@generations)) { cat(paste("* generation ", i, " = ", round(object@generations[i],digits=2)," %\n",sep="")) }
    }

    cat("********************************************************** \n")
  }
)
