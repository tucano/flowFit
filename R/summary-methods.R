## ==========================================================================
## summary methods provide useful summary statistics about an object. In
## flowCore, this will be mainly summaries about filtering opration, which
## are represented by their own class
## ==========================================================================

###############################
# parentFittingData: summary
###############################
setMethod(
  f="summary",
  signature="parentFittingData",
  def=function(object, ...)
  {
    cat("********* flowFit: Parent Population Data Object ********* \n")
    cat("* Data Range = "); cat(object@dataRange); cat("\n")
    cat("* Log Decades = "); cat(object@logDecades); cat("\n")
    cat("* Channel used for fitting = "); cat(object@channel); cat("\n")
    cat("* Number of events in flowFrame = "); cat(as.numeric(dim(object@data)[1])); cat("\n")

    if (length(object@lmOutput) > 0) {
      cat("******************** Fitting ********************** \n")
      cat(.getModelType(object));cat("\n")
      cat("* Parent Peak Position = "); cat(object@parentPeakPosition); cat("\n")
      cat("* Parent Peak Size = "); cat(object@parentPeakSize); cat("\n")
      cat("* Fitting Deviance = "); cat(object@fittingDeviance); cat("\n")
      cat("* Termination message:\n* "); cat(object@lmOutput$message); cat("\n")
      cat("* Iterations: "); cat(object@lmOutput$niter); cat("\n")
      cat("* Residual degrees-of-freedom: "); cat(df.residual(object@lmOutput)); cat("\n")
    }
    cat("********************************************************** \n")
  }
)

####################################
# proliferationFittingData: summary
####################################
setMethod(
  f="summary",
  signature="proliferationFittingData",
  def=function(object,...)
  {
    cat("********* flowFit: Final Population Data Object ********* \n")
    cat("* Data Range = "); cat(object@dataRange); cat("\n")
    cat("* Log Decades = "); cat(object@logDecades); cat("\n")
    cat("* Channel used for fitting = "); cat(object@channel); cat("\n")
    cat("* Number of events in flowFrame = "); cat(as.numeric(dim(object@data)[1])); cat("\n")
    if (length(object@lmOutput) > 0)
    {
      cat("******************** Fitting ********************** \n")
      cat(.getModelType(object));cat("\n")
      cat("* Number of Peaks in the model = "); cat(object@numberOfPeaks); cat("\n")
      cat("* Parent Peak Position = "); cat(object@parentPeakPosition); cat("\n")
      cat("* Parent Peak Size = "); cat(object@parentPeakSize); cat("\n")
      cat("* Estimated generations distance:"); cat(object@generationsDistance); cat("\n")
      cat("* Fitting Deviance = "); cat(object@fittingDeviance); cat("\n")
      cat("******************** Model ********************** \n")
      cat("* Termination message:\n * "); cat(object@lmOutput$message); cat("\n")
      cat("* Iterations: "); cat(object@lmOutput$niter); cat("\n")
      cat("* Residual degrees-of-freedom: "); cat(df.residual(object@lmOutput)); cat("\n")
      cat("******************** Parameters ******************** \n")
      for (i in 1:length(object@heights)) { cat(paste(names(object@heights)[i], " = ", object@heights[[i]]," \n",sep="")) }
      cat("*********************** Generations ********************** \n")
      for (i in 1:length(object@generations)) { cat(paste("* generation ", i, " = ", round(object@generations[i],digits=2)," %\n",sep="")) }
    }
    cat("********************************************************** \n")
  }
)