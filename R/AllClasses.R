## =========================================================================##
## =========================================================================##
##                    Class definitions, contructors and helpers            ##
## =========================================================================##
## =========================================================================##



## ===========================================================================
##  parentFittingData
## ---------------------------------------------------------------------------
## An object of class parentFittingData Provides S4 data structure and basic
## infrastructure and functions to store proliferation tracking data of the
## parent population
## ---------------------------------------------------------------------------
## Need this to be able to use the S3 class as a slot
setOldClass("nls.lm")

setClass (
  Class="parentFittingData",
  representation=representation
  (
    # DATA SLOTS
    data    = "flowFrame",              # store a flowFrame from flowCore
    channel = "character",              # FACS Stain (flowFrame column)
    estimatedPeakPosition = "numeric",  # ESTIMATE PEAK POSITION
    estimatedPeakSize = "numeric",      # ESTIMATE PEAK SIZE
    # FACS SLOTS
    dataRange  = "numeric",             # number of data points on the machine
    logDecades = "numeric",             # log decades dynamic range
    # FITTING
    parentPeakPosition = "numeric",     # the peak position after fitting
    parentPeakSize     = "numeric",     # the peak size after fitting
    fittingDeviance    = "numeric",     # the fitting deviance
    # PARAMS
    binning  = "logical",               # should I bin data?
    breaks   = "numeric",               # If I bin: how many breaks?
    dataSmooth = "logical",             # should I smooth data?
    smoothWindow = "numeric",           # Smoothing window
    fixedModel = "logical",             # should I keep parent mean, sd and distance fixed?
    fixedPars = "list",                 # Fixed MEAN, SD (for fixed model)
    # INTERNAL
    parStart = "list",                  # LM algorithm parameters start (guess)
    lmOutput = "nls.lm",                # lm output as list
    dataMatrix = "matrix",              # THE ORIGINAL DATA as a MATRIX
    dataPoints = "data.frame",          # THE ORIGINAL exprs DATA Tabulated (or binned)
    modelPoints = "data.frame",         # THE exprs DATA AFTER TRANSFORMATIONS
    residFun = "function"               # The function used in modeling
  ),
  prototype=prototype
  (
    estimatedPeakPosition = NA_real_,
    estimatedPeakSize = NA_real_,
    dataRange = NA_real_,
    logDecades = NA_real_,
    parentPeakPosition = NA_real_,
    parentPeakSize = NA_real_,
    fittingDeviance = NA_real_,
    binning = NA,
    breaks = NA_real_,
    dataSmooth = NA,
    smoothWindow = NA_real_,
    parStart = list(),
    dataMatrix = matrix(),
    dataPoints = data.frame(),
    modelPoints = data.frame(),
    lmOutput =structure(list(), class="nls.lm")
  )
)

setMethod (
  f = "initialize",
  signature = "parentFittingData",
  definition = function(.Object, flowframe, channel, ... )
  {
    # check flowFrame data
    .checkClass(flowframe, c("flowFrame"))
    .Object@data <- flowframe
    if(!(channel %in% colnames(.Object@data)))
    {
      stop("Invalid channel name not matching the flowFrame:\n    ", paste(channel, collapse=", "))
    }
    .Object@channel = channel
    # FIRST call basic initialize
    .Object <- callNextMethod(.Object, ...)
    return(.Object)
  }
)


## ===========================================================================
##  proliferationFittingData
## ---------------------------------------------------------------------------
## This class represents an S4 data structure with basic infrastructure
## and functions to store proliferation tracking data.
## ---------------------------------------------------------------------------
setClass(
  Class="proliferationFittingData",
  representation=representation(
    # DATA SLOTS
    data    = "flowFrame",              # store a flowFrame from flowCore
    channel = "character",              # FACS Stain (flowFrame column)
    estimatedPeakPosition = "numeric",  # ESTIMATE PEAK POSITION
    estimatedPeakSize = "numeric",      # ESTIMATE PEAK SIZE
    # FACS SLOTS
    dataRange  = "numeric",             # number of data points on the machine
    logDecades = "numeric",             # log decades dynamic range
    estimatedDistance = "numeric",      # estimated from formula (before fitting)
    # PARAMS
    binning = "logical",                # should I bin data?
    breaks = "numeric",                 # If I bin: how many breaks?
    dataSmooth = "logical",             # should I smooth data?
    smoothWindow = "numeric",           # Window used to smooth data
    fixedModel = "logical",             # should I keep parent mean, sd and distance fixed?
    fixedPars = "list",                 # Fixed MEAN, SD and DISTANCE (for fixed model)
    # INTERNAL
    parStart = "list",                  # LM algorithm parameters start (guess)
    lmOutput = "nls.lm",                # lm output as list
    dataMatrix = "matrix",              # THE ORIGINAL DATA as a MATRIX
    dataPoints = "data.frame",          # THE ORIGINAL exprs DATA Tabulated (or binned)
    modelPoints = "data.frame",         # THE exprs DATA AFTER TRANSFORMATIONS
    residFun = "function",              # The residual function used in modeling
    model    = "function",              # The function used in modeling
    numberOfPeaks = "numeric",          # number of PEaks in the model
    # FITTING
    parentPeakPosition  = "numeric",    # the peak position after fitting
    parentPeakSize      = "numeric",    # the peak size after fitting
    fittingDeviance     = "numeric",    # the fitting deviance
    generationsDistance = "numeric",    # final generationsDistance (after fitting)
    heights = "list",                   # list of peak heights
    generations = "numeric"             # vector of percentage of cells/generation
  ),
  prototype=prototype(
    estimatedPeakPosition = NA_real_,
    estimatedPeakSize = NA_real_,
    dataRange = NA_real_,
    logDecades = NA_real_,
    estimatedDistance = NA_real_,
    parentPeakPosition = NA_real_,
    parentPeakSize = NA_real_,
    fittingDeviance = NA_real_,
    binning = NA,
    breaks = NA_real_,
    dataSmooth = NA,
    smoothWindow = NA_real_,
    fixedModel = NA,
    fixedPars = list(),
    parStart = list(),
    dataMatrix = matrix(),
    dataPoints = data.frame(),
    modelPoints = data.frame(),
    numberOfPeaks = NA_real_,
    lmOutput =structure(list(), class="nls.lm"),
    generationsDistance = NA_real_,
    heights = list(),
    generations = NA_real_
  )
)

setMethod (
  f = "initialize",
  signature = "proliferationFittingData",
  definition = function(.Object, flowframe, channel, ... )
  {
    # check flowFrame data
    .checkClass(flowframe, c("flowFrame"))
    .Object@data <- flowframe
    if(!(channel %in% colnames(.Object@data)))
    {
      stop("Invalid channel name not matching the flowFrame:\n    ", paste(channel, collapse=", "))
    }
    .Object@channel = channel
    # FIRST call basic initialize
    .Object <- callNextMethod(.Object, ...)
    return(.Object)
  }
)
