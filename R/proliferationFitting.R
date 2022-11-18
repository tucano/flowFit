## =========================================================================##
## =========================================================================##
##                    proliferationFitting method                           ##
## =========================================================================##
## =========================================================================##

## TODO: Fixed Model, Fixed Distance

proliferationFitting <- function(
        flowframe,
        channel,
        estimatedParentPosition,
        estimatedParentSize,
        dataRange = NA,
        logDecades = NA,
        estimatedDistance = NA,
        binning = TRUE,
        breaks = 1024,
        dataSmooth = TRUE,
        smoothWindow = 2,
        fixedModel = FALSE,
        fixedPars = NA,
        verbose = FALSE
    )
{
    # CREATE
    fittingdata <- new("proliferationFittingData", flowframe, channel,
        estimatedPeakPosition = as.numeric(estimatedParentPosition), estimatedPeakSize = as.numeric(estimatedParentSize),
        binning=binning, breaks=breaks, dataSmooth=dataSmooth, smoothWindow=smoothWindow)

    # FIXED MODEL SETUP SLOTS HERE
    fittingdata@fixedModel = fixedModel

    # set dataRange and log decades
    if (is.na(dataRange)) {
        fittingdata@dataRange = .setDataRange(flowframe, channel)
        message(paste("parentFitting: No dataRange provided. Setting range from channel range:", fittingdata@dataRange, sep=" "))
    } else {
        fittingdata@dataRange = dataRange
    }

    if (is.na(logDecades)) {
        fittingdata@logDecades = .setLogDecades(flowframe, channel)
        message(paste("proliferationFitting: No logDecades provided. Setting LOG dynamic range from flowFrame: log decades: ", fittingdata@logDecades , sep="" ))
    } else {
        fittingdata@logDecades = logDecades
    }

    # SET DISTANCE
    if (is.na(estimatedDistance)) {
        fittingdata@estimatedDistance = generationsDistance(fittingdata@dataRange, fittingdata@logDecades)
        if (verbose)
            message(paste("Estimating generations distance with formulas:", fittingdata@estimatedDistance))
    } else {
        fittingdata@estimatedDistance = estimatedDistance
    }

    # get data points
    fittingdata@dataMatrix <- exprs(fittingdata@data[,fittingdata@channel])

    # FIXME: binning AUTOMATIC for more than 1024 points?!

    # BINNING OR TABULATE: dataPoints
    if (binning)
    {
        message(paste("Binning data into: ", breaks, " breaks", sep=""))
        .h <- hist(fittingdata@dataMatrix, breaks=breaks, plot=FALSE)
        .y.points <- .h$counts
        .x.points <- .h$mids
        fittingdata@dataPoints <- data.frame( x = .x.points, y = .y.points)
    }
    else
    {
        message(paste("Tabulating data in ", fittingdata@dataRange , " FACS data points", sep=""))
        .y.points <- tabulate(fittingdata@dataMatrix, nbins=fittingdata@dataRange)
        fittingdata@dataPoints <- data.frame(x = 0:(fittingdata@dataRange - 1), y = .y.points)
    }

    # SMOOTHING: modelPoints
    if (fittingdata@dataSmooth)
    {
        message("Smoothing data with Kolmogorov-Zurbenko low-pass linear filter.")
        .y.smooth <- kz(fittingdata@dataPoints$y, fittingdata@smoothWindow)
        fittingdata@modelPoints <- data.frame(x = fittingdata@dataPoints$x, y = .y.smooth)
    }
    else
    {
        fittingdata@modelPoints <- fittingdata@dataPoints
    }

    # NUMBER OF PEAKS
    # correct number of Peaks with LAST VALUE
    .lastPoint = fittingdata@modelPoints[1,]
    .realSpace = fittingdata@estimatedPeakPosition - .lastPoint$x
    fittingdata@numberOfPeaks = ceiling(.realSpace / fittingdata@estimatedDistance)
    if (fittingdata@numberOfPeaks > 20)
    {
        fittingdata@numberOfPeaks = 20
        if (verbose)
            message("Can't find more than 20 generations. Estimated number of peaks lowered to 20")
    }

    if (verbose)
    {
        message(paste("Loading data from FCS column: ", fittingdata@channel, sep=""))
        message(paste("Setting parStart and Building model. Estimated number of peaks:", fittingdata@numberOfPeaks ))
    }

    # generate model
    .coefficients = letters[1:fittingdata@numberOfPeaks]
    fittingdata@parStart = vector(mode="list",length(.coefficients))
    names(fittingdata@parStart) <- .coefficients
    fittingdata@parStart[1:length(.coefficients)] = 1
    fittingdata@parStart$M <- fittingdata@estimatedPeakPosition
    fittingdata@parStart$S <- fittingdata@estimatedPeakSize
    fittingdata@parStart$D <- fittingdata@estimatedDistance

    # FIXED MODEL
    if (fittingdata@fixedModel) fittingdata <- .setFixedValues(fittingdata, fixedPars, verbose)

    # MODEL FORMULA
    if (fittingdata@fixedModel) {
        fittingdata@model <- .buildFixedModel(fittingdata@numberOfPeaks)
        fittingdata@residFun <- function(p,N,xx,fixedPars) { N - fittingdata@model(p,xx,fixedPars) }
    } else {
        fittingdata@model <- .buildModel(fittingdata@numberOfPeaks)
        fittingdata@residFun <- function(p,N,xx) { N - fittingdata@model(p,xx) }
    }

    # FITTING on current data
    .verbosevar <- ifelse(verbose, 1, 0)

    if (verbose) message("FITTING MODEL...")

    if (fixedModel) {
        fittingdata@lmOutput <- nls.lm(
            par=fittingdata@parStart,
            fn = fittingdata@residFun,
            N = fittingdata@modelPoints$y,
            xx = fittingdata@modelPoints$x,
            fixedPars = fittingdata@fixedPars,
            control = nls.lm.control
            (
                nprint = .verbosevar,
                maxiter = 1024,
                factor = 0.1,
                maxfev=1000000,
                ptol=.Machine$double.xmin,
                gtol= 0,
                ftol=.Machine$double.xmin
            )
        )
        fittingdata@parentPeakPosition  = ifelse("M" %in% names(fittingdata@fixedPars), fittingdata@fixedPars$M, fittingdata@lmOutput$par$M)
        fittingdata@parentPeakSize = ifelse("S" %in% names(fittingdata@fixedPars), fittingdata@fixedPars$S, fittingdata@lmOutput$par$S)
        fittingdata@generationsDistance = ifelse("D" %in% names(fittingdata@fixedPars), fittingdata@fixedPars$D, fittingdata@lmOutput$par$D)
    } else {
        fittingdata@lmOutput <- nls.lm(
            par=fittingdata@parStart,
            fn = fittingdata@residFun,
            N = fittingdata@modelPoints$y,
            xx = fittingdata@modelPoints$x,
            control = nls.lm.control
            (
                nprint = .verbosevar,
                maxiter = 1024,
                factor = 0.1,
                maxfev=1000000,
                ptol=.Machine$double.xmin,
                gtol= 0,
                ftol=.Machine$double.xmin
            )
        )
        fittingdata@parentPeakPosition  = fittingdata@lmOutput$par$M
        fittingdata@parentPeakSize      = fittingdata@lmOutput$par$S
        fittingdata@generationsDistance = fittingdata@lmOutput$par$D
    }

    fittingdata@fittingDeviance     = fittingdata@lmOutput$deviance
    # I take the absolute values (it is a^2 ...)
    fittingdata@heights             = as.list(abs(coef(fittingdata@lmOutput))[1:fittingdata@numberOfPeaks])
    fittingdata@generations         = .generationsPercent(fittingdata)

    # RETURN
    return(fittingdata)
}