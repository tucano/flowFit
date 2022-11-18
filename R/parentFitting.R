## =========================================================================##
## =========================================================================##
##                       parentFitting method                               ##
##                     fitting a parent population                          ##
## =========================================================================##
## =========================================================================##

parentFitting <- function(
    flowframe,
    channel,
    estimatedPeakPosition = NA,
    estimatedPeakSize = NA,
    dataRange = NA,
    logDecades = NA,
    binning = TRUE,
    breaks = 1024,
    dataSmooth = TRUE,
    smoothWindow = 2,
    fixedModel = FALSE,
    fixedPars = NA,
    verbose = FALSE)
{
    parentdata <- new("parentFittingData", flowframe, channel, binning=binning, breaks=breaks, dataSmooth=dataSmooth, smoothWindow=smoothWindow)

    # FIXED MODEL SETUP SLOTS HERE
    parentdata@fixedModel = fixedModel

    #check for estimatedPeakPosition & estimatedPeakSize
    if (is.na(estimatedPeakPosition)) {
        parentdata@estimatedPeakPosition = mean(exprs(flowframe)[,channel])
        if (verbose)
            message(paste("parentFitting: No estimated peak position provided. Setting position from data mean:", parentdata@estimatedPeakPosition, sep=" "))
    } else {
        parentdata@estimatedPeakPosition = estimatedPeakPosition
    }
    if (is.na(estimatedPeakSize)) {
        parentdata@estimatedPeakSize = sd(exprs(flowframe)[,channel])
        if (verbose)
            message(paste("parentFitting: No estimated peak size provided. Setting size from data standard deviation:", parentdata@estimatedPeakSize, sep=" "))
    } else {
        parentdata@estimatedPeakSize = estimatedPeakSize
    }

    # set dataRange and log decades
    if (is.na(dataRange)) {
        parentdata@dataRange = .setDataRange(flowframe, channel)
        message(paste("parentFitting: No dataRange provided. Setting range from channel range:", parentdata@dataRange, sep=" "))
    } else {
        parentdata@dataRange = dataRange
    }

    if (is.na(logDecades)) {
        parentdata@logDecades = .setLogDecades(flowframe, channel)
        message(paste("proliferationFitting: No logDecades provided. Setting LOG dynamic range from flowFrame: log decades: ", parentdata@logDecades , sep="" ))
    } else {
        parentdata@logDecades = logDecades
    }

    # set parStart
    parentdata@parStart = list(a=1, M=parentdata@estimatedPeakPosition, S=parentdata@estimatedPeakSize)

    # FIXED MODEL
    if (parentdata@fixedModel) parentdata <- .setFixedValues(parentdata, fixedPars, verbose)

    # set the model and the data
    parentdata@dataMatrix <- exprs(parentdata@data[,channel])

    if (parentdata@fixedModel) {
        parentdata@residFun = .fixed.parent.model.residFun
    } else {
        parentdata@residFun = .parent.model.residFun
    }

    # BINNING AND SMOOTHING: dataPoints
    if (parentdata@binning) {
        if (verbose)
            message(paste("parentFitting: Binning data into:", parentdata@breaks, "breaks", sep=" "))
        h <- hist(parentdata@dataMatrix, breaks=parentdata@breaks, plot=FALSE)
        y.points <- h$counts
        x.points <- h$mids
        parentdata@dataPoints <- data.frame( x = x.points, y = y.points)
    } else {
        if (verbose)
            message(paste("parentFitting: Tabulating data into:", parentdata@dataRange, "FACS data points", sep=" "))
        .y.points <- tabulate(parentdata@dataMatrix, nbins=parentdata@dataRange)
        parentdata@dataPoints <- data.frame( x = 0:(parentdata@dataRange-1), y = .y.points)
    }

    # SMOOTHING: modelPoints
    if (parentdata@dataSmooth) {
        if (verbose)
            message("parentFitting: smoothing data with Kolmogorov-Zurbenko low-pass linear filter.")
        .y.smooth <- kz(parentdata@dataPoints$y, parentdata@smoothWindow)
        parentdata@modelPoints <- data.frame(x = parentdata@dataPoints$x, y = .y.smooth)
    } else {
        parentdata@modelPoints <- parentdata@dataPoints
    }

    .verbosevar <- ifelse(verbose, 1, 0)

    if (parentdata@fixedModel) {
        parentdata@lmOutput <- nls.lm(
        par=parentdata@parStart,
        fn = parentdata@residFun,
        N = parentdata@modelPoints$y,
        xx = parentdata@modelPoints$x,
        fixedPars = parentdata@fixedPars,
        control = nls.lm.control(
            nprint = .verbosevar,
            maxiter = 1024,
            factor = 0.1,
            maxfev=1000000,
            ptol=.Machine$double.xmin,
            gtol= 0,
            ftol=.Machine$double.xmin)
        )
        parentdata@parentPeakPosition = ifelse("M" %in% names(parentdata@fixedPars), parentdata@fixedPars$M, parentdata@lmOutput$par$M)
        parentdata@parentPeakSize = ifelse("S" %in% names(parentdata@fixedPars), parentdata@fixedPars$S, parentdata@lmOutput$par$S)
    } else {
        parentdata@lmOutput <- nls.lm(
        par=parentdata@parStart,
        fn = parentdata@residFun,
        N = parentdata@modelPoints$y,
        xx = parentdata@modelPoints$x,
        control = nls.lm.control(
            nprint = .verbosevar,
            maxiter = 1024,
            factor = 0.1,
            maxfev=1000000,
            ptol=.Machine$double.xmin,
            gtol= 0,
            ftol=.Machine$double.xmin)
        )
        parentdata@parentPeakPosition = parentdata@lmOutput$par$M
        parentdata@parentPeakSize = parentdata@lmOutput$par$S
    }

    parentdata@fittingDeviance = parentdata@lmOutput$deviance

    return(parentdata)
}
