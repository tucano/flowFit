## ===========================================================================
##  Check Class helper
##  FIXME or REMOVE is just a copy of flowFrame method
## ---------------------------------------------------------------------------
## Check for the class of object x and its length and cast error if wrong
## this is a copy and paste from flowFrame ....
.checkClass <- function(x, class, length=NULL, verbose=FALSE,
                       mandatory=TRUE)
{
    if(mandatory && missing(x))
        stop("Argument '", substitute(x), "' missing with no default",
             call.=verbose)
    msg <- paste("'", substitute(x), "' must be object of class ",
                 paste("'", class, "'", sep="", collapse=" or "), sep="")
    fail <- !any(sapply(class, function(c, y) is(y, c), x))
    if(!is.null(length) && length(x) != length)
    {
        if(!is.null(x))
        {
            fail <- TRUE
            msg <- paste(msg, "of length", length)
        }
    }
    if(fail) stop(msg, call.=verbose) else invisible(NULL)
}

## ===========================================================================
## check if we are in RStudio to access manipulate
## ---------------------------------------------------------------------------
.isRStudio <- function()
{
  return(any(grepl("RStudio", .libPaths())))
}

## ===========================================================================
## getKeywordValue
## get the keyword value from flowFrame
## ---------------------------------------------------------------------------
.getKeywordValue <- function(flowframe, channel, suffix)
{
    .channel.position <- match(channel, keyword(flowframe))
    # in some cases (flowJo), the channel name is without "<,>" so try again without "<"
    if (is.na(.channel.position))
    {
        .channel.position <- match(gsub("[<>]","", channel), keyword(flowframe))
    }

    # still na?
    if (is.na(.channel.position))
    {
        msg <- paste("Channel", channel, "is not present in flowFrame keywords!", sep=" ")
        stop(msg)
    }
    .channel.name <- names(keyword(flowframe))[.channel.position]

    # extract first part (remove N and add R)
    .channel.rname <- paste(sub("N","", .channel.name), suffix, sep="")

    # check if exists
    if (.channel.rname %in% names(keyword(flowframe))) {
        .my.value <- get(.channel.rname, keyword(flowframe))
        return(as.numeric(.my.value))
    } else {
        return(NA)
    }
}

## ===========================================================================
## setDataRange
##
## ---------------------------------------------------------------------------
.setDataRange <- function(flowframe, channel)
{
    .channel.position = which(colnames(flowframe) == channel)
    return(flowframe@parameters@data$maxRange[.channel.position])
}

## ===========================================================================
## setLogDecades
##
## ---------------------------------------------------------------------------
.setLogDecades <- function(flowframe, channel)
{
    .logDecades = .getKeywordValue(flowframe, channel, "M")
    .channel.position = which(colnames(flowframe) == channel)
    if (is.na(.logDecades)) {
        message("Setting logDecades using acquisition resolution")
        .logDecades = log10(flowframe@parameters@data$range[.channel.position])
    } else {
        message("Setting logDecades from keywords")
    }
    return(.logDecades)
}

## ===========================================================================
## FIXED MODEL HELPERS
## ---------------------------------------------------------------------------
.setFixedValues <- function(object, fixedPars, verbose)
{

    # SET FIXED PARS: [M,S]
    if (length(fixedPars) == 1 && is.na(fixedPars)) {
        stop("To model a fixed Model, I need some fixed parameters, fixedPars is empty!")
    } else {
        object@fixedPars  = fixedPars
    }
    object@parStart = object@parStart[!(names(object@parStart) %in% names(object@fixedPars))]

    # now examine values of fixedPars (M and S)
    if ("S" %in% names(object@fixedPars) && is.na(object@fixedPars$S))
    {
        if (verbose) message(paste("FixedModel: Using estimated peak Size: ", object@estimatedPeakSize))
        object@fixedPars$S = object@estimatedPeakSize
    }
    if ("M" %in% names(object@fixedPars) && is.na(object@fixedPars$M))
    {
        if (verbose) message(paste("FixedModel: Using estimated peak Position: ", object@estimatedPeakPosition))
        object@fixedPars$M = object@estimatedPeakPosition
    }
    if ("D" %in% names(object@fixedPars) && is.na(object@fixedPars$D))
    {
        if (verbose) message(paste("FixedModel: Using estimated generations Distance: ", object@estimatedDistance))
        object@fixedPars$D = object@estimatedDistance
    }
    return(object)
}

.getModelType <- function(object)
{
    .msgstr <- ""
    if (object@fixedModel) {
        .msgstr <- "* Fixed Model with locked parameters:"
        .parameters <- vector(mode="character")
        if ("M" %in% names(object@fixedPars)) {
            .parameters <- c(.parameters, "Peak Position")
        }
        if ("S" %in% names(object@fixedPars)) {
            .parameters <- c(.parameters, "Peak Size")
        }
        if ("D" %in% names(object@fixedPars)) {
            .parameters <- c(.parameters, "Generations Distance")
        }
        .msgstr <- c(.msgstr, paste(.parameters, collapse=", "))
    } else {
        .msgstr <- "* Unfixed model."
    }
    return(.msgstr)
}

## ===========================================================================
## MODEL FORMULAS
## ---------------------------------------------------------------------------

## ---------------------------------------------------------------------------
## GENERATIONS
## ---------------------------------------------------------------------------
.generation1  <- function(parSS,xx) { parSS$a^2*exp(-((xx-parSS$M)^2)/(2*parSS$S^2)) }
.generation2  <- function(parSS,xx) { parSS$b^2*exp(-((xx-(parSS$M - parSS$D))^2)/(2*parSS$S^2)) }
.generation3  <- function(parSS,xx) { parSS$c^2*exp(-((xx-(parSS$M - (2*parSS$D)))^2)/(2*parSS$S^2)) }
.generation4  <- function(parSS,xx) { parSS$d^2*exp(-((xx-(parSS$M - (3*parSS$D)))^2)/(2*parSS$S^2)) }
.generation5  <- function(parSS,xx) { parSS$e^2*exp(-((xx-(parSS$M - (4*parSS$D)))^2)/(2*parSS$S^2)) }
.generation6  <- function(parSS,xx) { parSS$f^2*exp(-((xx-(parSS$M - (5*parSS$D)))^2)/(2*parSS$S^2)) }
.generation7  <- function(parSS,xx) { parSS$g^2*exp(-((xx-(parSS$M - (6*parSS$D)))^2)/(2*parSS$S^2)) }
.generation8  <- function(parSS,xx) { parSS$h^2*exp(-((xx-(parSS$M - (7*parSS$D)))^2)/(2*parSS$S^2)) }
.generation9  <- function(parSS,xx) { parSS$i^2*exp(-((xx-(parSS$M - (8*parSS$D)))^2)/(2*parSS$S^2)) }
.generation10 <- function(parSS,xx) { parSS$j^2*exp(-((xx-(parSS$M - (9*parSS$D)))^2)/(2*parSS$S^2)) }
.generation11 <- function(parSS,xx) { parSS$k^2*exp(-((xx-(parSS$M - (10*parSS$D)))^2)/(2*parSS$S^2)) }
.generation12 <- function(parSS,xx) { parSS$l^2*exp(-((xx-(parSS$M - (11*parSS$D)))^2)/(2*parSS$S^2)) }
.generation13 <- function(parSS,xx) { parSS$m^2*exp(-((xx-(parSS$M - (12*parSS$D)))^2)/(2*parSS$S^2)) }
.generation14 <- function(parSS,xx) { parSS$n^2*exp(-((xx-(parSS$M - (13*parSS$D)))^2)/(2*parSS$S^2)) }
.generation15 <- function(parSS,xx) { parSS$o^2*exp(-((xx-(parSS$M - (14*parSS$D)))^2)/(2*parSS$S^2)) }
.generation16 <- function(parSS,xx) { parSS$p^2*exp(-((xx-(parSS$M - (15*parSS$D)))^2)/(2*parSS$S^2)) }
.generation17 <- function(parSS,xx) { parSS$q^2*exp(-((xx-(parSS$M - (15*parSS$D)))^2)/(2*parSS$S^2)) }
.generation18 <- function(parSS,xx) { parSS$r^2*exp(-((xx-(parSS$M - (16*parSS$D)))^2)/(2*parSS$S^2)) }
.generation19 <- function(parSS,xx) { parSS$s^2*exp(-((xx-(parSS$M - (17*parSS$D)))^2)/(2*parSS$S^2)) }
.generation20 <- function(parSS,xx) { parSS$t^2*exp(-((xx-(parSS$M - (18*parSS$D)))^2)/(2*parSS$S^2)) }

## ---------------------------------------------------------------------------
## PARENT MODEL
## ---------------------------------------------------------------------------
.parent.model <- function(parSS,xx) { .generation1(parSS,xx) }
.parent.model.residFun <- function(p,N,xx) { N - .parent.model(p,xx) }

## ---------------------------------------------------------------------------
## FIXED PARENT MODEL
## ---------------------------------------------------------------------------
.fixed.parent.model <- function(parSS,xx,fixedPars)
{
    pp <- c(parSS, fixedPars)
    return(.generation1(pp,xx))
}
.fixed.parent.model.residFun <- function(parSS,N,xx,fixedPars)
{
    N - .fixed.parent.model(parSS,xx,fixedPars)
}

## ---------------------------------------------------------------------------
## MODEL BUILDER
## ---------------------------------------------------------------------------
.buildModel <- function(generations)
{
    .snippets <- vector(mode="character", generations)
    for (i in 1:generations)
    {
        .fs  <- paste(".generation",i,sep="")
        .snippets[i] <- as.character(body(get(.fs)))[2]
    }
    .modelString <- paste(.snippets, collapse=" + ")
    # generate the parsed body function
    .modelParsed <- parse(text=.modelString)
    # make an empty function
    .model <- function(){}
    # all the vars
    .v <- c("parSS","xx")
    # construct arguments
    .call <- do.call("alist",as.list(rep(TRUE,length(.v))))
    names(.call) <- .v
    formals(.model) <- .call
    body(.model) <- .modelParsed
    return(.model)
}

.buildFixedModel <- function(generations)
{
    .model <- .buildModel(generations)
    .fixed.model <- function(parSS,xx,fixedPars)
    {
        pp <- c(parSS, fixedPars)
        return(.model(pp,xx))
    }
    return(.fixed.model)
}

## ---------------------------------------------------------------------------
## PERCENTAGE of CELL FOR GENERATION with INTEGRATE
## ---------------------------------------------------------------------------
.generationsPercent <- function(fittingdata)
{
    .mypars <- as.list(coef(fittingdata@lmOutput))

    if (fittingdata@fixedModel) {
        .all.cells <- integrate(fittingdata@model, lower=0, upper=(fittingdata@dataRange-1), parSS=.mypars, fixedPars=fittingdata@fixedPars)$value
    } else {
        .all.cells <- integrate(fittingdata@model, lower=0, upper=(fittingdata@dataRange-1), parSS=.mypars)$value
    }

    .generations <- vector(mode="numeric",fittingdata@numberOfPeaks)

    .m <- ifelse("M" %in% names(fittingdata@fixedPars), fittingdata@fixedPars$M, .mypars$M)
    .s <- ifelse("S" %in% names(fittingdata@fixedPars), fittingdata@fixedPars$S, .mypars$S)
    .d <- ifelse("D" %in% names(fittingdata@fixedPars), fittingdata@fixedPars$D, .mypars$D)

    .this.generation.pars <- c(.mypars[1], M=.m, S=.s)
    .generations[1] <- integrate(.generation1, lower=0, upper=(fittingdata@dataRange), parSS=.this.generation.pars)$value

    # Integrate the rest
    for (i in 2:fittingdata@numberOfPeaks)
    {
        .fx <- get(paste(".generation",i,sep=""))
        .this.generation.pars <- c(.mypars[i], M=.m, S=.s, D=.d)
        .generations[i] <- integrate(.fx, lower=0, upper=(fittingdata@dataRange), parSS=.this.generation.pars)$value
    }
    return((.generations/.all.cells)*100)
}

##############################################################################
# Generations Valid Checker
##############################################################################
.genIsValid <- function(mypar)
{
    if (mypar <= 0) {
        return(FALSE)
    } else {
        return(TRUE)
    }
}
