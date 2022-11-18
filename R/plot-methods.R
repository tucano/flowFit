## ---------------------------------------------------------------------------
## Generic function for plotting of R objects
## ---------------------------------------------------------------------------

#############################
# parentFittingData: plot
#############################
setMethod(
  "plot",
  signature(x = "parentFittingData"),
  function(x, y,
    main=NULL,
    xlab=NULL,
    ylab=NULL,
    legend=TRUE,
    logScale=TRUE,
    drawGrid=TRUE)
  {
    # legend is TRUE, we must set a layout
    if (legend) {
      def.par <- par(no.readonly = TRUE) # save default, for resetting...
      layout(rbind(1,2), heights=c(7,1))  # put legend on bottom 1/8th of the chart
    }

    # add main, xlab and ylab if NULL
    if (is.null(xlab)) xlab = x@channel
    if (is.null(ylab)) ylab = "Events"
    if (is.null(main)) main = "Parent Population"

    # generate a diagnostic message
    msg1 <- c(
      paste("Events: ", length(which(x@dataPoints$y > 0))),
      paste("Deviance:", round(x@fittingDeviance)),
      paste("Parent Peak Position:", round(x@parentPeakPosition, digits=2) ),
      paste("Parent Peak Size:", round(x@parentPeakSize, digits=2) )
    )

    # MAIN PLOT
    # changed to xlim=c(0,dataRange)  -1 is a problem for log values ...
    plot(x@modelPoints,
      main=main,
      xlab=xlab,
      ylab=ylab,
      type="l",
      lwd=2,
      ylim=c(0,max(x@modelPoints$y)),
      xlim=c(0,x@dataRange),
      axes=FALSE)

    # ADD MODEL POINTS
    if (x@fixedModel) {
      .fittedValues <- .fixed.parent.model(as.list(coef(x)), x@dataPoints$x, x@fixedPars)
    } else {
      .fittedValues <- .parent.model(as.list(coef(x)), x@dataPoints$x)
    }
    lines(x@dataPoints$x, .fittedValues, col="red", lwd=1, lty=1)

    # AXIS
    axis(2)
    if (logScale) {
      .ticks <- logTicks(x@dataRange,x@logDecades)
      axis(1,.ticks$all,label=FALSE)
      axis(1,at=.ticks$major,labels=.ticks$labels)
    } else {
      axis(1)
    }

    if (drawGrid)
      proliferationGrid(x@parentPeakPosition, dataRange = x@dataRange, logDecades = x@logDecades)

    if(legend) {
      # setup for no margins on the legend
      par(mar=c(0, 0, 0, 0))
      plot.new()
      legend("center", msg1, col=black, ncol=2, bty="n")
    }

    if (legend) par(def.par)#- reset to default
  }
)

#################################
# proliferationFittingData: plot
#################################
setMethod(
  "plot",
  signature(x="proliferationFittingData"),
  function(x,y,
    which="all",
    main=NULL,
    xlab=NULL,
    ylab=NULL,
    legend=TRUE,
    logScale = TRUE,
    drawGrid = TRUE)
  {
    # add main, xlab and ylab if NULL
    if (is.null(xlab)) xlab = x@channel
    if (is.null(ylab)) ylab = "Events"
    if (is.null(main)) main = "Proliferating Population"

    # create which list
    which.plot <- rep(FALSE,5)
    if (is.character(which) && which == "all") {
      which.plot[1:length(which.plot)] = TRUE
    } else {
      which.plot[which] = TRUE
    }

    # IF device is interactive or None ask for next plot
    if (dev.interactive(orNone = TRUE) && length(which.plot[which.plot] == TRUE) == length(which.plot)) par(ask = T)

    # colors FOR 20 generations
    mycolors <- rainbow(20)
    mycolors.alpha <- rainbow(20, alpha=0.6)
    # shrink vectors with the current number of peaks.
    # If we compare different samples: same color map
    mycolors <- mycolors[1:x@numberOfPeaks]
    mycolors <- mycolors.alpha[1:x@numberOfPeaks]
    # now generate labels
    mylabels <- paste("G", seq(1,x@numberOfPeaks), sep="")

    ## get the fitted values
    if (x@fixedModel) {
      .fittedValues <- x@model(x@lmOutput$par, x@modelPoints$x, x@fixedPars)
    } else {
      .fittedValues <- x@model(x@lmOutput$par, x@modelPoints$x)
    }

    if (legend) def.par <- par(no.readonly = TRUE) # save default, for resetting...

    # PLOT 1: THE DATA
    if (which.plot[1])
    {

      # legend is TRUE, we must set a layout
      if (legend) layout(rbind(1,2), heights=c(7,1))  # put legend on bottom 1/8th of the chart

      if (x@dataSmooth)
      {
        # DATA PLOT
        # changed to xlim=c(0,dataRange)  -1 is a problem for log values ...
        plot(x@dataPoints,
          main=main,
          xlab=xlab,
          ylab=ylab,
          type="l",
          lwd=2,
          lty=1,
          axes=FALSE,
          xlim=c(0,x@dataRange),
          ylim=c(0,max(x@modelPoints$y)))
        # SMOOTHED PLOT
        lines(x@modelPoints, col="blue", lwd=1, lty=1)

      } else {
        # changed to xlim=c(0,dataRange)  -1 is a problem for log values ...
        plot(x@dataPoints,
          main=main,
          xlab=xlab,
          ylab=ylab,
          type="l",
          axes=FALSE,
          xlim=c(0,x@dataRange),
          ylim=c(0,max(x@modelPoints$y)))
      }

      # AXIS
      axis(2)
      if (logScale) {
        .ticks <- logTicks(x@dataRange,x@logDecades)
        axis(1,.ticks$all,label=FALSE)
        axis(1,at=.ticks$major,labels=.ticks$labels)
      } else {
        axis(1)
      }

      # LEGEND
      if(legend) {
        # setup for no margins on the legend
        par(mar=c(0, 0, 0, 0))
        plot.new()
        if (x@dataSmooth) {
          legend("center", c("Input","Smooth"), col=c("black","blue"), lwd=2, lty=1, pch=-1, ncol=2, bty="n")
        } else {
          legend("center", c("Input"), col=c("black"), lwd=2, lty=1, pch=-1, ncol=2, bty="n")
        }
        # reset margins
        par(mar=def.par$mar)
      }

    }

    # PLOT 2: FITTING
    if (which.plot[2])
    {
      # legend is TRUE, we must set a layout
      if (legend) layout(rbind(1,2), heights=c(7,1))  # put legend on bottom 1/8th of the chart

      # MAIN FITTING PLOT
      # changed to xlim=c(0,dataRange)  -1 is a problem for log values ...
      # PLOT DATA
      plot(x@modelPoints,
        main = main,
        xlab=xlab,
        ylab=ylab,
        type="l",
        lwd=2,
        col="black",
        xlim=c(0,x@dataRange),
        ylim=c(0,max(x@modelPoints$y)),
        axes=FALSE)
      lines(x@modelPoints$x, .fittedValues, col="red", lwd=1, lty=1)

      # AXIS
      axis(2)
      if (logScale) {
        .ticks <- logTicks(x@dataRange,x@logDecades)
        axis(1,.ticks$all,label=FALSE)
        axis(1,at=.ticks$major,labels=.ticks$labels)
      } else {
        axis(1)
      }

      # GRID
      if (drawGrid) {
        proliferationGrid(x@parentPeakPosition,
          fittedDistance=x@generationsDistance,
          dataRange = x@dataRange,
          logDecades = x@logDecades)
      }

      # LEGEND
      if (legend) {
        # setup for no margins on the legend
        par(mar=c(0, 0, 0, 0))
        plot.new()
        legend("center", c("Input","Fitting"), col=c("black","red"), lwd=2, lty=1, pch=-1, ncol=2, bty="n")
        # reset margins
        par(mar=def.par$mar)
      }
    }

    # PLOT 3: SINGLE PEAKS PLOT
    if (which.plot[3])
    {
      # legend is TRUE, we must set a layout
      if (legend) layout(rbind(1,2), heights=c(7,1))  # put legend on bottom 1/8th of the chart

      # changed to xlim=c(0,dataRange)  -1 is a problem for log values ...
      plot(x@modelPoints,
        lwd=2,
        type="l",
        main=main,
        axes=FALSE,
        xlab=xlab,
        ylab=ylab,
        col="black",
        ylim=c(0,max(x@dataPoints$y)),
        xlim=c(0,x@dataRange))
      # FITTING
      lines(x@modelPoints$x, .fittedValues, col="red", lwd=1, lty=1)

      for (i in 1:length(x@heights))
      {
        if (.genIsValid(x@heights[[i]]))
        {
          .my.function.name <- paste('.generation',i,sep="")
          .my.pars <- c(x@heights[i], M=x@parentPeakPosition, S=x@parentPeakSize)
          if (i > 1) .my.pars <- c(.my.pars, D=x@generationsDistance)
          .my.args <- list(xx = x@dataPoints$x, parSS=.my.pars)
          .fittedGen <- do.call(.my.function.name, .my.args)
          # LINE AND POLY
          lines(x@modelPoints$x,.fittedGen, col=mycolors[i], lwd=1, lty=1)
          polygon(x@modelPoints$x,.fittedGen, col=mycolors.alpha[i], border=NA)
        }
      }

      # AXIS
      axis(2)
      if (logScale) {
        .ticks <- logTicks(x@dataRange,x@logDecades)
        axis(1,.ticks$all,label=FALSE)
        axis(1,at=.ticks$major,labels=.ticks$labels)
      } else {
        axis(1)
      }

      # grid lines
      if (drawGrid) {
        proliferationGrid(x@parentPeakPosition,
          fittedDistance=x@generationsDistance,
          dataRange = x@dataRange,
          logDecades = x@logDecades)
      }

      # LEGEND
      if (legend) {
        # setup for no margins on the legend
        par(mar=c(0, 0, 0, 0))
        plot.new()
        legend("center", mylabels,
          col=mycolors, lty=1, lwd=3, pch=-1, merge=TRUE, ncol=floor(x@numberOfPeaks/2), bty="n")
        # reset margins
        par(mar=def.par$mar)
      }
    }

    # PLOT 4: BARPLOT GENERATIONS
    if (which.plot[4])
    {
      # if legend is TRUE reset layout
      if (legend) layout(1)

      mp <- barplot2(x@generations,
        main=main,
        xlab="Generations",
        ylab="% of cells",
        col=mycolors,
        ylim=c(0,100),
        axes=FALSE,
        axisnames=FALSE)
      axis(2,seq(0,100,20))
      axis(1, at = mp, labels=1:x@numberOfPeaks, las=2)
      # add labels with percentage rounded only if legend = TRUE
      if (legend)
      {
        .positions <- x@generations
        .labels <- paste(round(x@generations, digits=2), "%", sep="")
        text(mp, .positions, .labels, pos=3, cex=0.7)
      }
    }

    # PLOT 5: MINPACK LM DIAG
    if (which.plot[5])
    {
      msg1 <- c(
        paste("Events: ", length(which(x@dataPoints$y > 0))),
        paste("Deviance:", round(x@lmOutput$deviance)),
        paste("[Distance] model: ", round(x@generationsDistance, digits=2), ", guess: " , round(x@estimatedDistance, digits=2)),
        paste("[Parent Peak Position] model: ", round(x@parentPeakPosition, digits=2), ", guess: " , round(x@estimatedPeakPosition, digits=2)),
        paste("[Parent Peak Size] model: ", round(x@parentPeakSize, digits=2), ", guess: " , round(x@estimatedPeakSize, digits=2) )
      )

      if (legend) layout(rbind(1,2), heights=c(5,3))

      plot(1:(x@lmOutput$niter+1), log(x@lmOutput$rsstrace),
        type="b",
        main=main,
        xlab="Iteration",
        ylab="Log residual sum of squares",
        pch=21,
        bg=2)

      if (legend) {
        # setup for no margins on the legend
        par(mar=c(0, 0, 0, 0))
        plot.new()
        legend("center", msg1, col=black, ncol=1, bty="n")
        # reset margins
        par(mar=def.par$mar)
      }
    }
    if (legend) par(def.par)#- reset to default
  }
)