## =========================================================================##
## =========================================================================##
##                       Generic methods definitions                        ##
## =========================================================================##
## =========================================================================##

## ===========================================================================
## Used to summarize the operation of a filter on a frame (already S3 in base)
## ---------------------------------------------------------------------------
setGeneric("summary", function(object,...) standardGeneric("summary"))

## ===========================================================================
## Generic to get to the data linked to a view or actionItem
## The simplest and most common situation is that name is already an ordinary non-generic non-primitive
## function, and you now want to turn this function into a generic. In this case you will most often supply
## only name
## ---------------------------------------------------------------------------
setGeneric("Data", function(object) standardGeneric("Data"))

## ---------------------------------------------------------------------------
## coef is a generic function which extracts model coefficients from objects returned by
## modeling functions. coefficients is an alias for it.
## ---------------------------------------------------------------------------
setGeneric("coef", function(object) standardGeneric("coef"))

## ---------------------------------------------------------------------------
## Computes confidence intervals for one or more parameters in a fitted model.
## There is a default and a method for objects inheriting from class "lm".
## ---------------------------------------------------------------------------
setGeneric("confint", function(object, parm, level = 0.95, ... ) standardGeneric("confint"))

## ---------------------------------------------------------------------------
## Generic function for plotting of R objects
## ---------------------------------------------------------------------------
setGeneric("plot",function(x,y,...) standardGeneric("plot"))
