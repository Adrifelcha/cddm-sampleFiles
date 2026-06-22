<<<<<<< HEAD
######################################################################
##
## Based on 
##  A stackoverflow answer posted by Gregor Thomas
##    (See https://stackoverflow.com/questions/13327489/vertical-histogram)
##  which was then perfected by user ozgen on Github
##    (See https://github.com/ozgeneral/VerticalHist) 
######################################################################

# Make up data to test these functions
######################################################################
N = 200
nLevels = 4
data = matrix(rnorm(N*nLevels,1:nLevels,0.1), byrow = TRUE, ncol=nLevels)
data = as.data.frame(data)


# Draw vertical histogram
######################################################################
VerticalHist <- function(x, xscale = NULL, xwidth, hist,
                         fillCol = "gray80", lineCol = "gray40") {
  ## x - the x position of each histogram
  ## xscale - "height" of the tallest bar (horizontally),
  ## xwidth - horizontal spacing between histograms
  ## hist - an object of type "histogram" (i.e., with $breaks and $density)
  binWidth <- hist$breaks[2] - hist$breaks[1]
  if (is.null(xscale)) xscale <- xwidth * 0.90 / max(hist$density)
  n <- length(hist$density)
  
  right_x.l <- rep(x, n)
  right_x.r <- right_x.l + hist$density * xscale
  left_x.r <- rep(x, n)
  left_x.l <- left_x.r - hist$density * xscale
  x.l = c(left_x.l,right_x.l)
  x.r = c(left_x.r,right_x.r)
  
  y.b <- hist$breaks[1:n]
  y.b = rep(y.b,2)
  y.t <- hist$breaks[2:(n + 1)]
  y.t = rep(y.t,2)
  
  rect(xleft = x.l, ybottom = y.b, xright = x.r, ytop = y.t,
       col = fillCol, border = fillCol)
}

# Main plotting function
######################################################################
makePlot = function(data, binWidth = 0.5, ylab="Parameter name", fillCol="indianred4"){
    numberOfHists <- ncol(data)
    binStarts <- 1:numberOfHists
    binMids <- binStarts + binWidth / 2
    axisCol <- "gray50"
      
    ## Data handling
    allValues <- unlist(as.list(data))
    DOYrange <- range( allValues, na.rm = TRUE )
    DOYrange <- c(floor(DOYrange[1]), ceiling(DOYrange[2]))
    DOYmean <- apply(data,2,mean)
    quantiles <- round(quantile( DOYrange, c(0.2, 0.4, 0.6, 0.8), na.rm=TRUE ), digits=0)
    gridlines <- round(quantile( DOYrange, c(0.1, 0.3, 0.5, 0.7, 0.9), na.rm=TRUE ), digits=0)

        ## Get the histogram obects
    histList <- apply( data, 2, function(x, hCol) hist(x, plot = FALSE))
      
    ## Plotting
    xlim <- c(0,numberOfHists+1)
    ylim <- DOYrange
    xlab=""
    
    plot(c(0, 5), DOYrange, type = "n", xlim=xlim, ylim=ylim,
         ann = FALSE, axes = FALSE, xaxs = "i", yaxs = "i")
    axis(1, 1:numberOfHists, paste("Level",1:numberOfHists,sep=" "),
         cex.axis = 1.2, col = axisCol)
    mtext(side = 1, outer = F, line = 3, xlab, cex = 1.2)
    y.seq = format(round(seq(DOYrange[1],DOYrange[2],length.out=10),digits = 2), nsmall = 1)
    axis(2, cex.axis = 0.95, las = 1, line = -.7, col = "white", tck = 0,
         at = y.seq, labels = y.seq, las=2)
    mtext(side = 2, outer = F, line = 2.8, ylab, cex = 1.2)
    box(bty = "o", col = axisCol)
      
    biggestDensity <- max(unlist(lapply(histList, function(h){max(h[[4]])})))
    xscale <- binWidth * .9 / biggestDensity
      
    ## Plot the histograms
    for (lengthBin in 1:numberOfHists) {
         VerticalHist(x = binStarts[lengthBin], 
                      xscale = xscale, 
                      xwidth = binWidth, 
                      hist= histList[[lengthBin]], 
                      fillCol = fillCol)
    }
}
  
makePlot(data)
=======
# Process posterior samples
######################################################################
source("./plotFunctions.R")

xlab <- "Levels"
legend.size <- 1.3
legend.col <- "gray40"
#########################################################
# Figure 1: Results of CARTESIAN recovery study
#########################################################
output.File <- "../output_cart_n200/"
load(paste(output.File,"sim_trueValues.RData",sep=""))
load(paste(output.File,"sim_meanPosteriors.RData",sep=""))

# retrievedValues_restricted <- rmv.true.value(trueValues,
#                                              output.set = retrievedValues,
#                                              bad.value = c(0.0,0.0),
#                                              rmv.parameter = c("mu1","mu2"))
# trueValues_restricted <- rmv.true.value(trueValues,
#                                         output.set = trueValues,
#                                         bad.value = c(0.0,0.0),
#                                         rmv.parameter = c("mu1","mu2"))


mu1 <- group.by.level(trueValues, retrievedValues,
                      parameter.name="mu1")
mu2 <- group.by.level(trueValues, retrievedValues,
                      parameter.name="mu2")
ndt <- group.by.level(trueValues, retrievedValues,
                      parameter.name="ndt")
bound <- group.by.level(trueValues, retrievedValues,
                        parameter.name="bound")

par(mfrow=c(2,2))
par(mar = c(2.5, 4, 1, 1))
makePlot(true.vector = bound$par.values,  
         retrieved.matrix=bound$groups, binWidth = 0.13, 
         ylab=expression(paste("Boundary radius (",eta,")")),
         #group.Names = c(expression(paste(eta, " = ", 1.5)),
         #                expression(paste(eta, " = ", 2.0)),
         #                expression(paste(eta, " = ", 2.5))),
         ylim=c(1,3.1), fillCol = "#3BA9A7", xlab=xlab,
         hor.col = legend.col)
legend(0.6,3.1, "True values", lty = 2, cex = legend.size, col=legend.col)

makePlot(true.vector = ndt$par.values,  
         retrieved.matrix=ndt$groups, binWidth = 0.008, 
         ylab=expression(paste("Nondecision time (",tau,")")),
         #group.Names = c(expression(paste(tau, " = ", 0.1)),
         #                expression(paste(tau, " = ", 0.2)),
         #                expression(paste(tau, " = ", 0.3))),
         ylim=c(-0.01,0.7), fillCol = "#3B77A9", 
         hor.col = legend.col)

makePlot(true.vector = mu1$par.values,  
         retrieved.matrix=mu1$groups, binWidth = 0.06, 
         ylab=expression(paste("Step size on X (",mu[x],")")),
         #group.Names = c(expression(paste(mu[1], " = ", -0.5)),
         #                expression(paste(mu[1], " = ", 0.0)),
         #                expression(paste(mu[1], " = ", 0.5)),
         #                expression(paste(mu[1], " = ", 1.0))),
         ylim=c(-1,1.5), fillCol = "#653BA9",
         hor.col = legend.col)

makePlot(true.vector = mu2$par.values,  
         retrieved.matrix=mu2$groups, binWidth = 0.04, 
         ylab=expression(paste("Step size on Y (",mu[y],")")),
         #group.Names = c(expression(paste(mu[2], " = ", -1.0)),
        #                 expression(paste(mu[2], " = ", -0.5)),
        #                 expression(paste(mu[2], " = ", 0.0)),
        #                 expression(paste(mu[2], " = ", 0.5))),
        ylim=c(-1.5,1), fillCol = "#3B3BA9",
        hor.col = legend.col)




#########################################################
# Figure 2: Results of POLAR recovery study
#########################################################
output.File <- "../output_polar_n200/"
load(paste(output.File,"sim_trueValues.RData",sep=""))
load(paste(output.File,"sim_meanPosteriors.RData",sep=""))

retrievedValues_restricted <- rmv.true.value(trueValues,
                                             output.set = retrievedValues,
                                             bad.value = 0.01,
                                             rmv.parameter = "driftLength")
trueValues_restricted <- rmv.true.value(trueValues,
                                        output.set = trueValues,
                                        bad.value = 0.01,
                                        rmv.parameter = "driftLength")

driftAngle <- group.by.level(trueValues_restricted, retrievedValues_restricted,
                             parameter.name="driftAngle")
driftLength <- group.by.level(trueValues, retrievedValues,
                              parameter.name="driftLength")
ndt <- group.by.level(trueValues, retrievedValues,
                      parameter.name="ndt")
bound <- group.by.level(trueValues, retrievedValues,
                        parameter.name="bound")


par(mfrow=c(2,2))
par(mar = c(2.5, 4, 1, 1))
makePlot(true.vector = bound$par.values,  
         retrieved.matrix=bound$groups, binWidth = 0.13, 
         ylab=expression(paste("Boundary radius (",eta,")")),
         #group.Names = c(expression(paste(eta, " = ", 1.5)),
         #                expression(paste(eta, " = ", 2.0)),
         #                expression(paste(eta, " = ", 2.5))),
         ylim=c(1,3), fillCol = "#3BA9A7",
         internal.margin.X = 0.7,
         xlab=xlab, hor.col = legend.col)
legend(0.375,3, "True values", lty = 2, cex = legend.size, col = legend.col)

makePlot(true.vector = ndt$par.values,  
         retrieved.matrix=ndt$groups, binWidth = 0.008, 
         ylab=expression(paste("Nondecision time (",tau,")")),
         #group.Names = c(expression(paste(tau, " = ", 0.1)),
         #                expression(paste(tau, " = ", 0.2)),
         #                expression(paste(tau, " = ", 0.3))),
         ylim=c(-0.01,0.5), fillCol = "#3B77A9",
         internal.margin.X = 0.7, hor.col = legend.col)

makePlot(true.vector = driftLength$par.values,  
         retrieved.matrix=driftLength$groups, binWidth = 0.08, 
         ylab=expression(paste("Drift length (",delta,")")),
         #group.Names = c(expression(paste(delta, " = ", 0.01)),
         #                expression(paste(delta, " = ", 1.0)),
         #                expression(paste(delta, " = ", 2.0))),
         ylim=c(0,3), fillCol = "#ECC857",
         internal.margin.X = 0.7, hor.col = legend.col)

driftAngle$par.values[3] <- driftAngle$par.values[3] - 2*pi
makePlot(true.vector = driftAngle$par.values,  
         retrieved.matrix=driftAngle$groups, binWidth = 0.15, 
         ylab=expression(paste("Drift angle (",theta,")")),
         #group.Names = c(expression(paste(theta, " = ", 0)),
         #                expression(paste(theta, " = ", 2)),
         #                expression(paste(theta, " = ", 4))),
         ylim=c(-3.2,3), fillCol = "#D68E51",
         internal.margin.X = 0.7, hor.col = legend.col)
>>>>>>> c908bf4d141e5acacbb6f73ab978d5fd2d5b09bf
