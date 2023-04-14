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
colnames(data) = paste("level", 1:nLevels, sep="")


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