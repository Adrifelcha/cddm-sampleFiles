######################################################################
######################################################################
#############        P L O T     F U N C T I O N S
######################################################################
######################################################################

# Process posterior samples
######################################################################
group.by.level <- function(trueValues.matrix, output.set,
                           parameter.name){
      par.levels <- table(trueValues.matrix[,parameter.name])
      par.values <- as.numeric(names(par.levels))
      groups <- matrix(NA, 
                       nrow = max(par.levels)*200,
                       ncol = length(par.values))
      for(i in 1:length(par.values)){
          a <- par.values[i]
          same.level <- trueValues.matrix[,parameter.name]==a
          groups[,i] <- output.set[,parameter.name,same.level]
      }
      output <- list("groups" = groups, "par.values" = par.values)
      return(output)
}

group.by.level.across.folders <- function(output.Folders,par.name){
  par.values <- NA
  groups <- NA
  for(f in 1:length(output.Folders)){
    output.File <- output.Folders[f]
    load(paste(output.File,"sim_trueValues.RData",sep=""))
    load(paste(output.File,"sim_meanPosteriors.RData",sep=""))
    
    if(par.name=="driftAngle"){
    retrievedValues <- rmv.true.value(trueValues,
                                      output.set = retrievedValues,
                                      bad.value = 0.01,
                                      rmv.parameter = "driftLength")
    trueValues_restricted <- rmv.true.value(trueValues,
                                            output.set = trueValues,
                                            bad.value = 0.01,
                                            rmv.parameter = "driftLength")
    trueValues <- trueValues_restricted
    }
    
    get <- group.by.level(trueValues, retrievedValues,
                          parameter.name=par.name)
    par.values <- c(par.values, get$par.values)
    groups <- cbind(groups, get$groups)
  }
  output <- list("groups" = groups[,-1], 
                 "par.values" = par.values[-1])
  return(output)
}

locate.par.sets <- function(trueValues, value, par.selected){
      output <- which(trueValues[,par.selected]==value)
      return(output)
}

rmv.true.value <- function(trueValues, output.set,
                           bad.value, rmv.parameter){
    no.RmvParameters <- length(rmv.parameter)
    bad.par.sets <- list()
    for(i in 1:no.RmvParameters){
        find <- locate.par.sets(trueValues,
                                value = bad.value[i],
                                par.selected = rmv.parameter[i])
        bad.par.sets <- append(bad.par.sets, find)
    }
    bad.par.sets <- unlist(bad.par.sets)
    joint.cases <- table(bad.par.sets)==no.RmvParameters
    joint.cases.locate <- as.logical(joint.cases)
    bad.par.sets <- as.numeric(names(joint.cases)[joint.cases.locate])
    
    d <- which(dim(output.set)==nrow(trueValues))
    w = length(dim(output.set))==length(dim(trueValues))
    if(w){
        x = unique(dim(output.set)==dim(trueValues))
        if(x){
            if(d==1){
              output <- output.set[-bad.par.sets,]
            }else{
              output <- output.set[,-bad.par.sets]      
            }
        }
    }else{
        if(d==1){
          output <- output.set[-bad.par.sets,,]
        }else{
          if(d==2){
            output <- output.set[,-bad.par.sets,]      
          }else{
            output <- output.set[,,-bad.par.sets]      
          }
        }
    }
    return(output)
}


#############################################################################
## Visualizing the distribution of parameter values recovered
##
## This code is based on a stackoverflow post by Gregor Thomas, 
##    (See https://stackoverflow.com/questions/13327489/vertical-histogram)
##    (See also: https://github.com/ozgeneral/VerticalHist) 
#############################################################################
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
  
  # Define horizontal limits of left/right bars
  right_x.l <- rep(x, n)
  right_x.r <- right_x.l + hist$density * xscale
  left_x.r <- rep(x, n)
  left_x.l <- left_x.r - hist$density * xscale
  x.l = c(left_x.l,right_x.l)
  x.r = c(left_x.r,right_x.r)
  
  # Define vertical limits of bars
  y.b <- hist$breaks[1:n]
  y.b = rep(y.b,2)
  y.t <- hist$breaks[2:(n + 1)]
  y.t = rep(y.t,2)
  
  rect(xleft = x.l, ybottom = y.b, xright = x.r, ytop = y.t,
       col = fillCol, border = fillCol)
}

# Main plotting function
######################################################################
makePlot = function(true.vector, retrieved.matrix, binWidth = 0.5,
                    ylim=NA, yaxis = NA, ylab="Parameter name", 
                    ylab.cex = 1.3, xlab = NA, xlab.cex = 1.2,
                    fillCol="indianred4", bg.color = "gray96",
                    hor.col1="gray60", hor.col2 = "black", 
                    output.Folders = NA, add.N = FALSE, cex.size = 1.2,
                    group.Names = NA, internal.margin.X = NA){
  data = retrieved.matrix
  numberOfHists <- ncol(data)
  binStarts <- 1:numberOfHists
  binMids <- binStarts + binWidth / 2
  axisCol <- "gray20"
  
  ## Data handling
  allValues <- unlist(as.list(data))
  DOYrange <- range( allValues, na.rm = TRUE )
  DOYrange <- c(floor(DOYrange[1]), ceiling(DOYrange[2]))
  DOYmean <- apply(data,2,mean)
  quantiles <- round(quantile( DOYrange, c(0.2, 0.4, 0.6, 0.8), na.rm=TRUE ), digits=0)
  gridlines <- round(quantile( DOYrange, c(0.1, 0.3, 0.5, 0.7, 0.9), na.rm=TRUE ), digits=0)
  ## Detect special cases
  multipleFolders <- sum(!is.na(output.Folders))>0
  noGroupLabelsProvided <- sum(is.na(group.Names))>0
  truevalueProvided <- sum(!is.na(true.vector))>0
  
  ## Get the histogram obects
  histList <- apply( data, 2, function(x, hCol) hist(x, plot = FALSE, breaks = 20))
  
  ## Set default values
  if(is.na(internal.margin.X)){ internal.margin.X <- 0.5}
  xlim <- c(1-internal.margin.X,numberOfHists+internal.margin.X)
  if(sum(is.na(ylim))>0){ylim <- DOYrange}
  if(is.na(xlab)){xlab <- ""}
  if(noGroupLabelsProvided){ 
    if(multipleFolders){
      group.Names <- rep(1:length(unique(true.vector)),
                         length(output.Folders))
    }else{
    group.Names <- paste(1:numberOfHists)
    }
  }
  
  ylim.Plot <- c(ylim[1],ylim[2]*1)
  plot(c(0, 5), DOYrange, type = "n", xlim=xlim+0.01, ylim=ylim.Plot,
       ann = FALSE, axes = FALSE, xaxs = "i", yaxs = "i")
  rect(xleft = xlim[1], ybottom = ylim.Plot[1], 
       xright = xlim[2], ytop = ylim.Plot[2],
       col = bg.color, border = "white")

  if(is.na(yaxis)){ 
      if(sum(is.na(true.vector))>0){
          yaxis = format(round(seq(ylim[1],ylim[2],length.out=8),digits = 2), nsmall = 1)
      }else{
          central <- sort(true.vector)
          lower <- ylim[1]
          upper <- ylim[2]
          yaxis = format(round(c(lower,central,upper),digits = 1), nsmall = 1)
      }
    }
  
  axis(1, at=1:numberOfHists, labels=rep("",numberOfHists),
       cex.axis = cex.size, col = axisCol, tck=-0.02)
  axis(1, at=1:numberOfHists, labels=group.Names,
       cex.axis = cex.size, tick = FALSE, line=-0.4)
  axis(2, at = yaxis, labels = rep("",length(yaxis)),
       cex.axis = cex.size, col = axisCol, tck = -.02)
  axis(2, at = yaxis, labels = yaxis, cex.axis = cex.size, 
       line = -0.45, tick=FALSE, las=2)
  mtext(side = 2, outer = F, line = 2.2, ylab, cex = ylab.cex)
  mtext(side = 1, outer = F, line = 2.2, xlab, cex = xlab.cex)
  
  #box(bty = "o", col = axisCol)
  if(truevalueProvided){
      nT <- length(true.vector)
      abline(h=true.vector, lty=3, col = hor.col1)
      for(t in 1:nT){
      tv <- true.vector[t]
      offset <- 0.5
      lines(c(t-offset,t+offset),c(tv,tv), lty=5, lwd = 1, col = hor.col2)
      }
      abline(v=1:nT+0.5, col="gray85")
  }
  
  if(sum(!is.na(output.Folders))>0){
      get.N <- as.numeric(gsub("\\D", "", output.Folders))
      nSim <- length(output.Folders)
      simSize <- length(true.vector)/nSim
      simSep <- seq(simSize, length(true.vector), simSize)+0.5
      for(sep in 1:length(simSep)){
          lines(c(simSep[sep],simSep[sep]),ylim, lty=1, col = "black")
      }
      if(add.N){
        text(simSep-1.5,ylim.Plot[2]+0.25,paste("n = ",get.N, sep=""),
             cex=xlab.cex*1.4, xpd=NA)
      }
  }
  
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