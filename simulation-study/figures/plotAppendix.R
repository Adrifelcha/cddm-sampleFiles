# Process posterior samples
######################################################################
source("./plotFunctions.R")

xlab <- "Levels"
legend.size <- 1.2
ylab.size <- 0.9
xlab.size <- ylab.size
legend.col <- "gray75"
#########################################################
# Appendix 1: CARTESIAN 
#########################################################
output.Folders <- c("../output_cart_n80/",
                    "../output_cart_n200/",
                    "../output_cart_n500/")

mu1 <- group.by.level.across.folders(output.Folders,"mu1")
mu2 <- group.by.level.across.folders(output.Folders,"mu2")
ndt <- group.by.level.across.folders(output.Folders,"ndt")
bound <- group.by.level.across.folders(output.Folders,"bound")

par(mfrow=c(4,1))
par(oma=c(0.1,0,3.5,0));  
par(mar = c(3.5, 4, 0, 0.5))

makePlot(true.vector = ndt$par.values, add.N = TRUE,
         retrieved.matrix=ndt$groups, binWidth = 0.0055, 
         ylab=expression(paste("Nondecision time (",tau,")")),
         ylim=c(-0.01,1), fillCol = "#3B77A9", bg.color="#E2EBF2", 
         output.Folders = output.Folders, cex.size = 1.1,
         hor.col1 = legend.col, ylab.cex = ylab.size)
legend(7.5,1.075, "True values   ", lty = 1, cex = legend.size, 
       col="black", bty = "n")

makePlot(true.vector = bound$par.values,  
         retrieved.matrix=bound$groups, binWidth = 0.1, 
         ylab=expression(paste("Boundary radius (",eta,")")),
         ylim=c(1,3.5), fillCol = "#3BA9A7", bg.color="#E4ECEC",
         hor.col1 = legend.col,  cex.size = 1.1,
         ylab.cex = ylab.size, output.Folders = output.Folders)

makePlot(true.vector = mu1$par.values,  
         retrieved.matrix=mu1$groups, binWidth = 0.04, 
         ylab=expression(paste("Step size on X (",mu[x],")")),
         ylim=c(-1,1.5), fillCol = "#653BA9", bg.color="#E7E4EC",
         output.Folders = output.Folders, cex.size = 1.1,
         hor.col1 = legend.col, ylab.cex = ylab.size)

makePlot(true.vector = mu2$par.values,  
         retrieved.matrix=mu2$groups, binWidth = 0.03, 
         ylab=expression(paste("Step size on Y (",mu[y],")")),
         ylim=c(-1.5,1), fillCol = "#3B3BA9", bg.color = "#E3E3ED",
         output.Folders = output.Folders, cex.size = 1.1,
         xlab.cex = xlab.size, xlab=xlab, 
         hor.col1 = legend.col, ylab.cex = ylab.size)


# #########################################################
# Appendix 2: POLAR
# #########################################################
output.Folders <- c("../output_polar_n80/",
                    "../output_polar_n200/",
                    "../output_polar_n500/")

driftAngle <- group.by.level.across.folders(output.Folders,"driftAngle")
driftLength <- group.by.level.across.folders(output.Folders,"driftLength")
ndt <- group.by.level.across.folders(output.Folders,"ndt")
bound <- group.by.level.across.folders(output.Folders,"bound")

par(mfrow=c(4,1))
par(mar = c(3.5, 4, 0, 0.5))

makePlot(true.vector = ndt$par.values, output.Folders = output.Folders,
         retrieved.matrix=ndt$groups, binWidth = 0.0025, bg.color="#E2EBF2",
         ylab=expression(paste("Nondecision time (",tau,")")), add.N = TRUE,
         ylim=c(-0.01,1), fillCol = "#3B77A9", ylab.cex = ylab.size,
         hor.col1 = legend.col, cex.size = 1.1)
legend(7.5,1.075, "True values   ", lty = 2, cex = legend.size, 
       col="black", bty = "n")

makePlot(true.vector = bound$par.values, ylab.cex = ylab.size,
         retrieved.matrix=bound$groups, binWidth = 0.15, 
         ylab=expression(paste("Boundary radius (",eta,")")),
         ylim=c(1,3.5), fillCol = "#3BA9A7", cex.size = 1.1,
         bg.color="#E4ECEC", hor.col1 = legend.col,  
         output.Folders = output.Folders)

makePlot(true.vector = driftLength$par.values, ylab.cex = ylab.size,
         retrieved.matrix=driftLength$groups, binWidth = 0.1,
         ylab=expression(paste("Drift length (",delta,")")),
         ylim=c(0,3), fillCol = "#E5B92D", bg.color = "#EDE8DE", 
         output.Folders = output.Folders, cex.size = 1.1,
         hor.col1 = legend.col)

angle.is.4 <- which(driftAngle$par.values==4)
driftAngle$par.values[angle.is.4] <- driftAngle$par.values[angle.is.4] - 2*pi
makePlot(true.vector = driftAngle$par.values, ylab.cex = ylab.size,
         retrieved.matrix=driftAngle$groups, binWidth = 0.06, cex.size = 1.1,
         ylab=expression(paste("Drift angle (",theta,")")), xlab=xlab,
         ylim=c(-3.2,3), fillCol = "#E1751C", xlab.cex = xlab.size,
         output.Folders = output.Folders, bg.color = "#E9E3DE",
         hor.col1 = legend.col)