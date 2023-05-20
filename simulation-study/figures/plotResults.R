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
