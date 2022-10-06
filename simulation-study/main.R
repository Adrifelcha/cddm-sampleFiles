# Load required R libraries
library(R2jags)
library(magrittr)

# Load Rscripts with functions required to generate data
source("../Functions/generateRandomParameterValues.R")
source("../Functions/simulateDataCDDM.R")

# Load JAGS modules
load.module("cddm")
load.module("vonmises")

##########################################################
# Lists of settings
##########################################################
iterations <- 200

# List of values to be used in simulation study
sampleSize.list   <-  c(100, 500, 2500)
driftAngle.list   <-  c(0.0, 1.0, 2.0)
driftLength.list  <-  c(0.5, 1.0, 1.5)
bound.list        <-  c(1.5, 2.0, 2.5)
nondecision.list  <-  c(0.1, 0.2, 0.3)

# Identify maximum valid index for every list
s.topIdx  <-  length(sampleSize.list)
a.topIdx  <-  length(driftAngle.list)
l.topIdx  <-  length(driftLength.list)
b.topIdx  <-  length(bound.list)
n.topIdx  <-  length(nondecision.list)

################################################################
# Define functions for Generating data and Recovering parameters
################################################################
generate <- function() { 
  X <-  cddm.simData(
			sampleSize,
			true.driftAngle,
			true.driftLength,	
			true.bound,
			true.nondecision)
  output <- t(X)
  return(output)
}

recover <- function(data) { 
  modelFile <- "cddm_JAGSmodel.bug"
  write('
      data {
            tmin <- 0.95 * min(X[2,])
            D <- dim(X)
            N <- D[1]
      }
      model{
            # Likelihood
              for (i in 1:N) {
                   X[1:2,i] ~ dcddm(drift, bound, ter0, theta0)
              }
                
            # Priors
              drift  ~ dnorm(0, 0.1)T(0,)
              theta0 ~ dnorm(0, 0.1)
              bound  ~ dgamma(2, 2)
              ter0   ~ dexp(1)T(, tmin)
            }',
        modelFile)
  
  n.chains =   4
  n.iter   = 1400 
  n.burnin = 400
  n.thin   =   1
  
  data <- list(X=data)
  parameters <- c("drift", "bound", "ter0", "theta0")
  
  samples <- jags(data=data, parameters.to.save=parameters, model=modelFile, 
                  n.chains=n.chains, n.iter=n.iter, n.burnin=n.burnin, 
                  n.thin=n.thin, DIC=T)
  
  est.driftLength  <- mean( samples$BUGSoutput$sims.array[,,"drift"]  )
  est.bound        <- mean( samples$BUGSoutput$sims.array[,,"bound"]  )
  est.nondecision  <- mean( samples$BUGSoutput$sims.array[,,"ter0"]   )
  est.driftAngle   <- mean( samples$BUGSoutput$sims.array[,,"theta0"] %% (2*pi) )
  
  output <- c(est.driftLength, est.bound, est.nondecision, est.driftAngle)
return(output)
}

### TESTING generate() and recover() functions over first set of true values.
#############################################################################
sampleSize        <-  sampleSize.list  [1]
true.driftAngle   <-  driftAngle.list  [1]
true.driftLength  <-  driftLength.list [1]
true.bound        <-  bound.list       [1]
true.nondecision  <-  nondecision.list [1]

Y <- matrix(nrow = 4, ncol = iterations)
for (k in 1:iterations) { 
  X <- generate()
  Y[,k] <- recover(X)
}

superMean <- rowMeans(Y)
trueVals <- c(true.driftLength, true.bound, true.nondecision, true.driftAngle)
X <- rbind(trueVals,superMean)
colnames(X) <- c("length","bound","ndt","angle")
rownames(X) <- c("true","retrieved")
X


run_sim_study <-function(){

for(n in 1:s.topIdx){
    for(){

}
}


}

