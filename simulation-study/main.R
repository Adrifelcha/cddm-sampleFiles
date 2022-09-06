# Load Rscripts with functions required to generate data
library(R2jags)
library(magrittr)

load.module("cddm")
load.module("vonmises")

source("../Functions/generateRandomParameterValues.R")
source("../Functions/simulateDataCDDM.R")



# Lists of settings

iterations <- 200

sampleSize.list   <-  c(100, 500, 2500)
driftAngle.list   <-  c(0.0, 1.0, 2.0)
driftLength.list  <-  c(1.0, 2.0, 3.0)
bound.list        <-  c(1.0, 2.0, 3.0)
nondecision.list  <-  c(0.1, 0.2, 0.3)

s.Idx  <-  1
a.Idx  <-  1
l.Idx  <-  1
b.Idx  <-  1
n.Idx  <-  1

sampleSize        <-  sampleSize.list  [s.Idx]
true.driftAngle   <-  driftAngle.list  [a.Idx]
true.driftLength  <-  driftLength.list [l.Idx]
true.bound        <-  bound.list       [b.Idx]
true.nondecision  <-  nondecision.list [n.Idx]

generate <- function() { 
  cddm.simData(
    sampleSize,
    true.driftAngle,
    true.driftLength,
    true.bound,
    true.nondecision) %>% t()
  }


recover <- function(data) { 
  modelFile <- "cddm.bug"
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
              drift  ~ dnorm(0, 0.1)T(-3.14, 3.14)
              theta0 ~ dnorm(0, 0.1)T(-3.14, 3.14)
              bound  ~ dgamma(2, 2)
              ter0   ~ dexp(1)T(, tmin)
            }',
        modelFile)
  
  n.chains =   4
  n.iter   = 900 
  n.burnin = 400
  n.thin   =   1
  
  data <- list(X=data)
  parameters <- c("drift", "bound", "ter0", "theta0")
  
  fileName <- "samples.RData"
  samples <- jags(data=data, parameters.to.save=parameters, model=modelFile, 
                  n.chains=n.chains, n.iter=n.iter, n.burnin=n.burnin, 
                  n.thin=n.thin, DIC=T)
  
  est.driftLength  <- samples$BUGSoutput$sims.array[,,"drift"]  %>% mean()
  est.bound        <- samples$BUGSoutput$sims.array[,,"bound"]  %>% mean()
  est.nondecision  <- samples$BUGSoutput$sims.array[,,"ter0"]   %>% mean()
  est.driftAngle   <- samples$BUGSoutput$sims.array[,,"theta0"] %>% mean()
  
  c(est.driftLength, est.bound, est.nondecision, est.driftAngle)
}


iterations <- 200
Y <- matrix(nrow = 4, ncol = iterations)
for (k in 1:iterations) { 
  X <- generate()
  Y[,k] <- recover(X)
}
Y
rowMeans(Y)
