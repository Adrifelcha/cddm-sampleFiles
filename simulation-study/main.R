# Load required R libraries
library(R2jags)
library(magrittr)
library(rstan)
library(readr)

# Load Rscripts with functions required to generate data
source("../Functions/generateRandomParameterValues.R")
source("../Functions/simulateDataCDDM.R")

# Load JAGS modules
load.module("cddm")
#load.module("vonmises")

##########################################################
# Simulation settings
##########################################################
iterations <- 200

# List of values to be used in simulation study
sampleSize.list   <-  c(100, 500, 2500)
driftAngle.list   <-  c(0.0, 2.0, 4.0)
driftLength.list  <-  c(0.0, 1.0, 2.0)
bound.list        <-  c(1.5, 2.0, 2.5)
nondecision.list  <-  c(0.1, 0.2, 0.3)

# Identify maximum valid index for every list
s.topIdx  <-  length(sampleSize.list)
a.topIdx  <-  length(driftAngle.list)
m.topIdx  <-  length(driftLength.list)   # 'm' for magnitude
b.topIdx  <-  length(bound.list)
n.topIdx  <-  length(nondecision.list)

# Count the number of possible combinations
possible.combinations <- s.topIdx * a.topIdx * m.topIdx * b.topIdx * n.topIdx

################################################################
# Define functions for Generating data and Recovering parameters
################################################################
generate <- function(seed) {
  set.seed(seed)
  
  X <-  cddm.simData(
			sampleSize,
			true.driftAngle,
			true.driftLength,	
			true.bound,
			true.nondecision)$data
  
  # if(any(is.na(X[,2]))){
  #   na.found <- which(is.na(X[,2]))
  #   validValues <- X[-na.found,]
  #   randomlyChosen <- sample(1:nrow(validValues),length(na.found), replace=FALSE)
  #   X[na.found,] <- X[randomlyChosen,]
  #   print(paste("NAs found: ",length(na.found),sep=""))
  #   print(paste("Seed: ", seed,
  #               " | Angle: "      , true.driftAngle,
  #               " | Length: "     , true.driftLength,
  #               " | Bound: "      , true.bound, 
  #               " | NDT: "        , true.nondecision, 
  #               " | Sample size: ", sampleSize, 
  #                                   sep=""))
  # }
  
  output <- t(X)
  return(output)
}

recover <- function(data) {
  # Write JAGS model 
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
              drift  ~ dnorm(0, 1)T(0,)
              theta0 ~ dnorm(0, 0.1)
              bound  ~ dgamma(3, 2)
              ter0   ~ dexp(1)T(, tmin)
            }',
        modelFile)
  # General sampling settings
  n.chains =   4
  n.iter   = 2500 
  n.burnin = 500
  n.thin   =   1
  
  data <- list(X=data)
  parameters <- c("drift", "bound", "ter0", "theta0")

  # Run model
  
  tic()
  
  samples <- jags(data=data, parameters.to.save=parameters, model=modelFile, 
                  n.chains=n.chains, n.iter=n.iter, n.burnin=n.burnin, 
                  n.thin=n.thin, DIC=T)
  
  elapsed.time <- toc()$callback_msg
  clock <- parse_number(elapsed.time)
  
  # Compute mean posterior  
  est.driftLength  <- mean( samples$BUGSoutput$sims.array[,,"drift"]  )
  est.bound        <- mean( samples$BUGSoutput$sims.array[,,"bound"]  )
  est.nondecision  <- mean( samples$BUGSoutput$sims.array[,,"ter0"]   )
  est.driftAngle   <- mean( samples$BUGSoutput$sims.array[,,"theta0"] )
  mean.est <- c(est.driftLength, est.bound, est.nondecision, est.driftAngle)
  
  
  # Compute standard deviation
  sd.driftLength  <- sd( samples$BUGSoutput$sims.array[,,"drift"]  )
  sd.bound        <- sd( samples$BUGSoutput$sims.array[,,"bound"]  )
  sd.nondecision  <- sd( samples$BUGSoutput$sims.array[,,"ter0"]   )
  sd.driftAngle   <- sd( samples$BUGSoutput$sims.array[,,"theta0"] )
  sd.est <- c(sd.driftLength, sd.bound, sd.nondecision, sd.driftAngle)
  
  # Keep whole theta0 chains
  theta0.samples <- samples$BUGSoutput$sims.array[,,"theta0"]
  
  # Compute Rhats
  object <- samples$BUGSoutput$sims.array
  Rhats <- apply(object,3,Rhat)
  
  output <- list(theta0.samples,
                 mean.est, 
                 sd.est, 
                 Rhats,
                 clock)
  
  names(output) <- c("theta0.samples",
                     "mean.estimates", 
                     "std.estimates", 
                     "Rhats",
                     "time.elapsed")
return(output)
}

###################################################################################
# A function to run the simulation study
###################################################################################
run_sim_study <-function(){
  par.labels <- c("driftLength","bound","ndt","driftAngle")
  
  trueValues <- array(NA, dim=c(possible.combinations,5))
  colnames(trueValues) <- c("trials",par.labels)
  
  theta0.samples <- array(NA,dim=c(2000,4,iterations,possible.combinations))
  
  retrievedValues <- array(NA,dim=c(iterations,4,possible.combinations))
  colnames(retrievedValues) <- par.labels
  
  retrievedValues_sd <- array(NA,dim=c(iterations,4,possible.combinations))
  colnames(retrievedValues_sd) <- par.labels
  
  rhats <- array(NA,dim=c(iterations,5,possible.combinations))
  timers <- array(NA,dim=c(iterations,possible.combinations))

  page <- 1
  for(n in 1:n.topIdx){
      for(a in 1:a.topIdx){
          for(b in 1:b.topIdx){
              for(m in 1:m.topIdx){
                  for(s in 1:s.topIdx){
                      sampleSize       <-   sampleSize.list  [s]
	                    true.driftAngle  <-   driftAngle.list  [a] 
		                  true.driftLength <-   driftLength.list [m]
		                  true.bound       <-   bound.list       [b]
		                  true.nondecision <-   nondecision.list [n]
                       
                      this.truth <- c(sampleSize,true.driftLength,true.bound,true.nondecision,true.driftAngle)
                      trueValues[page,] <- this.truth
		                  
                      for(i in 1:iterations){
			                    X <- generate(seed = i)
			                    Y <- recover(data = X) 
			                    retrievedValues[i,,page] <- Y$mean.estimates
			                    retrievedValues_sd[i,,page] <- Y$std.estimates
			                    timers[i,page] <- Y$time.elapsed
			                    rhats[i,,page] <- Y$Rhats
			                    theta0.samples[,,i,page] <- Y$theta0.samples
                      }
                      
                      this.matrix.means <- retrievedValues[,,page]
                      this.matrix.sd <- retrievedValues_sd[,,page]
                      this.matrix.rhats <- rhats[,,page]
                      this.timer <- timers[,page]
                      these.thetas <- theta0.samples[,,,page]
                      
                      Z <- list(this.truth,
                                this.matrix.means,
                                this.matrix.sd,
                                this.matrix.rhats,
                                this.timer,
                                these.thetas)
                      
                      names(Z) <- c("current.truth",
                                    "current.matrix.means",
                                    "current.matrix.sd",
                                    "current.matrix.rhats",
                                    "current.timer",
                                    "current.thetas")
                       
                      fileName <- paste("./output_new/subset",page,"_",
                                         "n",n,"-",
                                         "a",a,"-",
                                         "b",b,"-",
                                         "m",m,"-",
                                         "s",s,".RData",sep="")
                      save(Z, file=fileName)
                       
                      page <- page+1
		                  }
	               }
	           }
        }
    }

  colnames(rhats) <- names(Y$Rhats)
  save(rhats,file="./output2/simStudy_Rhats.RData")
  
  save(trueValues,file="./output2/simStudy_trueValues.RData")
  save(retrievedValues, file="./output2/simStudy_meanPosteriors.RData")
  save(retrievedValues_sd, file="./output2/simStudy_std.RData")
  save(timers, file="./output2/simStudy_timers.RData")
}


##############################################################
# RUN SIMULATION STUDY
##############################################################
run_sim_study()