##########################################################
# LOAD FUNCTIONS/PACKAGES
##########################################################
######## Install/Load required R packages
    # Required packages
    list.of.packages <- c("R2jags", "magrittr", "rstan", "tictoc", "readr")
    # Install if not already installed
    new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
    if(length(new.packages) > 0){  install.packages(new.packages, dep=TRUE) }
    # Load packages
    for(package.i in list.of.packages){
        suppressPackageStartupMessages(library(package.i, character.only = TRUE))
    }

###### Load R scripts with functions required to generate data
    source("../Functions/generateRandomParameterValues.R")
    source("../Functions/simulateDataCDDM.R")
    source("../Functions/processJAGSsamples.R")

###### Load JAGS modules
    load.module("cddm")
    
##########################################################
# SIMULATION SETTINGS
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
possible.combinations <- s.topIdx * a.topIdx * m.topIdx * b.topIdx * n.topIdx

# Settings to be passed into JAGS
n.iter    = 2500
n.burnin  = 500
n.chains  = 4
n.thin    = 1

# Output setup
output.folder <- "./output_dec22/"
simstudy.Name <- "simStudy_"

################################################################
# MAIN FUNCTIONS: Generating data and Recovering parameters
################################################################

######## Function to GENERATE DATA
generate <- function(seed) {
        set.seed(seed)
        X <- cddm.simData(
                          trials = sampleSize,
                          drift.Angle = true.driftAngle,
                          drift.Length = true.driftLength,	
                          boundary = true.bound,
                          ndt = true.nondecision)$data
        output <- t(X)
        return(output)
}

######## Function to RECOVER parameter values
recover <- function(data) {
      ### Write JAGS model  #########################
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

      ### General sampling settings
      n.chains =  n.chains
      n.iter   =  n.iter
      n.burnin =  n.burnin
      n.thin   =  n.thin

      data <- list(X=data)
      parameters <- c("drift", "bound", "ter0", "theta0")

      ########## RUN MODEL ######################################################
      tic()  # Set timer
      samples <- jags.parallel(data = data, parameters.to.save =parameters,
                               model.file = modelFile, jags.module = 'cddm',
                               n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin=n.thin, DIC=T)
      elapsed.time <- toc()$callback_msg   # Record time
      clock <- parse_number(elapsed.time)  # Extract the seconds
      ###########################################################################
      # PROCESS SAMPLES
      ### Extract samples per parameter
      driftLength.samples  <- samples$BUGSoutput$sims.array[,,"drift"]
      bound.samples        <- samples$BUGSoutput$sims.array[,,"bound"]
      ndt.samples          <- samples$BUGSoutput$sims.array[,,"ter0"]
      theta0.samples       <- samples$BUGSoutput$sims.array[,,"theta0"]

      ### Compute mean posterior per parameter
      est.driftLength  <- mean( driftLength.samples )
      est.bound        <- mean( bound.samples       )
      est.nondecision  <- mean( ndt.samples         )
      est.driftAngle   <- mean( theta0.samples      )
      mean.est <- c(est.driftLength, est.bound, est.nondecision, est.driftAngle)

      # Compute MAP per parameter
      map.driftLength  <- myJAGS.MAP( driftLength.samples )
      map.bound        <- myJAGS.MAP( bound.samples       )
      map.nondecision  <- myJAGS.MAP( ndt.samples         )
      map.driftAngle   <- myJAGS.MAP( theta0.samples      )
      map.est <- c(map.driftLength, map.bound, map.nondecision, map.driftAngle)

      # Compute standard deviation
      sd.driftLength  <- sd( driftLength.samples )
      sd.bound        <- sd( bound.samples       )
      sd.nondecision  <- sd( ndt.samples         )
      sd.driftAngle   <- sd( theta0.samples      )
      sd.est <- c(sd.driftLength, sd.bound, sd.nondecision, sd.driftAngle)

      # Compute Rhats
      object <- samples$BUGSoutput$sims.array
      Rhats <- apply(object,3,Rhat)

      # Prepare recovery output
      output <- list(theta0.samples,
                     driftLength.samples,
                     bound.samples,
                     ndt.samples,
                     mean.est,
                     map.est,
                     sd.est,
                     Rhats,
                     clock)
      names(output) <- c("theta0.samples",
                         "driftLength.samples",
                         "bound.samples",
                         "ndt.samples",
                         "mean.estimates",
                         "map.estimates",
                         "std.estimates",
                         "Rhats",
                         "time.elapsed")
  
return(output)
}

###################################################################################
##  Run the simulation
###################################################################################

possible.combinations <- s.topIdx * a.topIdx * m.topIdx * b.topIdx * n.topIdx
  
################ Create empty arrays to store the simulation output
  par.labels <- c("driftLength","bound","ndt","driftAngle") # Shared labels
  # True parameter values used to generate data
  trueValues <- array(NA, dim=c(possible.combinations,5))
  colnames(trueValues) <- c("trials",par.labels)
  # Posterior chains
  theta0.samples <- array(NA,dim=c(2000,4,iterations,possible.combinations))
  ndt.samples <- array(NA,dim=c(2000,4,iterations,possible.combinations))
  driftL.samples <- array(NA,dim=c(2000,4,iterations,possible.combinations))
  bound.samples <- array(NA,dim=c(2000,4,iterations,possible.combinations))
  # Mean posteriors
  retrievedValues <- array(NA,dim=c(iterations,4,possible.combinations))
  colnames(retrievedValues) <- par.labels
  # Standard deviation
  retrievedValues_sd <- array(NA,dim=c(iterations,4,possible.combinations))
  colnames(retrievedValues_sd) <- par.labels
  # MAP
  mapValues <- array(NA,dim=c(iterations,4,possible.combinations))
  colnames(mapValues) <- par.labels
  # R hats
  rhats <- array(NA,dim=c(iterations,5,possible.combinations))
  # Seconds elapsed per simulation
  timers <- array(NA,dim=c(iterations,possible.combinations))
  
  ################ Run simulation
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
                       
                      this.truth <- c(sampleSize,
                                      true.driftLength,
                                      true.bound,
                                      true.nondecision,
                                      true.driftAngle)
                      trueValues[page,] <- this.truth
		                  
                      for(i in 1:iterations){
			                    X <- generate(seed = i)
			                    Y <- recover(data = X) 
			                    retrievedValues[i,,page] <- Y$mean.estimates
			                    retrievedValues_sd[i,,page] <- Y$std.estimates
			                    mapValues[i,,page] <- Y$map.estimates
			                    timers[i,page] <- Y$time.elapsed
			                    rhats[i,,page] <- Y$Rhats
			                    theta0.samples[,,i,page] <- Y$theta0.samples
			                    ndt.samples[,,i,page] <- Y$ndt.samples
			                    driftL.samples[,,i,page] <- Y$driftLength.samples
			                    bound.samples[,,i,page] <- Y$bound.samples
                      }
                      
                      this.matrix.means <- retrievedValues[,,page]
                      this.matrix.sd <- retrievedValues_sd[,,page]
                      this.matrix.MAPs <- mapValues[,,page]
                      this.matrix.rhats <- rhats[,,page]
                      this.timer <- timers[,page]
                      these.thetas <- theta0.samples[,,,page]
                      these.ndts <- ndt.samples[,,,page]
                      these.lengths <- driftL.samples[,,,page]
                      these.bounds <- bound.samples[,,,page]
                      
                      Z <- list(this.truth,
                                this.matrix.means,
                                this.matrix.sd,
                                this.matrix.MAPs,
                                this.matrix.rhats,
                                this.timer,
                                these.thetas,
                                these.ndts,
                                these.lengths,
                                these.bounds)
                      names(Z) <- c("current.truth",
                                    "current.matrix.means",
                                    "current.matrix.sd",
                                    "current.matrix.MAPs",
                                    "current.matrix.rhats",
                                    "current.timer",
                                    "current.thetas",
                                    "current.ndts",
                                    "current.lengths",
                                    "current.bounds")
                      
                      if(is.na(run.id)){ 
                                        parallel.run.id <- "_"
                      }else{
                            parallel.run.id <- paste("_",run.id,"_",sep="")
                      }
                       
                      fileName <- paste(output.folder,
                                        "subset",page,
                                        parallel.run.id,
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
save(rhats, file = paste(output.folder,simstudy.Name,"Rhats.RData",sep=""))
save(trueValues, file = paste(output.folder,simstudy.Name,"trueValues.RData",sep=""))
save(retrievedValues, file = paste(output.folder,simstudy.Name,"meanPosteriors.RData",sep=""))
save(retrievedValues_sd, file = paste(output.folder,simstudy.Name,"std.RData",sep=""))
save(mapValues, file = paste(output.folder,simstudy.Name,"MAPs.RData",sep=""))
save(timers, file = paste(output.folder,simstudy.Name,"timers.RData",sep=""))
save(theta0.samples, file = paste(output.folder,simstudy.Name,"theta0.RData",sep=""))
save(ndt.samples, file = paste(output.folder,simstudy.Name,"ndt.RData",sep=""))
save(driftL.samples, file = paste(output.folder,simstudy.Name,"driftLength.RData",sep=""))
save(bound.samples, file = paste(output.folder,simstudy.Name,"bound.RData",sep=""))


######################################################################################
#####  CONVERGENCE CHECKS
######################################################################################
hist(rhats)
abline(v=1.05, col="red", lty=2)
legend("top","Rhat = 1.05", lty=2, col="red", cex=0.4)

bad.Rhat <- which(rhats>1.04,arr.ind = TRUE)
colnames(rhats[,bad.Rhat[,2],])

######################################################################################
######            P L O T S                                  #########################
######################################################################################
