##########################################################
# LOAD FUNCTIONS/PACKAGES
##########################################################
######## Install/Load required R packages
    # Required packages
    list.of.packages <- c("R2jags", "magrittr", "rstan", "doParallel", 
                          "tictoc", "readr", "foreach")
    # Install if not already installed
    new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
    if(length(new.packages) > 0){
       install.packages(new.packages, dep=TRUE) 
      }
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
output.folder <- "./parallel_seeds/"
simstudy.Name <- "simStudy_"

################################################################
# MAIN FUNCTIONS: Generating data and Recovering parameters
################################################################

######## Function to GENERATE DATA
generate <- function(seed) {
        set.seed(seed)
        X <-  cddm.simData(
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
              data{
                    tmin <- 0.95 * min(X[2,])
                    D <- dim(X)
                    N <- D[2]
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
      samples <- jags(data = data, parameters.to.save = parameters, 
                      model = modelFile, n.chains = n.chains, n.iter = n.iter, 
                      n.burnin = n.burnin, n.thin = n.thin, DIC = T)
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
######## Define size variables
    # Parameter labels
      par.labels <- c("driftLength","bound","ndt","driftAngle")
    # No. of parameters
      npar <- length(par.labels)
    # Extensive no. of columns (for True value matrix and Rhats)
      ncols <- npar+1
    # Number of samples kept per chain
      nrows <- n.iter-n.burnin
    # Total number of samples per parameter
      total.samples <- nrows*n.chains

######## Create empty arrays to store the simulation output 
    # True parameter values used to generate data
      trueValues            <- array(NA, dim=c(possible.combinations,ncols))
      colnames(trueValues)  <- c("trials",par.labels)
    # Posterior chains
      posterior.chains  <- array(NA,dim=c(nrows,npar,iterations,possible.combinations))
      theta0.samples    <- posterior.chains
      ndt.samples       <- posterior.chains
      driftL.samples    <- posterior.chains
      bound.samples     <- posterior.chains
    # Mean posteriors
      retrievedValues           <- array(NA,dim=c(iterations,npar,possible.combinations))
      colnames(retrievedValues) <- par.labels
    # Standard deviation
      retrievedValues_sd            <- array(NA,dim=c(iterations,npar,possible.combinations))
      colnames(retrievedValues_sd)  <- par.labels
    # MAP
      mapValues           <- array(NA,dim=c(iterations,npar,possible.combinations))
      colnames(mapValues) <- par.labels
    # R hats
      rhats   <- array(NA,dim=c(iterations,ncols,possible.combinations))
    # Seconds elapsed per simulation
      timers  <- array(NA,dim=c(iterations,possible.combinations))

######## Locate parameter samples in parallel output
    theta.1st <- 19
    theta.last <- (theta.1st-1)+total.samples
    ndt.1st <- 1+theta.last
    ndt.last <-(ndt.1st-1)+total.samples
    length.1st <- ndt.last+1
    length.last <- (length.1st-1)+total.samples
    bound.1st <- length.last+1
    bound.last <- (bound.1st-1)+total.samples

######## Define clusters to run in parallel
    my.cluster <- parallel::makeCluster(
                          4, 
                          type = "FORK"
                  )
    doParallel::registerDoParallel(cl = my.cluster)
    
######## Run simulation
page <- 1
for(m in 1:m.topIdx){
    for(b in 1:b.topIdx){
        for(n in 1:n.topIdx){
            for(a in 1:a.topIdx){
                for(s in 1:s.topIdx){
                    sampleSize       <-   sampleSize.list  [s]
                    true.driftAngle  <-   driftAngle.list  [a] 
	                  true.driftLength <-   driftLength.list [m]
	                  true.bound       <-   bound.list       [b]
	                  true.nondecision <-   nondecision.list [n]
                     
                    this.truth        <-  c(sampleSize,
                                            true.driftLength,
                                            true.bound,
                                            true.nondecision,
                                            true.driftAngle)
                    trueValues[page,] <-  this.truth
	                  
                    output <- foreach(i = 1:iterations, .combine = 'rbind') %dopar% {
                                      X <-  generate(seed = i)
                                      Y <-  recover(data = X)
                                      output <- c(Y$mean.estimates,
                                                  Y$std.estimates,
                                                  Y$map.estimates,
                                                  Y$time.elapsed,
                                                  Y$Rhats,
                                                  Y$theta0.samples,
                                                  Y$ndt.samples,
                                                  Y$driftLength.samples,
                                                  Y$bound.samples)
                              }
                    
                    this.matrix.means <- output[,1:4]
                    this.matrix.sd    <- output[,5:8]
                    this.matrix.MAPs  <- output[,9:12]
                    this.timer        <- output[,13]
                    this.matrix.rhats <- output[,14:18]
                    these.thetas      <- matrix(output[,theta.1st:theta.last],byrow=TRUE,nrow=nrows,ncol=n.chains)
                    these.ndts        <- matrix(output[,ndt.1st:ndt.last],byrow=TRUE,nrow=nrows,ncol=n.chains)
                    these.lengths     <- matrix(output[,length.1st:length.last],byrow=TRUE,nrow=nrows,ncol=n.chains)
                    these.bounds      <- matrix(output[,bound.1st:bound.last],byrow=TRUE,nrow=nrows,ncol=n.chains)
                    
                        Z     <- list(this.truth,
                                        this.matrix.means,
                                        this.matrix.sd,
                                        this.matrix.MAPs,
                                        this.matrix.rhats,
                                        this.timer,
                                        these.thetas,
                                        these.ndts,
                                        these.lengths,
                                        these.bounds)
                  names(Z)    <- c("current.truth",
                                   "current.matrix.means",
                                   "current.matrix.sd",
                                   "current.matrix.MAPs",
                                   "current.matrix.rhats",
                                   "current.timer",
                                   "current.thetas",
                                   "current.ndts",
                                   "current.lengths",
                                   "current.bounds")
                    
                    fileName  <- paste(output.folder,
                                       "set",
                                       page,
                                       "_",
                                       "n",n,
                                       "a",a,
                                       "b",b,
                                       "m",m,
                                       "s",s,".RData",sep="")
                    save(Z, file = fileName)
                    page <- page+1
                }
            }
        }
    }
}

parallel::stopCluster(cl = my.cluster)
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
