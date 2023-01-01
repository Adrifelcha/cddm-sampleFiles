##########################################################
# LOAD FUNCTIONS/PACKAGES
##########################################################
######## Load required R packages for parallel processing
  library(foreach)
  library(doParallel)

##########################################################
# SIMULATION SETTINGS
##########################################################
settings <- list(
              iterations = 200,
              sampleSize.list  =  c(100, 500, 2500),
              driftAngle.list  =  c(0.0, 2.0, 4.0),
              driftLength.list =  c(0.0, 1.0, 2.0),
              bound.list       =  c(1.5, 2.0, 2.5),
              nondecision.list =  c(0.1, 0.2, 0.3))
settings <- c(settings,
              s.topIdx  = length(settings$sampleSize.list),
              a.topIdx  = length(settings$driftAngle.list),
              m.topIdx  = length(settings$driftLength.list),   # 'm' for magnitude
              b.topIdx  = length(settings$bound.list),
              n.topIdx  = length(settings$nondecision.list))
settings <- c(settings,
              possible.combinations = settings$s.topIdx * settings$a.topIdx * settings$m.topIdx * settings$b.topIdx * settings$n.topIdx,
              n.iter    = 2500,
              n.burnin  = 500,
              n.chains  = 4,
              n.thin    = 1,
              output.folder = "./parallel_seeds/",
              simstudy.Name = "simStudy_")

################################################################
# Define simulation functions
################################################################
######  1. Function to GENERATE DATA 
####################################
generate <- function(seed, sampleSize, true.driftAngle, true.driftLength, true.bound, true.nondecision) {
    set.seed(seed)
    X <-  cddm.simData(trials = sampleSize,
                       drift.Angle = true.driftAngle,
                       drift.Length = true.driftLength,	
                       boundary = true.bound,
                       ndt = true.nondecision)$data
    output <- t(X)
    return(output)
}

######  2. Function to RECOVER PARAMETER VALUES
###############################################
recover <- function(data, settings) {
      ########## Write JAGS model  #############################################
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
      data <- list(X=data)
      parameters <- c("drift", "bound", "ter0", "theta0")
    
      ########## RUN MODEL ######################################################
      tic <- clock::date_now(zone="UTC")
      samples <- jags(data = data, 
                      parameters.to.save = parameters,
                      model = modelFile, 
                      n.chains = settings$n.chains, 
                      n.iter   = settings$n.iter, 
                      n.burnin = settings$n.burnin, 
                      n.thin   = settings$n.thin, 
                      DIC = T)
      toc <- clock::date_now(zone="UTC")
      clock <- as.numeric(toc-tic, units="secs")  # Record time
      
      ########## PROCESS SAMPLES ################################################
      ### Extract samples per parameter
      driftLength.samples  <- samples$BUGSoutput$sims.array[,,"drift"]
      bound.samples        <- samples$BUGSoutput$sims.array[,,"bound"]
      ndt.samples          <- samples$BUGSoutput$sims.array[,,"ter0"]
      theta0.samples       <- samples$BUGSoutput$sims.array[,,"theta0"]
      theta0.samples       <- theta0.samples %% (2*pi)
      
      ### Compute circular mean
      mu1 <- array(NA, dim=dim(theta0.samples))  # Create empty vectors
      mu2 <- array(NA, dim=dim(theta0.samples))
      for(i in 1:n.chains){
          mu1[,i] <- cddm.polarToRect(theta0.samples[,i],driftLength.samples[,i])$mu1
          mu2[,i] <- cddm.polarToRect(theta0.samples[,i],driftLength.samples[,i])$mu2
          }
      Mu1 <- mean(mu1)
      Mu2 <- mean(mu2)
      est.driftAngle   <- cddm.getVectorAngle(Mu1,Mu2)
      ### Compute mean posterior per parameter
      est.driftLength  <- mean( driftLength.samples )
      est.bound        <- mean( bound.samples       )
      est.nondecision  <- mean( ndt.samples         )
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
      output <- list("mean.estimates" = mean.est, 
                     "map.estimates"  = map.est,
                     "std.estimates"  = sd.est,
                     "Rhats" = Rhats,
                     "time.elapsed" = clock)
return(output)
}

######  3. Function to RUN the SIMULATION
#########################################
run_simulation <- function(seed=i, settings){
      ########## Load packages/scripts ##########################################
      ### R packages
      library(R2jags)
      library(magrittr)
      library(rstan)
      library(readr)
      ### R scripts with custom R functions
      source("../Functions/generateRandomParameterValues.R")
      source("../Functions/simulateDataCDDM.R")
      source("../Functions/processJAGSsamples.R")
      ### JAGS CDDM module
      load.module("cddm")

      ########## Load relevant settings ##########################################
      # List of values to be used in simulation study
      sampleSize.list   <-  settings$sampleSize.list
      driftAngle.list   <-  settings$driftAngle.list
      driftLength.list  <-  settings$driftLength.list
      bound.list        <-  settings$bound.list
      nondecision.list  <-  settings$nondecision.list
      # Identify maximum valid index for every list
      s.topIdx  <-  settings$s.topIdx
      a.topIdx  <-  settings$a.topIdx
      m.topIdx  <-  settings$m.topIdx
      b.topIdx  <-  settings$b.topIdx
      n.topIdx  <-  settings$n.topIdx
      possible.combinations <- settings$possible.combinations
      
      ########## Run simulation for every Parameter combination #################
      out <- list()
      for(s in 1:s.topIdx){
          for(m in 1:m.topIdx){
              for(b in 1:b.topIdx){
                  for(n in 1:n.topIdx){
                      for(a in 1:a.topIdx){
                            sampleSize       <-   sampleSize.list  [s]
                            true.driftAngle  <-   driftAngle.list  [a] 
                            true.driftLength <-   driftLength.list [m]
                            true.bound       <-   bound.list       [b]
                            true.nondecision <-   nondecision.list [n]
                
                            this.truth        <-  c(true.driftLength,
                                                    true.bound,
                                                    true.nondecision,
                                                    true.driftAngle)
                            
                            X   <-  generate(seed = i,
                                           sampleSize = sampleSize,
                                           true.driftAngle = true.driftAngle,
                                           true.driftLength = true.driftLength,
                                           true.bound = true.bound,
                                           true.nondecision = true.nondecision)
                            Y   <-  recover(data = X, settings)
                            out <- rbind(out,list(seed = i,
                                                  n = sampleSize,
                                                  true.values    = this.truth,
                                                  mean.estimates = Y$mean.estimates,
                                                  std.estimates  = Y$std.estimates,
                                                  map.estimates  = Y$map.estimates,
                                                  elapsed.time   = Y$time.elapsed,
                                                  rhats = Y$Rhats))
                      }
                  }
              }
          }
      }
return(out)      
}

######  3. Function to STORE the OUTPUT
#######################################
store_output <- function(output, settings){
      ########## Load relevant Settings #########################################
      ### R packages
      iterations            <- settings$iterations
      possible.combinations <- settings$possible.combinations
      output.folder         <- settings$output.folder
      studyName             <- settings$studyName
      n.iter    <-  settings$n.iter
      n.burnin  <-  settings$n.burnin
      
      ########## Create empty arrays to save output #############################
      ### Size variables
      # Parameter labels
      par.labels <- c("driftLength","bound","ndt","driftAngle")
      # No. of parameters
      npar <- length(par.labels)
      # Extensive no. of columns (for True value matrix and Rhats)
      ncols <- npar+1
      # Number of samples kept per chain
      nrows <- n.iter-n.burnin
      ### Array 1: True parameter values used to generate data
      trueValues            <- array(NA, dim=c(possible.combinations,npar))
      colnames(trueValues)  <- par.labels
      ### Array 2: Mean posteriors
      retrievedValues           <- array(NA,dim=c(iterations,npar,possible.combinations))
      colnames(retrievedValues) <- par.labels
      ### Array 3: Standard deviation
      retrievedValues_sd            <- array(NA,dim=c(iterations,npar,possible.combinations))
      colnames(retrievedValues_sd)  <- par.labels
      ### Array 4: MAP
      mapValues           <- array(NA,dim=c(iterations,npar,possible.combinations))
      colnames(mapValues) <- par.labels
      ### Array 5: R hats
      rhats   <- array(NA,dim=c(iterations,ncols,possible.combinations))
      ### Array 6: Seconds elapsed per simulation
      timers  <- array(NA,dim=c(iterations,possible.combinations))
      ### Array 7: Record seeds
      seeds <- array(NA,dim=c(iterations,possible.combinations))
      ### Array 8: Record sample sizes
      size <- array(NA,dim=c(iterations,possible.combinations))
      
      out.size <- possible.combinations * iterations
      for(set in 1:possible.combinations){
          J <- seq(set,out.size,possible.combinations) 
          for(i in 1:iterations){
              j <- J[i]
              S <- output[j,]
              seeds[i,set] <- S$seed
              size[i,set] <- S$n
              timers[i,set] <- S$elapsed.time
              rhats[i,,set] <- S$rhats
              trueValues[set,] <- S$true.values
              retrievedValues[i,,set] <- S$mean.estimates
              retrievedValues_sd[i,,set] <- S$std.estimates
              mapValues[i,,set] <-S$map.estimates
          }
      }
      
      save(timers, file = paste(output.folder,studyName,"_timers.RData",sep=""))
      save(seeds, file = paste(output.folder,studyName,"_seeds.RData",sep=""))
      save(size, file = paste(output.folder,studyName,"_sizes.RData",sep=""))
      colnames(rhats) <- names(S$Rhats)
      save(rhats, file = paste(output.folder,simstudy.Name,"Rhats.RData",sep=""))
      save(trueValues, file = paste(output.folder,simstudy.Name,"trueValues.RData",sep=""))
      save(retrievedValues, file = paste(output.folder,simstudy.Name,"meanPosteriors.RData",sep=""))
      save(retrievedValues_sd, file = paste(output.folder,simstudy.Name,"std.RData",sep=""))
      save(mapValues, file = paste(output.folder,simstudy.Name,"MAPs.RData",sep=""))
}

################################################################
# Define simulation functions
################################################################
cores       <-  detectCores()
my.cluster  <-  makeCluster(cores[1]-6)

registerDoParallel(cl = my.cluster)
output <- foreach(i = 1:settings$iterations, 
                  .errorhandling = "pass",
                  .combine = 'rbind'
                  ) %dopar% {
                  Z <- run.Sim(seed = i, settings)
                  }
stopCluster(cl = my.cluster)

store_output(output,settings)

######################################################################################
#####  CONVERGENCE CHECKS
######################################################################################
hist(rhats)
abline(v=1.05, col="red", lty=2)
legend("top","Rhat = 1.05", lty=2, col="red", cex=0.4)

bad.Rhat <- which(rhats>1.04,arr.ind = TRUE)
colnames(rhats[,bad.Rhat[,2],])
