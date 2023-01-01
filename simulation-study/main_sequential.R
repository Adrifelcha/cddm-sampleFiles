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

# Determine whether simulation output is already available
test <- !file.exists(file=paste(output.folder,simstudy.Name,"meanPosteriors.RData",sep=""))

################################################################
# Define functions to GENERATE DATA and RECOVER PARAMETERS
################################################################
######  1. Function to GENERATE DATA 
####################################
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

######  2. Function to RECOVER PARAMETER VALUES
###############################################
recover <- function(data) {
      ########## Write JAGS model  #############################################
      modelFile <- "cddm_JAGSmodel.bug"
      write('
              data {
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
      samples <- jags.parallel(data = data, 
                               parameters.to.save = parameters,
                               model.file = modelFile, 
                               jags.module = 'cddm',
                               n.chains = n.chains, 
                               n.iter = n.iter, 
                               n.burnin = n.burnin, 
                               n.thin=n.thin, 
                               DIC=T)
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

      ### Compute MAP per parameter
      map.driftLength  <- myJAGS.MAP( driftLength.samples )
      map.bound        <- myJAGS.MAP( bound.samples       )
      map.nondecision  <- myJAGS.MAP( ndt.samples         )
      map.driftAngle   <- myJAGS.MAP( theta0.samples      )
      map.est <- c(map.driftLength, map.bound, map.nondecision, map.driftAngle)

      ### Compute standard deviation
      sd.driftLength  <- sd( driftLength.samples )
      sd.bound        <- sd( bound.samples       )
      sd.nondecision  <- sd( ndt.samples         )
      sd.driftAngle   <- sd( theta0.samples      )
      sd.est <- c(sd.driftLength, sd.bound, sd.nondecision, sd.driftAngle)

      ### Compute Rhats
      object <- samples$BUGSoutput$sims.array
      Rhats <- apply(object,3,Rhat)

      ### Prepare recovery output
      output <- list("mean.estimates" = mean.est,
                     "map.estimates"  = map.est,
                     "std.estimates"  = sd.est,
                     "Rhats" = Rhats,
                     "time.elapsed" = clock)
return(output)
}

###################################################################################
##  Run the simulation
###################################################################################
if(test){
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
        colnames(trueValues)  <- c(par.labels,"trials")
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
  
######## Run simulation
    page <- 1
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
                                                true.driftAngle,
                                                sampleSize)
                        trueValues[page,] <- this.truth
    		                  
                          for(i in 1:iterations){
    			                    X   <-  generate(seed = i)
    			                    Y   <-  recover(data = X) 
    			                    retrievedValues[i,,page]    <- Y$mean.estimates
    			                    retrievedValues_sd[i,,page] <- Y$std.estimates
    			                    mapValues[i,,page]          <- Y$map.estimates
    			                    timers[i,page]              <- Y$time.elapsed
    			                    rhats[i,,page]              <- Y$Rhats
    			                    theta0.samples[,,i,page]    <- Y$theta0.samples
    			                    ndt.samples[,,i,page]       <- Y$ndt.samples
    			                    driftL.samples[,,i,page]    <- Y$driftLength.samples
    			                    bound.samples[,,i,page]     <- Y$bound.samples
                          }
                          
                          this.matrix.means <- retrievedValues[,,page]
                          this.matrix.sd    <- retrievedValues_sd[,,page]
                          this.matrix.MAPs  <- mapValues[,,page]
                          this.matrix.rhats <- rhats[,,page]
                          this.timer        <- timers[,page]
                          these.thetas      <- theta0.samples[,,,page]
                          these.ndts        <- ndt.samples[,,,page]
                          these.lengths     <- driftL.samples[,,,page]
                          these.bounds      <- bound.samples[,,,page]
                          
                          Z <- list("current.truth" = this.truth,
                                    "current.matrix.means" = this.matrix.means,
                                    "current.matrix.sd" = this.matrix.sd,
                                    "current.matrix.MAPs" = this.matrix.MAPs,
                                    "current.matrix.rhats" = this.matrix.rhats,
                                    "current.timer" = this.timer)
                          fileName <- paste(output.folder,
                                            "subset",page,
                                            "_m",n,
                                            "b",a,
                                            "n",b,
                                            "a",m,
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
}else{
    load(paste(output.folder,simstudy.Name,"trueValues.RData",sep=""))
    load(paste(output.folder,simstudy.Name,"Rhats.RData",sep=""))
    load(paste(output.folder,simstudy.Name,"meanPosteriors.RData",sep=""))
    load(paste(output.folder,simstudy.Name,"std.RData",sep=""))
    load(paste(output.folder,simstudy.Name,"MAPs.RData",sep=""))
    load(paste(output.folder,simstudy.Name,"timers.RData",sep=""))
}

######################################################################################
#####  CONVERGENCE CHECKS
######################################################################################
rule <- 1.05
bad.Rhat <- which(rhats>rule,arr.ind = TRUE)

hist(rhats)
    abline(v=rule, col="red", lty=2)
    legend("top",paste("Rhat = ",rule," | ",
                       round(nrow(bad.Rhat)/prod(dim(rhats)),5),
                       "% of chains"), lty=2, col="red", cex=0.4)
    
bad.Rhat.ID <- colnames(rhats[,bad.Rhat[,2],])
table(bad.Rhat.ID)

######################################################################################
######            P L O T S                                  #########################
######################################################################################
boxplot.perPar <- function(parameter.name, color="blue"){
        par.levels <- table(trueValues[,parameter.name])
        par.values <- as.numeric(names(par.levels))
        group.by.level <- matrix(NA,nrow=max(par.levels)*200,
                                 ncol=length(par.values))
        colors <- paste(color,2:4,sep="")
        for(i in 1:length(par.values)){
              a <- par.values[i]
              same.level <- trueValues[,parameter.name]==a
              group.by.level[,i] <- retrievedValues[,parameter.name,same.level]
            }
        boxplot(group.by.level, col=colors, pch=16, cex=0.5)
        abline(h=par.values, col=colors, lty=2)
        mtext(paste("Parameter:",parameter.name),3)
}

boxplot.perPar("driftLength", "blue")
boxplot.perPar("ndt", "green")
boxplot.perPar("bound", "purple")
boxplot.perPar("driftAngle", "indianred")
