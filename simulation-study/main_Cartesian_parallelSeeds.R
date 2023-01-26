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
              iterations = 12,
              sampleSize.list  =  c(200),
              mu1.list  =  c(0.0, 0.5, 0.5,-1.0),
              mu2.list  =  c(0.5, 0.0, 1.0, 0.7),
              bound.list       =  c(1.5, 2.0, 2.5),
              nondecision.list =  c(0.1, 0.2, 0.3))
settings <- c(settings,
              s.topIdx  = length(settings$sampleSize.list),
              m1.topIdx  = length(settings$mu1.list),
              m2.topIdx  = length(settings$mu2.list),   # 'm' for magnitude
              b.topIdx  = length(settings$bound.list),
              n.topIdx  = length(settings$nondecision.list))
settings <- c(settings,
              possible.combinations = settings$s.topIdx * settings$m1.topIdx * settings$m2.topIdx * settings$b.topIdx * settings$n.topIdx,
              n.iter    = 2500,
              n.burnin  = 500,
              n.chains  = 4,
              n.thin    = 1,
              output.folder = "./output_cart_test/",
              simstudy.Name = "sim_")

################################################################
# Define simulation functions
################################################################
######  1. Function to GENERATE DATA 
####################################
generate <- function(seed, sampleSize, true.mu1, true.mu2, true.bound, true.nondecision) {
    set.seed(seed)
    X <-  cddm.simData.cart(trials = sampleSize,
                            mu1 = true.mu1,
                            mu2 = true.mu2,	
                            boundary = true.bound,
                            ndt = true.nondecision)$data
    return(X)
}

######  2. Function to RECOVER PARAMETER VALUES
###############################################
recover <- function(data, settings, seed) {
      ########## Write JAGS model  #############################################
      modelFile <- paste("cddm-C_seed-", seed, ".bug",sep="")
      write('
              data{
                    tmin <- 0.95 * min(X[,2])
                    D <- dim(X)
                    N <- D[1]
                    }
              model{
                    # Likelihood
                      for (i in 1:N) {
                           X[i,1:2] ~ dcddm(xdrift, ydrift, bound, ter0)
                      }
                    # Priors
                      xdrift  ~ dnorm(0, .1)
                      ydrift ~ dnorm(0, .1)
                      bound  ~ dgamma(3, 2)
                      ter0   ~ dexp(1)T(, tmin)
                    }', 
            modelFile)
      data <- list(X=data)
      parameters <- c("xdrift", "ydrift","bound", "ter0")
    
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
      mu1.samples  <- samples$BUGSoutput$sims.array[,,"xdrift"]
      mu2.samples       <- samples$BUGSoutput$sims.array[,,"ydrift"]
      bound.samples        <- samples$BUGSoutput$sims.array[,,"bound"]
      ndt.samples          <- samples$BUGSoutput$sims.array[,,"ter0"]
      
      ### Compute mean posterior per parameter
      est.mu1         <- mean( mu1.samples    )
      est.mu2         <- mean( mu2.samples    )
      est.bound       <- mean( bound.samples  )
      est.nondecision <- mean( ndt.samples    )
      mean.est        <- c(est.mu1, est.mu2, est.bound, est.nondecision)
      
      # Compute MAP per parameter
      map.mu1           <- myJAGS.MAP( mu1.samples    )
      map.mu2           <- myJAGS.MAP( mu2.samples    )
      map.bound         <- myJAGS.MAP( bound.samples  )
      map.nondecision   <- myJAGS.MAP( ndt.samples    )
      map.est           <- c(map.mu1, map.mu2, map.bound, map.nondecision)
      
      # Compute standard deviation
      sd.mu1          <- sd( mu1.samples   )
      sd.mu2          <- sd( mu2.samples   )
      sd.bound        <- sd( bound.samples )
      sd.nondecision  <- sd( ndt.samples   )
      sd.est <- c(sd.mu1, sd.mu2, sd.bound, sd.nondecision)
    
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
run_simulation <- function(seed, settings){
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
      mu1.list          <-  settings$mu1.list
      mu2.list          <-  settings$mu2.list
      bound.list        <-  settings$bound.list
      nondecision.list  <-  settings$nondecision.list
      # Identify maximum valid index for every list
      s.topIdx  <-  settings$s.topIdx
      m1.topIdx  <-  settings$m1.topIdx
      m2.topIdx  <-  settings$m2.topIdx
      b.topIdx  <-  settings$b.topIdx
      n.topIdx  <-  settings$n.topIdx
      output.folder         <- settings$output.folder
      possible.combinations <- settings$possible.combinations
      
      ########## Print a file to indicate the start of a new seed ###############
      init.File <- paste(output.folder,"seed-",seed,"_start.txt",sep="")
      write('Seed has been initiated', init.File)
      
      ########## Run simulation for every Parameter combination #################
      out <- list()
      for(s in 1:s.topIdx){
          for(b in 1:b.topIdx){
              write(paste("b =",b), file=init.File, append=TRUE)
              for(m1 in 1:m1.topIdx){
                  for(m2 in 1:m2.topIdx){
                      for(n in 1:n.topIdx){
                            sampleSize       <-   sampleSize.list  [s]
                            true.mu1         <-   mu1.list         [m1] 
                            true.mu2         <-   mu2.list         [m2]
                            true.bound       <-   bound.list       [b]
                            true.nondecision <-   nondecision.list [n]
                
                            this.truth        <-  c(true.mu1,
                                                    true.mu2,
                                                    true.bound,
                                                    true.nondecision)
                            
                            X   <-  generate(seed = seed,
                                             sampleSize = sampleSize,
                                             true.mu1 = true.mu1,
                                             true.mu2 = true.mu2,
                                             true.bound = true.bound,
                                             true.nondecision = true.nondecision)
                            Y   <-  recover(data = X, settings = settings, seed = seed)
                            out <- rbind(out,list(seed = seed,
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
      
      ########## Print a file to indicate the end of a the seed ###############
      end.File <- paste(output.folder,"seed-",seed,"_end.txt",sep="")
      write('Seed has ended running', init.File)
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
      studyName             <- settings$simstudy.Name
      n.iter    <-  settings$n.iter
      n.burnin  <-  settings$n.burnin
      
      ########## Create empty arrays to save output #############################
      ### Size variables
      # Parameter labels
      par.labels <- c("mu1","mu2", "bound","ndt")
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
   
      out.size <- possible.combinations * iterations
      for(set in 1:possible.combinations){
          J <- seq(set,out.size,possible.combinations) 
          for(i in 1:iterations){
              j <- J[i]
              S <- output[j,]
              seeds[i,set] <- S$seed
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
      colnames(rhats) <- names(S$rhats)
      save(rhats, file = paste(output.folder,studyName,"Rhats.RData",sep=""))
      save(trueValues, file = paste(output.folder,studyName,"trueValues.RData",sep=""))
      save(retrievedValues, file = paste(output.folder,studyName,"meanPosteriors.RData",sep=""))
      save(retrievedValues_sd, file = paste(output.folder,studyName,"std.RData",sep=""))
      save(mapValues, file = paste(output.folder,studyName,"MAPs.RData",sep=""))
      save(settings, file = paste(output.folder,studyName,"settings.RData",sep=""))
}

################################################################
# Define simulation functions
################################################################
cores       <-  detectCores()
my.cluster  <-  makeCluster(cores[1]-4)

registerDoParallel(cl = my.cluster)
output <- foreach(i = 1:settings$iterations, 
                    .errorhandling = "pass",
                    .combine = 'rbind'
                  ) %dopar% {
                    Z <- run_simulation(seed = i, settings)
                  }
stopCluster(cl = my.cluster)

store_output(output,settings)

######################################################################################
#####  LOAD OUTPUT
######################################################################################
output.File <- settings$output.folder

load(paste(output.File,"/sim_Rhats.RData",sep=""))
load(paste(output.File,"/sim_trueValues.RData",sep=""))
load(paste(output.File,"/sim_std.RData",sep=""))
load(paste(output.File,"/sim_meanPosteriors.RData",sep=""))
load(paste(output.File,"/sim_MAPs.RData",sep=""))
load(paste(output.File,"/sim__timers.RData",sep=""))
load(paste(output.File,"/sim__seeds.RData",sep=""))

######################################################################################
#####  CONVERGENCE CHECKS
######################################################################################
rule <- 1.05
bad.Rhat <- which(rhats>rule,arr.ind = TRUE)

test.rhat <- nrow(bad.Rhat)==0

if(!test.rhat){
      hist(rhats[rhats>rule],breaks=1000)
      abline(v=rule, col="red", lty=2)
      
      hist(rhats,breaks=1000)
      abline(v=rule, col="red", lty=2)
      legend("top",paste("Rhat = ",rule," | ",
                         (round(nrow(bad.Rhat)/prod(dim(rhats)),5))*100,
                         "% of chains"), lty=2, col="red", cex=0.4)
      
      bad.Rhat.ID <- colnames(rhats[,bad.Rhat[,2],])
      table(bad.Rhat.ID)
}else{
      paste("No Rhat greater than ", rule, sep="")      
}

######################################################################################
######            P L O T S                                  #########################
######################################################################################
boxplot.perPar <- function(trueValues,parameter.name, color="blue"){
  par.levels <- table(trueValues[,parameter.name])
  par.values <- as.numeric(names(par.levels))
  group.by.level <- matrix(NA,nrow=max(par.levels)*settings$iterations,
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

boxplot.perPar("mu1", "blue")
boxplot.perPar("mu2", "indianred")
boxplot.perPar("ndt", "green")
boxplot.perPar("bound", "purple")
