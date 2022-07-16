###################################################
###################################################
#####   Functions for running the CDDM module in  
#####                                          JAGS 
###################################################
###########################   by Adriana F. Chavez

library(R2jags)  #Load up packages
library(brms)
library(rstan)
library(readr)

###################################################
##### SAMPLING EX GAUSSIAN PARAMETERS
###################################################
myJAGSsampling.CDDM <- function(sampling.Settings,model,fileName){
      n.chains = as.numeric(sampling.Settings[1])
      n.iter = as.numeric(sampling.Settings[2])
      n.burnin = as.numeric(sampling.Settings[3])
      n.thin = as.numeric(sampling.Settings[4])
      perParticipant = sampling.Settings[5]
      perTask = sampling.Settings[6]
      
      if(!perParticipant){
          # data to be feed to JAGS
          data <- list("RT","C","N") 
          # parameters to be monitored:	
          parameters <- c("alpha", "sigma", "mu1","mu2")
      }else{
          if(!perTask){
            # data to be feed to JAGS
            data <- list("RT","N", "C", "P" ,"sub") 
          }else{
            # data to be feed to JAGS
            data <- list("RT","N", "C", "P", "sub","J","task") 
          }
          # parameters to be monitored:	
          parameters <- c("alpha", "sigma", "mu1","mu2",
                          "mu.alpha","sigma.alpha",
                          "mu.sigma","sigma.sigma",
                          "mu.mu1","sigma.mu1",
                          "mu.mu2","sigma.mu2")
      }
      # Start chain
  samples <- jags(data, parameters.to.save=parameters,
                  model.file=model, n.chains=n.chains,
                  n.iter=n.iter, n.burnin=n.burnin, n.thin=n.thin, DIC=T)
  save(samples,file=fileName)
  return(samples)
}


###################################################
##### EXTRACTING SAMPLES
###################################################
myJAGSsampling.extractSamples <- function(parameter.name, samples){
  # We use the sims.array and NOT sims.list and sims.matrix
  # https://sourceforge.net/p/mcmc-jags/discussion/610037/thread/cc61b820/
  postParam.Array <- samples$BUGSoutput$sims.array
  
  # Identify to which parameter each page corresponds to
  samplesID <- names(postParam.Array[1,1,])
  
  # Locate pages that contain samples for the Parameter of interest
  locateParameter <- which(grepl(parameter.name,samplesID))
  
  # Now, for individual/task level parameters...
  if(length(locateParameter)>1){
    
    # We isolate the pages ID that contain the parameter of interest
    samplesRelated <- samplesID[locateParameter]
    # We evaluate whether these vary across individuals AND tasks
    doubleIndex <- which(grepl(",",samplesRelated))
    
    # IF these parameters are estimated per individual AND task:  
    if(length(doubleIndex)>0){
      # We isolate every sample set with a double index
      locateParameter <- samplesRelated[doubleIndex]
      # If these parameters only vary across individuals OR tasks:    
    }else{
      # We identify the maximum index
      locateParticipantID <- parse_number(samplesRelated)
      nP <- max(locateParticipantID,na.rm = TRUE)
      # We locate the page ID containing the last index
      lastP <- which(grepl(nP,samplesRelated))
      # And isolate all IDs from there
      locateParameter <- samplesRelated[(lastP-nP)+1:lastP]
    }
  }
  
  # We retrieve only the pages containing the parameter
  x <- postParam.Array[,,locateParameter]
  return(x)
}


###################################################
##### EVALUATING SAMPLES: DEVIANCE AND RHAT
###################################################

myJAGSsampling.extractDeviance <- function(samples){
    x <- samples$BUGSoutput$sims.array
    dev <- x[,,"deviance"]
    return(dev)
}


myJAGSsampling.Rhat <- function(samples){
    x <- samples$BUGSoutput$sims.array
    Rhats <- apply(x,3,Rhat)
    return(Rhats)
}


myJAGSsampling.Rhat.max <- function(samples,maxRhat = 1.05){
      Rhats <- myJAGSsampling.Rhat(samples)
      exceedingRhats <- which(Rhats>maxRhat)
      if(length(exceedingRhats)>0){
        return(Rhats[exceedingRhats])
      }
      max <- max(Rhats)
      maxChain <- names(Rhats[which(Rhats==max(Rhats))])
      return(paste("The maximum value of Rhat observed was ", round(max,4), " which corresponds to parameter ", maxChain))
}
