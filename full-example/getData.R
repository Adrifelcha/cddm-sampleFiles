#####################################################################
#####################################################################
#####  An Rscript for loading or generating random parameter values 
#####     and running a simulation under the CDDM framework.                                        
#####################################################################
#####################################################################

# Load Rscripts with functions required to generate data
source("../Functions/generateRandomParameterValues.R")
source("../Functions/simulateDataCDDM.R")
library(readr)

##########################################################################
# General settings

  if(!exists("ForceSample")){
      ForceSample <- FALSE  # By default, do not overwrite
    }
  
  if(!exists("trials")){
      trials = 200 #Standard number of trials
    }

dataFile.name <- "data.csv"
parameterObject.name <- "trueParValues.RData"
simulation.ID <- round(abs(runif(1,100,999)),0)

###########################################################################
# To prevent overwritting, we test if there's already a data file and 
# an R object containing the parameter True values

  # First, if a datafile does not exist yet, we initiate the simulation
  identify.dataFile <- file.exists(dataFile.name)
  if(!identify.dataFile){
     needSimulation <- TRUE
  }else{
        needSimulation <- FALSE
        identify.parFile <- file.exists(parameterObject.name)
           if(!identify.parFile){
              print("True parameter values unknown.")
           }
  }
  
  needSimulation <- needSimulation|ForceSample
  # If data doesn't exist, we generate new parameters and simulate data
  if(needSimulation){
      # Generate random parameter values
      par <- cddm.generateParameters()
      save(par,file=parameterObject.name)
      
      # Load parameters
      drift.Angle <- par["true.theta0"]
      drift.Length <- par["true.drift"]
      thresh <- par["true.bound"]
      ndt <- par["true.ter0"]
      
      # Generate data
      data <- cddm.simData(trials,drift.Angle,drift.Length,thresh,ndt)
      write.csv(data,dataFile.name,row.names = FALSE)
  }else{
      data <- read.csv(dataFile.name)
      load(parameterObject.name)
      
      if(trials!=nrow(data)){
        print("Warning: Loading an existing datafile with a different no. of trials. To overwrite this data object with the specified no. of trials, load a variable `ForceSample <- TRUE`")
      }
  }