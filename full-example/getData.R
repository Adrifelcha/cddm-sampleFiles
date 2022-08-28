#####################################################################
#####################################################################
#####  An Rscript for loading or generating random parameter values 
#####     and running a simulation under the CDDM framework.                                        
#####################################################################
#####################################################################

# Load Rscripts with functions required to generate data
source("../Functions/fun_genParameters.R")
source("../Functions/fun_simData.R")
library(readr)

##########################################################################
# General settings
trials = 200
ForceSample <- FALSE

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
      drift.Angle <- par$driftAngle
      drift.Length <- par$driftLength
      thresh <- par$thresh
      ndt <- par$ndt
      
      # Generate data
      data <- cdd.simData(trials,drift.Angle,drift.Length,thresh,ndt)
      write.csv(data,dataFile.name,row.names = FALSE)
  }else{
      data <- read.csv(dataFile.name)
      load(parameterObject.name)
  }


# Print and plot the data
head(data)