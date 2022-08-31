# Load Rscripts with functions required to generate data
source("../Functions/fun_genParameters.R")
source("../Functions/fun_simData.R")

# Parameter values used in the simulation
trials = 200
true.bound = 2.45
mu1 <- 1.5
mu2 <- 1.25
true.ndt <- 0.1
true.drift.Angle <- cddm.getVectorAngle(mu1,mu2)
true.drift.Length <- cddm.getVectorLength(mu1,mu2)

# Evaluate if a data file already exists
searchFile <- "./toyData.csv"
      #If file does not exist, we simulate data and write a file
      if(!file.exists(searchFile)){
         simulateData <- cdd.simData(trials,true.drift.Angle,true.drift.Length,true.bound,true.ndt)
         write.csv(simulateData,"./toyData.csv",row.names = FALSE)
      }

# Read data
datos <- read.csv("./toyData.csv")    
# Print and plot the data
cddm.plotData(datos)
