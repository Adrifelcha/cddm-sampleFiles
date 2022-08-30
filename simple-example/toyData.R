# Load Rscripts with functions required to generate data
source("../Functions/fun_genParameters.R")
source("../Functions/fun_simData.R")

# Parameter values used in the simulation
trials = 200
thresh=2.45
mu1 <- 1.5
mu2 <- 1.25
ndt <- 0.1
drift.Angle <- cddm.getVectorAngle(mu1,mu2)
drift.Length <- cddm.getVectorLength(mu1,mu2)

# Evaluate if a data file already exists
searchFile <- "./toyData.csv"
      #If file does not exist, we simulate data and write a file
      if(!file.exists(searchFile)){
         simulateData <- cdd.simData(trials,drift.Angle,drift.Length,thresh,ndt)
         write.csv(simulateData,"./toyData.csv",row.names = FALSE)
      }

# Read data
y <- read.csv("./toyData.csv")    
# Print and plot the data
head(y,200)
cddm.plotData(y)
