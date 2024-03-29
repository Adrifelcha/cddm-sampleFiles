############################################################################################
###  The toyData.R script loads/generates the toyData.csv data set as needed.
###  The data is generated by simulating the random walk process as prescribed by the
###      function files fun_genParameters.R and fun_simData.R contained in the Functions/
###      folder.
###  The data is generated from a fixed set of parameters
############################################################################################
# Load customized functions
source("../Functions/generateRandomParameterValues.R")
source("../Functions/simulateDataCDDM.R")

# Step 0. Load the fixed parameter values used in the simulation
trials = 200
true.bound = 2.45
mu1 <- 1.5
mu2 <- 1.25
true.ndt <- 0.1
true.drift.Angle <- cddm.getVectorAngle(mu1,mu2)
true.drift.Length <- cddm.getVectorLength(mu1,mu2)

# Step 1. Determine whether toyData.csv already exists in this folder
searchFile <- "./toyData.csv"
file.not.found <- !file.exists(searchFile)

# Step 2. Load/generate toyData.csv data set
      # If file does not exist, we simulate data and save it into a brand new toyData.csv file
      if(file.not.found){
		# Run simulation
		      simulateData <- cddm.simData(trials = trials, 
						   drift.Angle = true.drift.Angle,
						   drift.Length = true.drift.Length,
						   thresh = true.bound,
						   ndt = true.ndt)
		# Store data into a toyData.csv file
		       write.csv(simulateData,"./toyData.csv",row.names = FALSE)
      }
      # Read data
      datos <- read.csv("./toyData.csv")    

# Step 3. Plot the data
cddm.plotData(datos)
