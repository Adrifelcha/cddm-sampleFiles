################################################################################
################################### Load data
################################################################################
# Establish no. of trials
trials <- 350
# Call Rscript to generate simulated data / load it  if already existing
source("./getData.R")
head(data)
# Plot data
cddm.plotData(data)
# Print parameter values used to generate this data
par

################################################################################
################################### Model preparation
################################################################################
modelFile <- "cddm.bug"
write('
            model{
                  # Likelihood
                    for (i in 1:N) {
                         X[1:2,i] ~ dcddm(drift, bound, ter0, theta0)
                    }
                      
                  # Priors
                    drift  ~ dnorm(0, 0.1)T(0,)
                    theta0 ~ dnorm(0, 0.1)T(-3.14, 3.14)
                    bound  ~ dgamma(3, 2)
                    ter0   ~ dexp(1)T(, tmin)
                  }',
      modelFile)

# Prepare Settings to be passed to JAGS
n.chains = 4
n.iter = 2500
n.burnin = 500
n.thin = 1
perParticipant = FALSE
perTask = FALSE

sampling.Settings <- list(n.chains,n.iter,n.burnin,n.thin,perParticipant,perTask)
names(sampling.Settings) <- c("n.chains","n.iter","n.burnin","n.thin","perParticipant","perTask")

################################################################################
################################### Run JAGS model 
################################################################################
source("../Functions/runCDDMjags.R")
source("../Functions/plotJAGSsamples.R")
source("../Functions/processJAGSsamples.R")

samplesFile <- "samples.RData"


if(!file.exists(samplesFile)){ myJAGSsampling.CDDM(sampling.Settings,modelFile,samplesFile,data) }
  
load(file=samplesFile)
  
if(file.exists(samplesFile)){ samples }

# Extract posterior samples and check them!
#########################################################
drift  <- myJAGSsampling.extractSamples("drift",samples)
bound  <- myJAGSsampling.extractSamples("bound",samples)
ter0   <- myJAGSsampling.extractSamples("ter0",samples)
theta0 <- myJAGSsampling.extractSamples("theta0",samples)

plot.ShowAllChains(samples)

myJAGSsampling.Rhat.max(samples)

plot.PosteriorDensity(drift,par["true.drift"])
plot.PosteriorDensity(bound,par["true.bound"])
plot.PosteriorDensity(ter0,par["true.ter0"])
plot.PosteriorDensity(theta0,par["true.theta0"])

################################################################################
################################### Check posterior samples against EZ-CDDM
################################################################################
source("../Functions/ezcdm.R")

EZ <- ezcdm.fit(data$Choice,data$RT)

# First, compare EZ estimates against true parameter values
###########################################################
EZ
par

# Now, get point descriptors of posterior samples
###########################################################
map.theta0 <- JAGSoutput.maxDensity(theta0)
map.drift <- JAGSoutput.maxDensity(drift)
map.bound <- JAGSoutput.maxDensity(bound)
map.ter0 <- JAGSoutput.maxDensity(ter0)
MAPS <- c(map.theta0,map.drift,map.bound,map.ter0)
names(MAPS) <- c("map.theta0","map.drift","map.bound","map.ter0")

mean.theta0 <- mean(theta0)
mean.drift <- mean(drift)
mean.bound <- mean(bound)
mean.ter0 <- mean(ter0)
means <- c(mean.theta0,mean.drift,mean.bound,mean.ter0)
names(means) <- c("mean.theta0","mean.drift","mean.bound","mean.ter0")

# compare EZ estimates and true parameter values against point descriptors
##########################################################################
EZ
par
MAPS
means
