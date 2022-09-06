# Establish no. of trials
trials <- 350
# Call Rscript to generate simulated data / load it  if already existing
source("./getData.R")
dim(data)
# Plot data
cddm.plotData(data)

# Print parameter values used to generate this data
par

modelFile <- "cddm.bug"
write('
            model{
                  # Likelihood
                    for (i in 1:N) {
                         X[1:2,i] ~ dcddm(drift, bound, ter0, theta0)
                    }
                      
                  # Priors
                    drift  ~ dnorm(0, 1)
                    theta0 ~ dnorm(0, 1)T(-3.14, 3.14)
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

# Run JAGS model
source("../Functions/runCDDMjags.R")

samplesFile <- "samples.RData"
myJAGSsampling.CDDM(sampling.Settings,modelFile,samplesFile,data)
