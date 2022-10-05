################################################################################
##  The toyExample.R file runs the simple-example 
##  The simple-example is meant to be self-contained:
##      1. Calls the Rscript toyData.R to load/generate the dataset
##      2. Writes a simple JAGS model to implement the CDDM
##      3. Runs the model on the data
##      4. Plots the posterior chains
##  None of the function files contained in the Functions/ folder is used in this
##      first, simple-example.
################################################################################

# Step 0. Load packages
################################################################################
library(R2jags)  
load.module("cddm")

# Step 1. Load and prepare data
################################################################################
source("./toyData.R")    # Run Rscript that loads/generates data file as needed
datos <- t(datos)        # Data file needs to be transposed (2-by-trials)

X <- datos
N <- ncol(X)

# Step 2. Write simple JAGS model
################################################################################
modelFile <- "cddm_JAGSmodel.bug"
write('
        data {
              tmin <- 0.95 * min(X[2,])
        }
  
        model{
              # Likelihood
                for (i in 1:N) {
                     X[1:2,i] ~ dcddm(drift, bound, ter0, theta0)
                }
                  
              # Priors
                drift ~ dnorm(0,0.1)T(0,)
                bound ~ dgamma(3,2)
                ter0 ~ dexp(1)T(,tmin)
                theta0 ~ dnorm(0,0.1)T(-3.14,3.14)
              }',
      modelFile)

# Steo 3. Define settings to be passed to JAGS
################################################################################
n.chains = 3
n.iter = 10500 
n.burnin = 500
n.thin = 1

data <- list("X","N")
parameters <- c("drift", "bound", "ter0", "theta0")

# Step 4. Run JAGS and save samples
################################################################################
fileName <- "toyExample_samples.RData"
samples <- jags(data=data, parameters.to.save=parameters, model=modelFile, 
                n.chains=n.chains, n.iter=n.iter, n.burnin=n.burnin, 
                n.thin=n.thin, DIC=T)
save(samples,file=fileName)
load(fileName)

# Step 5. Show chains obtained per parameter
################################################################################
posterior.samples <- samples$BUGSoutput$sims.array

labels <- names(posterior.samples[1,1,])
for(i in 1:dim(posterior.samples)[3]){
    plot(posterior.samples[,1,i], type="l", main=labels[i], xlab="Iteration",
         ylab="Value sampled")  
    if(n.chains>1){
         for(a in 2:n.chains){
         lines(posterior.samples[,a,i],col=a)
         }
      }
}

# Step 6. Extract chains specific to every parameter in the CDDM
################################################################################
drift <- samples$BUGSoutput$sims.array[,,"drift"]
bound <- samples$BUGSoutput$sims.array[,,"bound"]
ter0 <- samples$BUGSoutput$sims.array[,,"ter0"]
theta0 <- samples$BUGSoutput$sims.array[,,"theta0"]

# Step 7. Plot posterior densities against true values
################################################################################
plotFunction <- function(samples,true.value){
    support <- round(seq(min(samples),max(samples),length.out = 10),2)
    plot(density(c(samples)), main="",axes=F,ann=F)
    abline(v=true.value, col="indianred4",lwd=2)
    legend("topright","true value",lwd=2,col="indianred4",cex=0.8,bty="n")
    mtext("Posterior values",1,line=2,f=2)
    mtext("Density",2,line=0,f=2)
    axis(1,support,support)
    mtext(paste("Posterior density - ",substitute(samples),sep=""),3, line=0.5,f=2)
}

plotFunction(drift,true.drift.Length)
plotFunction(bound,true.bound)
plotFunction(ter0,true.ndt)
plotFunction(theta0,true.drift.Angle)
