################################################################################
# Load packages
################################################################################
library(R2jags)  
load.module("cddm")

################################################################################
# Load data
################################################################################
source("./toyData.R")
datos <- t(datos)

X <- datos
N <- ncol(X)

################################################################################
# Write simple JAGS model
################################################################################
modelFile <- "cddm.bug"
write('
            model{
                  # Likelihood
                    for (i in 1:N) {
                         X[1:2,i] ~ dcddm(drift, bound, ter0, theta0)
                    }
                      
                  # Priors
                    drift ~ dnorm(0,1)T(0,)
                    bound ~ dunif(0,5)
                    ter0 ~ dexp(1)T(,0.4)
                    theta0 ~ dunif(0,6.283185)
                  }',
      modelFile)

################################################################################
# Settings to be passed to JAGS
################################################################################
n.chains = 1
n.iter = 1000 
n.burnin = 0
n.thin = 1

data <- list("X","N")
parameters <- c("drift", "bound", "ter0", "theta0")

################################################################################
# Run JAGS
################################################################################
fileName <- "toyExample_samples.RData"
samples <- jags(data=data, parameters.to.save=parameters, model=modelFile, 
                n.chains=n.chains, n.iter=n.iter, n.burnin=n.burnin, 
                n.thin=n.thin, DIC=T)
save(samples,file=fileName)
load(fileName)
################################################################################
# Plot chains
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

################################################################################
# Get parameter values
################################################################################
drift <- samples$BUGSoutput$sims.array[,,"drift"]
bound <- samples$BUGSoutput$sims.array[,,"bound"]
ter0 <- samples$BUGSoutput$sims.array[,,"ter0"]
theta0 <- samples$BUGSoutput$sims.array[,,"theta0"]

################################################################################
# Plot posterior densities
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