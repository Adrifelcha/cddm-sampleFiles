library(R2jags)  #Load up packages
load.module("cddm")

datos <- read.csv("./toyData.csv")
datos <- t(datos)

n.chains = 1
n.iter = 100 
n.burnin = 0
n.thin = 1

X <- datos[,1:10]
N <- ncol(X)

modelFile <- "cddm.bug"
    write('
            model{
                  # Likelihood
                    for (i in 1:N) {
                         X[1:2,i] ~ dcddm(drift, bound, ter0, theta0)
                    }
                      
                  # Priors
                    drift ~ dnorm(0,1)T(0,)
                    bound ~ dunif(0,4)
                    ter0 ~ dexp(1)
                    theta0 ~ dunif(0,6.283185)
                  }',
                      modelFile)

data <- list("X","N")
parameters <- c("drift", "bound", "ter0", "theta0")

fileName <- "toyExample_samples.RData"
samples <- jags(data=data, parameters.to.save=parameters, model=modelFile, n.chains=n.chains,
                n.iter=n.iter, n.burnin=n.burnin, n.thin=n.thin, DIC=T)
save(samples,file=fileName)
