library(R2jags)  #Load up packages
load.module("cddm")

datos <- read.csv("./toyData.csv")
datos <- t(datos)

n.chains = 2
n.iter = 1000 
n.burnin = 0
n.thin = 1

X <- datos[,1:10]
N <- length(X)

modelFile <- "cddm.bug"
    write('
            model{
                  # Likelihood
                    for (i in 1:N) {
                         X[1:2,] ~ dcddm(drift, bound, ter0, theta0)
                    }
                      
                  # Priors
                    drift ~ dnorm(0,1)
                    bound ~ dunif(0,1000)
                    ter0 ~ dexp(1)
                    theta0 ~ dunif(0,6.283185)
                  }',
                      modelFile)

data <- list("X","N")
parameters <- c("drift", "bound", "ter0", "theta0")

samples <- jags(data, parameters, model.file=model, n.chains=n.chains,
                n.iter=n.iter, n.burnin=n.burnin, n.thin=n.thin, DIC=T)
save(samples,file=fileName)