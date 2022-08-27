########################################################
########################################################
#####   A set of functions to generate some data under
#####   the Circular Drift Diffusion Model
########################################################
###########################   by Adriana F. Chavez   ###
source("./fun_genParameters.R")

# Variable dictionary: ##################################################
# mu1 and mu2 - Individual drift rates for the motion on the x and y axes
# drift.Angle - Direction of the drift vector
# drift.Length - Magnitude of the drift vector
# thresh - Boundary (radius)
# ndt - Non decision time
# drift.Coeff - Within-trial variability on the sampling process
# dt - Step size ("delta-t")
# state - rectangular coordinates recorded during the random walk
#########################################################################

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulate the full random walk across many trials (for each trial, 
# keeps the full chain of coordinates visited and response times)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
myCDDMsimData.randomWalk <- function(trials, mu1, mu2, thresh, ndt=0.1, drift.Coeff=1, dt=0.15){
      sqDT <- sqrt(dt)
      s.init <- c(0,0) 
      iter <- 20000  # Maximum number of iterations on the random walk process
      state <- array(NA, dim = c(iter, 2, trials))   #States are saved in a 3dimensional array
      finalT <- NA #An additional empty vector to store RT (a.k.a. total number of iterations)
          
      for(a in 1:trials){  
              state[1,,a] <- s.init   # States are stored in 3D array
            
              for(t in 2:iter){
                  d1 <- rnorm(1,mu1*dt,drift.Coeff*sqDT)
                  d2 <- rnorm(1,mu2*dt,drift.Coeff*sqDT)
                  state[t,,a] <- state[t-1,,a]+c(d1,d2)
                  pass <- sqrt(sum(state[t,,a]^2))
              
                  if(pass >= thresh){
                     finalT[a] <- t+(ndt/dt)   #Total no. of iterations required on each trial
                     break
                  }
              }
            
              if(pass > thresh){
                    digits.in.threshold <- nchar(sub('.*\\.', '', thresh))
                    d <- digits.in.threshold
                    A <- state[t-1,,a];  B <- state[t,,a]
                    last.step.x <- seq(A[1],B[1],length.out = 10^(d+1))
                    last.step.y <- seq(A[2],B[2],length.out = 10^(d+1))
                     
                    for(i in 1:length(last.step.x)){
                        pass2 <- sqrt(last.step.x[i]^2+last.step.y[i]^2)
                        if(round(pass2,digits.in.threshold) == thresh){
                             circunf <- c(last.step.x[i],last.step.y[i])
                             break
                        }
                    }
                    state[t,,a] <- circunf
              }
      }
  finalT <- finalT*ndt
  output <- list(state,finalT)
  names(output) <- c("state","RT")
  return(output)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Take full random walk coordinates and extract final response
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
myCDDMsimData.getCoordinates <- function(randomWalk){
   K <- nrow(randomWalk)
   I <- dim(randomWalk)[3]
   
   coord <- matrix(NA, ncol=2,nrow=I)
     for(i in 1:I){
           for(k in 1:K){
               if(!is.na(randomWalk[k,1,i])){
                   a <- k
               }else{
                   break
               }
           }
         coord[i,] <- randomWalk[a,,i]
     }
   return(coord)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Transform rectangular coordinates to degrees
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
myCDDMsimData.getDegrees <-  function(coord){
    x <- coord[,1]
    y <- coord[,2]
    theta <- atan2(y,x) * 180 / pi  #Angle with respect of y=0
        if(sum(theta<0)>0){
           theta.0 <- which(theta<0)
           theta[theta.0] <- theta[theta.0]+360   #Correction for whole circle (360)
        }   
  return(theta)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Transform degrees into radians
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
myCDDMsimData.getRadians <- function(theta.d){  
    theta <-  theta.d * pi /180  #Transform to radians
    return(theta)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulate data from the 4 parameters used to implement the cddm jags 
# module (with default values for the drift.Coefficient and dt)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
myCDDMsimData.simData <- function(trials, drift.Angle, drift.Length, thresh, ndt=0.1, drift.Coeff=1, dt=0.15){
  
      Mu <- myCDDM.getRectCoordinates(drift.Angle,drift.Length)
      mu1 <- Mu$mu1
      mu2 <-Mu$mu2
  
      randomWalk <-  myCDDMsimData.randomWalk(trials,mu1=mu1,mu2=mu2,drift.Coeff,thresh)
      RT <- randomWalk$RT
      randomWalk <- randomWalk$state
      coord <- myCDDMsimData.getCoordinates(randomWalk)
      degrees <- myCDDMsimData.getDegrees(coord)
      radians <- myCDDMsimData.getRadians(degrees)
      
      data <- as.data.frame(cbind(radians,RT))
      colnames(data) <- c("Choice","RT")
      
  return(data)
}
