########################################################
########################################################
#####   A set of functions to generate some data under
#####   the Circular Drift Diffusion Model
########################################################
###########################   by Adriana F. Chavez   ###
library(plotrix) #Library to draw circles

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
cddm.randomWalk <- function(trials, mu1, mu2, thresh, ndt=0.1, drift.Coeff=1, dt=0.0015){
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
      
  finalT <- finalT*dt
  output <- list(state,finalT)
  names(output) <- c("state","RT")
  return(output)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Take full random walk coordinates and extract final response
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cddm.getFinalState <- function(randomWalk.states){
   randomWalk <- randomWalk.states
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
cddm.coordToDegrees <-  function(coord){
    x <- coord[,1]
    y <- coord[,2]
    theta <- atan2(y,x) * 180 / pi  #Angle with respect of y=0
        if(sum(theta<0)>0){
           theta.0 <- which(theta<0)
           theta[theta.0] <- theta[theta.0]+360   #Correction for whole circle (360)
        }   
  return(theta)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulate data from the 4 parameters used to implement the cddm jags 
# module (with default values for the drift.Coefficient and dt)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cdd.simData <- function(trials, drift.Angle, drift.Length, thresh, ndt=0.1, drift.Coeff=1, dt=0.015){
  
      Mu <- cddm.polarToRect(drift.Angle,drift.Length)
      mu1 <- Mu$mu1
      mu2 <-Mu$mu2
  
      randomWalk <-  cddm.randomWalk(trials=trials,mu1=mu1,mu2=mu2,thresh=thresh,
                                     ndt=ndt,drift.Coeff=drift.Coeff,dt=dt)
      RT <- randomWalk$RT
      randomWalk <- randomWalk$state
      coord <- cddm.getFinalState(randomWalk)
      degrees <- cddm.coordToDegrees(coord)
      radians <- cddm.degToRad(degrees)
      radians <- round(radians,2)
      
      data <- as.data.frame(cbind(radians,RT))
      colnames(data) <- c("Choice","RT")
      
  return(data)
}

############################################################################
#####  Plotting functions
#####  Note: The margins of the plotting space may need to be adjusted to 
#####        fully appreciate the symmetry of the circle drawn on screen.
############################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot the random walk (and RT distribution) from cddm.randomWalk()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cddm.plotRW <- function(randomWalk){
  state  <- randomWalk$state
  finalT <- randomWalk$RT
  trials <- length(finalT)
  choices <- cddm.getFinalState(state)
  thresh <- cddm.getVectorLength(choices[1,1],choices[1,2])
  thresh <- round(thresh,2)
  
  par(mfrow = c(1,2))  # Open space for 2 plots
  pm <- thresh+0.5 #Plot margin
  plot(-10:10,-10:10,type="n", ann = FALSE, axes = FALSE,
       xlim=c(-pm,pm),ylim=c(-pm,pm))
  for(b in 1:trials){
      points(state[,,b], type = "l", col=rgb(1,0,0.5,0.1))
  }
  draw.circle(0,0,radius = thresh, border = "black")
  abline(h = 0, lty=2, col="gray50")
  abline(v = 0, lty=2, col="gray50")
  legend("topright",paste("No. trials =", trials), 
         pch=16, col="white",bty = "n", cex=0.8)
  for(b in 1:trials){
      points(choices[b,1],choices[b,2], type = "p", pch =16, cex=0.9,
             col=rgb(0.75,0.25,0.5,0.2))
  }
  
  maxRT <- max(finalT)+5
  x.axis <- round(c(0,seq(0,maxRT,length.out=10)),2)
  hist(finalT, col = "darkorchid4", breaks = 50, ann=FALSE, axes=FALSE)
  mtext("Response Times", 1, line=2, f=2)
  mtext("Frequency", 2, line = 2.5, cex=0.8)
  axis(2, seq(0,trials,5), seq(0,trials,5), las=2)
  axis(1, x.axis,x.axis)
  
  par(mfrow = c(1,1)) #As a precaution, go back to single plot spaces
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot  observed choices and RT
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cddm.plotData <- function(simData){
      choice <- simData$Choice
      RT <- simData$RT
      trials <- length(RT)
      
      direction <- cddm.radToDeg(choice) # Transform radian choices into degrees
      thresh <- 9 # Arbitrary radius, used to define magnitude
      magnitude <- rep(thresh,length(choice)) 
      coord.on.circumference <- cddm.polarToRect(choice,magnitude) #get Rectangular coordinates
      
      par(mfrow = c(1,2))  # Open space for 2 plots
      
      plot(-10:10,-10:10,type="n", ann = FALSE, axes = FALSE)
      for(b in 1:trials){
        points(coord.on.circumference[b,1],coord.on.circumference[b,2], 
               type = "p", pch =16, cex=0.9,
               col=rgb(0.75,0.25,0.5,0.2))
      }
      draw.circle(0,0,radius = thresh, border = "black")
      abline(h = 0, lty=2, col="gray50")
      abline(v = 0, lty=2, col="gray50")
      legend("topright",paste("No. trials =", trials), 
             pch=16, col="white",bty = "n", cex=0.8)
      
      maxRT <- max(RT)+5
      x.axis <- round(c(0,seq(0,maxRT,length.out=10)),2)
      hist(RT, col = "darkorchid4", breaks = 50, ann=FALSE, axes=FALSE)
      mtext("Response Times", 1, line=2, f=2)
      mtext("Frequency", 2, line = 2.5, cex=0.8)
      axis(2, seq(0,trials,5), seq(0,trials,5), las=2)
      axis(1, x.axis,x.axis)
      
      par(mfrow = c(1,1)) #As a precaution, go back to single plot spaces
}
