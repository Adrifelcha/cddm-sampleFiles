###################################################
###################################################
#####   Functions for simulating observations of a 
#####     CDDM process, given some parameter values 
###################################################
###########################   by Adriana F. Chavez

myCDDMsimData.randomWalk <- function(trials, mu.1, mu.2, std.driftrate, thresh){
    s.init <- c(0,0) 
    iter <- 1000
    state <- array(NA, dim = c(iter, 2, trials))   #States are saved in a 3dimensional array
    finalT <- NA #An additional empty vector to store RT (a.k.a. total number of iterations)
    
    for(a in 1:trials){  
        state[1,,a] <- s.init   # States are stored in 3D array
      
        for(t in 2:iter){
            d1 <- rnorm(1,mu1.drift,std.driftrate)
            d2 <- rnorm(1,mu2.drift,std.driftrate)
            state[t,,a] <- state[t-1,,a]+c(d1,d2)
            pass <- sqrt(sum(state[t,,a]^2))
        
            if(pass >= thresh){
               finalT[a] <- t   #Total no. of iterations required on each trial
               break
            }
        }
      
        if(pass > thresh){
          A <- state[t-1,,a];   B <- state[t,,a]
          last.step.x <- seq(A[1],B[1],length.out = 100)
          last.step.y <- seq(A[2],B[2],length.out = 100)
           
          for(i in 1:length(last.step.x)){
              pass2 <- sqrt(last.step.x[i]^2+last.step.y[i]^2)
            
              if(round(pass2,1) == thresh){
                 circunf <- c(last.step.x[i],last.step.y[i])
                 break
                }
          }
          state[t,,a] <- circunf
        }
    }
    
    output <- list(state,finalT)
    names(output) <- c("state","RT")
    return(output)
}

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

myCDDMsimData.getRadians <- function(theta.d){  
    theta <-  theta.d * pi /180  #Transform to radians
    return(theta)
}


myCDDMsimData.simData <- function(trials, mu.1, mu.2, std.driftrate, thresh){
      randomWalk <-  myCDDMsimData.randomWalk(trials,mu.1,mu.2,std.driftrate,thresh)
      RT <- randomWalk$RT
      randomWalk <- randomWalk$state
      coord <- myCDDMsimData.getCoordinates(randomWalk)
      degrees <- myCDDMsimData.getDegrees(coord)
      radians <- myCDDMsimData.getRadians(degrees)
      
      output <- list(randomWalk,coord,degrees,radians)
      names(output) <- c("randomWalk","coord","degrees","radians")
      
      return(output)
}