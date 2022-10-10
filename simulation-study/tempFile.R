##################################################################
##
##################################################################
iterations=200
possible.combinations = 3^5

load("test_IndividualEstimates.RData")

trueValues <- array(NA, dim=c(1,5,possible.combinations))
retrievedValues <- array(NA,dim=c(iterations,4,possible.combinations))

par.labels <- c("driftLength","bound","ndt","driftAngle")
colnames(trueValues) <- c("trials",par.labels)
colnames(retrievedValues) <- par.labels

trueValues[,1,] <- rep(c(100,500,1000), each=3^4)
trueValues[,2,] <- rep(c(0.5,0.6,0.7), each=3^3)
trueValues[,3,] <- rep(c(1,1.5,2), each=3^2)
trueValues[,4,] <- rep(c(0.1,0.2,0.3), each=3)   
trueValues[,5,] <- rep(c(0,0.5,1),3^4)
retrievedValues[,,] <- Y


array.True <- trueValues
array.Retrieved <- retrievedValues
makeBoxplot <- function(array.True, array.Retrieved, fix.by="driftAngle"){

  fixPar.rec <- array.Retrieved[,fix.by,]
  fixPar.levels <- unique(array.True[,fix.by,])
  
  for(level in fixPar.levels){
    this.level <- which(array.True[1,fix.by,]==level)
    boxplot(fixPar.rec[,this.level])
    
  }
  
    
  nPars <- ncol(array.Retrieved)
  sizeLevels <- unique(array.True[,"trials",])
  
  for(par in 1:nPars){
    for(size in sizeLevels){
      this.size <- which(array.True[1,1,]==size)
      
      this.par <- colnames(array.Retrieved)[par]
      if(this.par==fix.by){
      
      }else{
        boxplot(array.Retrieved[,par,this.size])  
      }
      
    }
  }
}