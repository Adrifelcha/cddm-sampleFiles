##################################################################
##
##################################################################
iterations=200
possible.combinations = 3^5

#load("trueValues.RData")
#load("retrievedValues.RData")
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

  size.Levels <- unique(array.True[,"trials",])
  ndt.Levels <- unique(array.True[,"ndt",])
  fixPar.Levels <- unique(array.True[,fix.by,])
  
  pars <- colnames(array.Retrieved)
  rotatingPars <- pars[-which(pars==fix.by)]
  
    for(ter in ndt.Levels){
        same.ndt <- array.True[1,"ndt",]==ter
        for(size in size.Levels){
            same.size <- array.True[1,"trials",]==size
            for(level in fixPar.levels){
                same.fixedlevel <- array.True[1,fix.by,]==level
                
                keep <- same.ndt & same.size & same.fixedlevel
                selected.Samples <- array.Retrieved[,,keep]
                selected.Values <- array.True[,,keep]
                
                boxplot(selected.Samples[,fix.by,], col="gray80")
                abline(h=level, col="indianred", lty=3, lwd=2)
                mtext(paste("n = ",size,"; ndt =", ter, sep=""),3)
                mtext(fix.by,2,line=2.25,f=2,cex=1.5)
                
                for(par in rotatingPars){
                    boxplot(selected.Samples[,par,], col="gray80")
                    abline(h=selected.Values[par,], col="indianred", lty=3)
                    mtext(paste(fix.by," = ",level,"; n = ",size,"; ndt =", ter, sep=""),3)
                    mtext(par,2,line=2.25,f=2,cex=1.5)
                }
              
           }
        }  
    }
}