##################################################################
##
##################################################################
iterations=200
possible.combinations = 3^5


makeBoxplot <- function(array.True, array.Retrieved, fix.by="driftAngle"){

  size.Levels <- unique(array.True[,"trials"])
  ndt.Levels <- unique(array.True[,"ndt"])
  fixPar.Levels <- unique(array.True[,fix.by])
  
  pars <- colnames(array.Retrieved)
  rotatingPars <- pars[-which(pars==fix.by)]
  
    for(ter in ndt.Levels){
        same.ndt <- array.True[,"ndt"]==ter
        for(size in size.Levels){
            same.size <- array.True[,"trials"]==size
            for(level in fixPar.Levels){
                same.fixedlevel <- array.True[,fix.by]==level
                
                keep <- same.ndt & same.size & same.fixedlevel
                selected.Samples <- array.Retrieved[,,keep]
                selected.Values <- array.True[keep,]
                
                boxplot(selected.Samples[,fix.by,], col="gray80")
                abline(h=level, col="indianred", lty=3, lwd=2)
                mtext(paste("n = ",size,"; ndt =", ter, sep=""),3)
                mtext(fix.by,2,line=2.25,f=2,cex=1.5)
                
                for(par in rotatingPars){
                    boxplot(selected.Samples[,par,], col="gray80")
                    abline(h=selected.Values[,par], col="indianred", lty=3)
                    mtext(paste(fix.by," = ",level,"; n = ",size,"; ndt =", ter, sep=""),3)
                    mtext(par,2,line=2.25,f=2,cex=1.5)
                }
              
           }
        }  
    }
}

makeBoxplot(trueValues, retrievedValues)