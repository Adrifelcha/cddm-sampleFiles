##########################################################
# Load simulation study settings and output
##########################################################
# Hard-code the location of the simulation settings
load(paste("../output_ene23-2/settings.RData",sep=""))
output.File <- paste(".",settings$output.folder,sep="")
# Load the remaining output
load(paste(output.File,"/sim_Rhats.RData",sep=""))
load(paste(output.File,"/sim_trueValues.RData",sep=""))
load(paste(output.File,"/sim_std.RData",sep=""))
load(paste(output.File,"/sim_meanPosteriors.RData",sep=""))
load(paste(output.File,"/sim_MAPs.RData",sep=""))
load(paste(output.File,"/sim__timers.RData",sep=""))
load(paste(output.File,"/sim__seeds.RData",sep=""))


##########################################################
#  Functions
##########################################################
group.by.level <- function(trueValues.matrix, output.set,
                           parameter.name){
        par.levels <- table(trueValues.matrix[,parameter.name])
        par.values <- as.numeric(names(par.levels))
        groups <- matrix(NA, 
                         nrow = max(par.levels)*settings$iterations,
                         ncol = length(par.values))
        for(i in 1:length(par.values)){
            a <- par.values[i]
            same.level <- trueValues.matrix[,parameter.name]==a
            groups[,i] <- output.set[,parameter.name,same.level]
        }
        output <- list("groups" = groups,
                       "par.values" = par.values)
        return(output)
}


locate.par.sets <- function(trueValues, 
                            value, par.selected){
        output <- which(trueValues[,par.selected]==value)
        return(output)
}

rmv.true.value <- function(trueValues, output.set,
                           bad.value, rmv.parameter){
        bad.par.sets <- locate.par.sets(trueValues,
                                        value = bad.value,
                                        par.selected = rmv.parameter)
        d <- which(dim(output.set)==nrow(trueValues))
        
        if(unique(dim(output.set)==dim(trueValues))=="TRUE"){
              if(d==1){
                        output <- output.set[-bad.par.sets,]
              }else{
                        output <- output.set[,-bad.par.sets]      
              }
        }else{
              if(d==1){
                        output <- output.set[-bad.par.sets,,]
              }else{
                    if(d==2){
                              output <- output.set[,-bad.par.sets,]      
                    }else{
                              output <- output.set[,,-bad.par.sets]      
                    }
              }
        }
        return(output)
}

##########################################################
############                  P L O T S
##########################################################

# Plot 1:
# 
##################################################
par(mfrow=c(2,2))
par(mar = c(2.5, 3.5, 2, 2.5))

ndt <- group.by.level(trueValues, retrievedValues,
                      parameter.name="ndt")
ndt.groups <- ndt$groups
boxplot(ndt.groups, pch=16, cex=0.5, ann=F, axes=F,
        col=c("#88E8E8","#C3E7E7","#26F3F3"),
        ylim=c(0,0.5), notch=FALSE)
for(i in 1:3){
    lines(c(i-0.55,i+0.55),rep(ndt$par.values[i],2), col="indianred4", lty=1,lwd=2)
    }
mtext(paste("Non-decision time"),3, line = 0.2, f=1, cex=0.8)
axis(1,c(1:3),paste("Level ",1:3,sep=""))
axis(2,seq(0,0.5,0.05),
     format(round(seq(0,0.5,0.05), 2), nsmall = 2), 
     las=2, tck=-0.02,mgp=c(3, .5, 0))
legend("topleft","True values", lty=1, lwd=2, bty = "n", col="indianred4", cex=0.8)


bound <- group.by.level(trueValues, retrievedValues,
                        parameter.name="bound")
bound.groups <- bound$groups
boxplot(bound.groups, pch=16, cex=0.5, ann=F, axes=F,
        col=c("#FDB040","#FFB954","#FFD395"),
        ylim=c(1.4,3))
for(i in 1:3){
  lines(c(i-0.55,i+0.55),rep(bound$par.values[i],2), col="indianred4", lty=1,lwd=2)
}
mtext(paste("Boundary"),3, line = 0.2, f=1, cex=0.8)
axis(1,c(1:3),paste("Level ",1:3,sep=""))
axis(2,seq(1.4,3,0.2),
     format(round(seq(1.4,3,0.2), 2), nsmall = 1), 
     las=2, tck=-0.02,mgp=c(3, .5, 0))
legend("topleft","True values", lty=1, lwd=3, bty = "n", col="indianred4", cex=0.8)


dlength <- group.by.level(trueValues, retrievedValues,
                          parameter.name="driftLength")
dl.groups <- dlength$groups
boxplot(dl.groups, pch=16, cex=0.5, ann=F, axes=F,
        col=c("#C494E1","#C57AF5","#E8D0F8"),
        ylim=c(0,2.4))
for(i in 1:3){
  lines(c(i-0.55,i+0.55),rep(dlength$par.values[i],2), col="indianred4", lty=1,lwd=2)
}
mtext(paste("Drift length"),3, line = 0, f=1, cex=0.8)
axis(1,c(1:3),paste("Level ",1:3,sep=""))
axis(2,seq(0,2.4,0.4),
     format(round(seq(0,2.4,0.4), 2), nsmall = 1), 
     las=2, tck=-0.02,mgp=c(3, .5, 0))
legend("topleft","True values", lty=1, lwd=3, bty = "n", col="indianred", cex=0.8)


dangle <- group.by.level(trueValues, retrievedValues,
                         parameter.name="driftAngle")
da.groups <- dangle$groups
boxplot(da.groups, pch=16, cex=0.5, ann=F, axes=F,
        col=c("#FF988A","#FF5D46","#FFB0A5"),
        ylim=c(-3,3))
        for(i in 1:3){
            lines(c(i-0.55,i+0.55), rep(c(0,2,-2)[i],2), 
                  col="indianred4", lty=1,lwd=2)
        }
mtext(paste("Drift angle"),3, line = 0, f=1, cex=0.8)
axis(1,c(1:3),paste("Level ",1:3,sep=""))
axis(2,seq(-3,3,1), format(round(seq(-3,3,1), 2), nsmall = 1), 
     las=2, tck=-0.02,mgp=c(3, .5, 0))



# Plot 2:
# 
##################################################
par(mfrow=c(2,2))
par(mar = c(2.5, 3.5, 2, 2.5))

retrievedValues_restricted <- rmv.true.value(trueValues,
                                             output.set = retrievedValues,
                                             bad.value = 0.01,
                                             rmv.parameter = "driftLength")
trueValues_restricted <- rmv.true.value(trueValues,
                                        output.set = trueValues,
                                        bad.value = 0.01,
                                        rmv.parameter = "driftLength")

par(mfrow=c(2,2))
par(mar = c(2.5, 3.5, 2, 2.5))

ndt <- group.by.level(trueValues_restricted, retrievedValues_restricted,
                      parameter.name="ndt")
ndt.groups <- ndt$groups
boxplot(ndt.groups, pch=16, cex=0.5, ann=F, axes=F,
        col=c("#88E8E8","#C3E7E7","#26F3F3"),
        ylim=c(0,0.5), notch=FALSE)
for(i in 1:3){
  lines(c(i-0.55,i+0.55),rep(ndt$par.values[i],2), col="indianred4", lty=1,lwd=2)
}
mtext(paste("Non-decision time"),3, line = 0.2, f=1, cex=0.8)
axis(1,c(1:3),paste("Level ",1:3,sep=""))
axis(2,seq(0,0.5,0.05),
     format(round(seq(0,0.5,0.05), 2), nsmall = 2), 
     las=2, tck=-0.02,mgp=c(3, .5, 0))
legend("topleft","True values", lty=1, lwd=2, bty = "n", col="indianred4", cex=0.8)


bound <- group.by.level(trueValues_restricted, retrievedValues_restricted,
                        parameter.name="bound")
bound.groups <- bound$groups
boxplot(bound.groups, pch=16, cex=0.5, ann=F, axes=F,
        col=c("#FDB040","#FFB954","#FFD395"),
        ylim=c(1.4,3))
for(i in 1:3){
  lines(c(i-0.55,i+0.55),rep(bound$par.values[i],2), col="indianred4", lty=1,lwd=2)
}
mtext(paste("Boundary"),3, line = 0.2, f=1, cex=0.8)
axis(1,c(1:3),paste("Level ",1:3,sep=""))
axis(2,seq(1.4,3,0.2),
     format(round(seq(1.4,3,0.2), 2), nsmall = 1), 
     las=2, tck=-0.02,mgp=c(3, .5, 0))
legend("topleft","True values", lty=1, lwd=3, bty = "n", col="indianred4", cex=0.8)


dlength <- group.by.level(trueValues_restricted, retrievedValues_restricted,
                          parameter.name="driftLength")
dl.groups <- dlength$groups
boxplot(dl.groups, pch=16, cex=0.5, ann=F, axes=F,
        col=c("#C494E1","#C57AF5","#E8D0F8"),
        ylim=c(0.6,2.4))
for(i in 1:3){
  lines(c(i-0.55,i+0.55),rep(dlength$par.values[i],2), col="indianred4", lty=1,lwd=2)
}
mtext(paste("Drift length"),3, line = 0, f=1, cex=0.8)
axis(1,c(1:3),paste("Level ",1:3,sep=""))
axis(2,seq(0.6,2.4,0.3),
     format(round(seq(0.6,2.4,0.3), 2), nsmall = 1), 
     las=2, tck=-0.02,mgp=c(3, .5, 0))
legend("topleft","True values", lty=1, lwd=3, bty = "n", col="indianred", cex=0.8)


dangle <- group.by.level(trueValues_restricted, retrievedValues_restricted,
                         parameter.name="driftAngle")
da.groups <- dangle$groups
boxplot(da.groups, pch=16, cex=0.5, ann=F, axes=F,
        col=c("#FF988A","#FF5D46","#FFB0A5"),
        ylim=c(-3,3))
for(i in 1:3){
  lines(c(i-0.55,i+0.55), rep(c(0,2,-2)[i],2), 
        col="indianred4", lty=1,lwd=2)
}
mtext(paste("Drift angle"),3, line = 0, f=1, cex=0.8)
axis(1,c(1:3),paste("Level ",1:3,sep=""))
axis(2,seq(-3,3,1), format(round(seq(-3,3,1), 2), nsmall = 1), 
     las=2, tck=-0.02,mgp=c(3, .5, 0))



# Plot 3:
#
######################################################
par(mfrow=c(2,2))
par(mar = c(2.5, 3.5, 2, 2.5))

hist(retrievedValues[,4,which(trueValues[,4]==0.00&trueValues[,1]==0.01)],
     main="Angle = 0; Length = 0")
hist(retrievedValues[,4,which(trueValues[,4]==0.00&trueValues[,1]!=0.01)],
     main="Angle = 0; Length != 0")
hist(retrievedValues[,4,which(trueValues[,4]==4&trueValues[,1]==0.01)],
     main="Angle = -2; Length = 0")
hist(retrievedValues[,4,which(trueValues[,4]==4&trueValues[,1]!=0.01)],
     main="Angle = -2; Length != 0")

# Plot 4:
#
######################################################
layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
hist(retrievedValues[,4,which(trueValues[,4]&trueValues[,1]==0.01)])

retrievedValues_restricted <- rmv.true.value(trueValues,
                                             output.set = retrievedValues,
                                             bad.value = 0.01,
                                             rmv.parameter = "driftLength")
trueValues_restricted <- rmv.true.value(trueValues,
                                        output.set = trueValues,
                                        bad.value = 0.01,
                                        rmv.parameter = "driftLength")

dangle1 <- group.by.level(trueValues, retrievedValues,
                          parameter.name="driftAngle")
da.groups <- dangle1$groups
boxplot(da.groups, pch=16, cex=0.5, ann=F, axes=F,
        col=c("#FF988A","#FF5D46","#FFB0A5"),
        ylim=c(-3,3))
for(i in 1:3){
  lines(c(i-0.55,i+0.55),rep(c(0,2,-2)[i],2), col="indianred4", lty=1,lwd=2)
}
mtext(paste("Drift angle"),3, line = 0, f=1, cex=0.8)
axis(1,c(1:3),paste("Level ",1:3,sep=""))
axis(2,seq(-3,3,1),
     format(round(seq(-3,3,1), 2), nsmall = 1), 
     las=2, tck=-0.02,mgp=c(3, .5, 0))

dangle2 <- group.by.level(trueValues_restricted, 
                          retrievedValues_restricted,
                          parameter.name="driftAngle")
da.groups <- dangle2$groups
boxplot(da.groups, pch=16, cex=0.5, ann=F, axes=F,
        col=c("#FF988A","#FF5D46","#FFB0A5"),
        ylim=c(-3,3))
for(i in 1:3){
  lines(c(i-0.55,i+0.55),rep(c(0,2,-2)[i],2), col="indianred4", lty=1,lwd=2)
}
mtext(paste("Drift angle"),3, line = 0, f=1, cex=0.8)
axis(1,c(1:3),paste("Level ",1:3,sep=""))
axis(2,seq(-3,3,1),
     format(round(seq(-3,3,1), 2), nsmall = 1), 
     las=2, tck=-0.02,mgp=c(3, .5, 0))