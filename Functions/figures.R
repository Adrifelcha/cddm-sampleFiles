source("./simulateDataCDDM.R")
source("./generateRandomParameterValues.R")

set.seed(123)

#################################
## OVERALL SET UP
#################################
trials <- 100
drift.Angle <- 45
drift.Length <- 2
boundary <- 2
ndt <- 0.1
Mu <- cddm.polarToRect(drift.Angle,drift.Length)
mu1 <- Mu$mu1
mu2 <-Mu$mu2

#################################
## CREATE DATA
#################################
randomWalk <- cddm.randomWalk(trials, mu1, mu2, boundary, ndt)
state  <- randomWalk$state
finalT <- randomWalk$RT
trials <- length(finalT)
choices <- cddm.getFinalState(state)
boundary <- cddm.getVectorLength(choices[1,1],choices[1,2])
boundary <- round(boundary,2)

##################################################
## FIGURES 1 & 2: Illustrating single-trial
##################################################
par(mfrow = c(1,1))  

circle <- cddm.polarToRect(all.Angles,boundary)
plot(-10:10,-10:10,type="n", ann = FALSE, axes = FALSE,
     xlim=c(-pm,pm),ylim=c(-pm,pm))
points(circle[,1],circle[,2], type="l")
abline(h = 0, lty=2, col="gray50")
abline(v = 0, lty=2, col="gray50")

plot(-10:10,-10:10,type="n", ann = FALSE, axes = FALSE,
     xlim=c(-pm,pm),ylim=c(-pm,pm))
points(state[,,1], col="pink", type="l")
points(circle[,1],circle[,2], type="l")
abline(h = 0, lty=2, col="gray50")
abline(v = 0, lty=2, col="gray50")
points(choices[1,1],choices[1,2], type = "p", pch =16, 
       cex=1.2, col="purple")

##################################################
## FIGURE 3: Side-By-side Circle - RT
##################################################
par(mfrow = c(1,2))  # Open space for 2 plots
pm <- boundary+0.5 #Plot margin
plot(-10:10,-10:10,type="n", ann = FALSE, axes = FALSE,
     xlim=c(-pm,pm),ylim=c(-pm,pm))
for(b in 1:trials){
  points(state[,,b], type = "l", col=rgb(1,0,0.5,0.1))
}
points(circle[,1],circle[,2], type="l")
abline(h = 0, lty=2, col="gray50")
abline(v = 0, lty=2, col="gray50")
for(b in 1:trials){
  points(choices[b,1],choices[b,2], type = "p", pch =16, cex=0.9,
         col=rgb(0.75,0.25,0.5,0.3))
}

minRT <- min(finalT)
maxRT <- max(finalT)+0.05
x.axis <- round(seq(minRT,maxRT,length.out=10),2)
hist(finalT, col = "darkorchid4", breaks = 40, 
     ann=FALSE, axes=FALSE, ylim=c(0,1.75), 
     freq = FALSE)
#mtext("Response Time distributions", 3, line=0.5, f=2)
axis(1, x.axis,x.axis)

##################################################
## FIGURE 4: Side-By-side Circle - Angle choice
##################################################
par(mfrow = c(1,2))  # Open space for 2 plots
pm <- boundary+0.5 #Plot margin
plot(-10:10,-10:10,type="n", ann = FALSE, axes = FALSE,
     xlim=c(-pm,pm),ylim=c(-pm,pm))
for(b in 1:trials){
  points(state[,,b], type = "l", col=rgb(1,0,0.5,0.1))
}
points(circle[,1],circle[,2], type="l")
abline(h = 0, lty=2, col="gray50")
abline(v = 0, lty=2, col="gray50")
for(b in 1:trials){
  points(choices[b,1],choices[b,2], type = "p", pch =16, cex=0.9,
         col=rgb(0.75,0.25,0.5,0.3))
}

degrees <- cddm.coordToDegrees(choices)
radians <- cddm.degToRad(degrees)
radians <- round(radians,4)
x.axis <- round(seq(0,2*pi,length.out=10),2)
hist(radians, col = "blue4", breaks = 40, 
     ann=FALSE, axes=FALSE, ylim=c(0,0.9), 
     freq = FALSE, xlim=c(0,2*pi))
#mtext("Angle choice", 3, line=0.5, f=2)
axis(1, x.axis,x.axis)





##################################################
## Make new data with Angle = 0
##################################################
drift.Angle <- 0
Mu <- cddm.polarToRect(drift.Angle,drift.Length)
mu1 <- Mu$mu1
mu2 <-Mu$mu2

randomWalk <- cddm.randomWalk(trials, mu1, mu2, boundary, ndt)
state  <- randomWalk$state
finalT <- randomWalk$RT
trials <- length(finalT)
choices <- cddm.getFinalState(state)
boundary <- cddm.getVectorLength(choices[1,1],choices[1,2])
boundary <- round(boundary,2)

circle <- cddm.polarToRect(all.Angles,boundary)
par(mfrow = c(1,2))  # Open space for 2 plots
pm <- boundary+0.5 #Plot margin
plot(-10:10,-10:10,type="n", ann = FALSE, axes = FALSE,
     xlim=c(-pm,pm),ylim=c(-pm,pm))
    for(b in 1:trials){
        points(state[,,b], type = "l", col=rgb(1,0,0.5,0.1))
    }
points(circle[,1],circle[,2], type="l")
abline(h = 0, lty=2, col="gray50")
abline(v = 0, lty=2, col="gray50")
    for(b in 1:trials){
        points(choices[b,1],choices[b,2], type = "p", pch =16, cex=0.9,
               col=rgb(0.75,0.25,0.5,0.3))
    }

degrees <- cddm.coordToDegrees(choices)
radians <- cddm.degToRad(degrees)
radians <- round(radians,4)
x.axis <- round(seq(0,2*pi,length.out=10),2)
hist(radians, col = "blue4", breaks = 40, 
     ann=FALSE, axes=FALSE, ylim=c(0,0.8), 
     freq = FALSE, xlim=c(0,2*pi))
#mtext("Angle choice", 3, line=0.5, f=2)
axis(1, x.axis,x.axis)



#######################################
###
#######################################
trials <- 2

randomWalk <- cddm.randomWalk(trials, mu1, mu2, boundary, ndt,dt = 0.015)
state  <- randomWalk$state
finalT <- randomWalk$RT
trials <- length(finalT)
choices <- cddm.getFinalState(state)
boundary <- cddm.getVectorLength(choices[1,1],choices[1,2])
boundary <- round(boundary,2)

par(mfrow = c(1,1))  
plot(-10:10,-10:10,type="n", ann = FALSE, axes = FALSE,
     xlim=c(-pm,pm),ylim=c(-pm,pm))
points(state[,,1], col="pink3", type="l")
points(circle[,1],circle[,2], type="l")
abline(h = 0, lty=2, col="gray50")
abline(v = 0, lty=2, col="gray50")
points(choices[1,1],choices[1,2], type = "p", pch =16, 
       cex=1.2, col="purple")