load(file="./output_dec22/subset1_m1b1n1a1s1.RData")

npar <- ncol(Z$current.matrix.means)
for(i in 1:npar){
      boxplot(Z$current.matrix.means[,i])
      abline(h=Z$current.truth[i])
      
      plot(density(Z$current.matrix.means[,i]))
      abline(v=Z$current.truth[i], col="red")
}

boxplot(Z$current.matrix.means)
abline(h=Z$current.truth)

X <- Z$current.thetas %% (2*pi)
#X <- cddm.radToDeg(X)
Y <- Z$current.lengths

mu1 <- X
mu2 <- Y

for(a in 1:dim(X)[3]){
    for(i in 1:ncol(X)){
        mu1[,i,a] <- cddm.polarToRect(X[,i,a],Y[,i,a])$mu1
        mu2[,i,a] <- cddm.polarToRect(X[,i,a],Y[,i,a])$mu2
    }
}

Mu1 <- apply(mu1,3,mean)
Mu2 <- apply(mu2,3,mean)

est <- cddm.getVectorAngle(Mu1,Mu2)

boxplot(est)
abline(h=0)

source("../Functions/generateRandomParameterValues.R")
z <- cddm.polarToRect(x,y)

mu <- apply(z,2,mean)
