# This function evaluates log[ theta^r p*(theta|y) ] where p*(.|.) is the 
# unnormalized posterior.
# This function is used by most of the implementation methods described
# below.
#
gfun <- function(theta,y,Expect,a,b,r) {
  const <- choose(N,y)
  if (r==0) gfun <- log(const) 
      (theta-a)^2/(2*b*b)
  if (r!=0) {
    # Integrals are defined wrt exp( r*log(theta) + log p*(theta|y) ) so
    # we must exclude theta<0 pts if we write in this form -- not a problem 
    # for these data (only used by Laplace and GH approaches)
    if (theta < 0) gfun <- -10000 
    else gfun <- log(const) + r*log(theta)-Expect*exp(theta)+y*theta-
        (theta-a)^2/(2*b*b)
  }
  gfun
}

gfun <- function(N=3330,y=50,delta=0.8, gamma=0.995, pi_){
  const <- choose(N, y)
  #(pi_*(delta+gamma-1) + (1-gamma))^y * (1-(pi_*(delta+gamma-1) + (1-gamma)))^(N-y)
  log(const) + y*log(pi_*(delta+gamma-1) + (1-gamma)) + (N-y)*log(1-(pi_*(delta+gamma-1) + (1-gamma)))
}


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Metropolis-Hastings algorithm
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
n <- 51000 # no if iterations
thetaMH <- NULL
thetahat <- exp(theta_hat)/(1+exp(theta_hat))
sethetahat <- sqrt(  p*(1-p) / (N*(delta + gamma -1)^2) )
thetaMH[1] <- thetahat
accept <- 0
for (i in 2:n){
  thetaMH[i] <- thetaMH[i-1]
  thetaprop <- max(rnorm(1,m=thetaMH[i-1],s=sqrt(3)*sethetahat),0)
  print(thetaprop)
  if (log(runif(1)) < gfun(pi_=thetaprop) - gfun(pi_=thetaMH[i-1])) {
    thetaMH[i] <- thetaprop; accept <- accept+1
  }
}
burnin <- 1000; indcalc <- seq(burnin+1,n)
# Estimate se's using batching
cat("MH:\n")
cat("Acceptance prob = ",accept/n,"\n")
B <- 1000 # length of each batch 
pmeanMH <- mean(thetaMH[indcalc])
pvarMH <- var(thetaMH[indcalc])
K <- length(indcalc)/B # K is the number of batches
batchmean <- batchvar <- NULL
for (i in 1:K){
  batchmean[i] <- mean(thetaMH[seq(burnin+1+(i-1)*K,burnin+i*K)])
  batchvar[i] <- var(thetaMH[seq(burnin+1+(i-1)*K,burnin+i*K)])
}
batchmeanvar <- var(batchmean)
batchvarvar <- var(batchvar)
batchmeanse <- sqrt(batchmeanvar/K)
batchvarse <- sqrt(batchvarvar/K)
par(mfrow=c(1,2))
#plot(thetaMH); abline(v=burnin)
plot(batchmean)
plot(batchvar)
cat("Posterior mean estimate: ",pmeanMH,pmeanMH-1.96*batchmeanse,pmeanMH+1.96*batchmeanse,"\n")
cat("Posterior var estimate: ",pvarMH,pvarMH-1.96*batchvarse,pvarMH+1.96*batchvarse,"\n")


#cat("Posterior variance estimate from batching: ",sum(batchest)/K,"\n")
#
# Lung cancer and radon example in Section 3.8
#
poisloglik <- function(x,y,Expected,beta0,beta1){
  mu <- rep(0,length(y))
  mu <- Expected*exp(beta0+beta1*x)
  poisloglik <- -sum(mu)+sum(y*log(mu))
  return(poisloglik)
}
# Read in data
lung <- read.table("MNlung.txt", header=TRUE, sep="\t")
radon <- read.table("MNradon.txt", header=TRUE)
Obs <- apply(cbind(lung[,3], lung[,5]), 1, sum)
Exp <- apply(cbind(lung[,4], lung[,6]), 1, sum)
rad.avg <- rep(0, length(lung$X))
for(i in 1:length(lung$X)) {
  rad.avg[i] <- mean(radon[radon$county==i,2])
}
x <- rad.avg
rad.avg[26]<-0
rad.avg[63]<-0
x[26] <- NA
x[63] <- NA
newy <- Obs[is.na(x)==F]
newx <- x[is.na(x)==F]
newE <- Exp[is.na(x)==F]


library(MASS)
#	Univariate random walk M-H 
nit <- 10000
beta0 <- rep(0,nit)
beta1 <- rep(0,nit)
mod <- glm(newy~offset(log(newE))+newx,family="poisson")
beta0[1] <- mod$coeff[1]
beta1[1] <- mod$coeff[2]
sd0 <- sqrt(vcov(mod)[1,1])
sd1 <- sqrt(vcov(mod)[2,2])
c <- .1
prevloglik <- poisloglik(x=newx,y=newy,Expected=newE,beta0[1],beta1[1])
count1 <- count2 <- count3 <- 0
for( i in 2:nit){
  beta0new <- rnorm(1,beta0[i-1],sd0*c)
  beta1new <- rnorm(1,beta1[i-1],sd1*c)
  newloglik <- poisloglik(newx,newy,newE,beta0new,beta1new)
  u <- runif(1)
  if (log(u) < newloglik-prevloglik){
    beta0[i] <- beta0new
    beta1[i] <- beta1new
    prevloglik <-  newloglik
    count1 <- count1 + 1
  }
  else {
    beta0[i] <- beta0[i-1]
    beta1[i] <- beta1[i-1]
  }
}
#
# Now run again with a different proposal (larger c)
#
iteration <- seq(1,1000,1)
gamma0 <- beta0
gamma1 <- beta1
beta0 <- rep(0,nit)
beta1 <- rep(0,nit)
mod <- glm(newy~offset(log(newE))+newx,family="poisson")
beta0[1] <- mod$coeff[1]
beta1[1] <- mod$coeff[2]
sd0 <- sqrt(vcov(mod)[1,1])
sd1 <- sqrt(vcov(mod)[2,2])
c <- 2
prevloglik <- poisloglik(newx,newy,newE,beta0[1],beta1[1])
for( i in 2:nit){
  beta0new <- rnorm(1,beta0[i-1],sd0*c)
  beta1new <- rnorm(1,beta1[i-1],sd1*c)
  newloglik <- poisloglik(newx,newy,newE,beta0new,beta1new)
  u <- runif(1)
  if (log(u) < newloglik-prevloglik){
    beta0[i] <- beta0new
    beta1[i] <- beta1new
    prevloglik <-  newloglik
    count2 <- count2 + 1
  }
  else {
    beta0[i] <- beta0[i-1]
    beta1[i] <- beta1[i-1]
  }
}
delta0 <- beta0
delta1 <- beta1



