# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Student: Raul
# Course: STAT 570
# Assignment: midterm
# Date: 2022/11/08
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
library(dplyr)
library(ggplot2)
rm(list = ls())


# ~~~~~~~~~~#
# Problem 1 #
# ~~~~~~~~~~#

## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## (b)  
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

N <- 3330
delta <- 0.8
gamma <- 0.995
y <- 50

# p stuff
p <- y/N
I_p <- N*(p-p^2)
partial_loglike <- y*(1/p) - ((N-y)/(1-p))

# pi stuff
pi_hat <- (y-N*(1-gamma))/(N*(delta+gamma-1)) 
I_pi <- (N*(delta + gamma -1)^2) / p*(1-p)
var_pi_hat <- p*(1-p) / (N*(delta + gamma -1)^2)
CI_wald <- c(pi_hat+qnorm(c(0.05,0.95), mean=0, sd=1)*sqrt(var_pi_hat))

## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## (d)
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
theta_hat <- log(pi_hat/(1-pi_hat))
I_theta <- I_pi * pi_hat^2*(1-pi_hat^2)
var_theta_hat <- 1/I_theta
exp(c(theta_hat -1.96 * sqrt(1/I_theta), theta_hat +1.96 * sqrt(1/I_theta)))



# ~~~~~~~~~ #
# Problem 2 #
# ~~~~~~~~~ #

## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## (a)  Posterior
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Grid method to gain intuition for what a solution should look like:
pi_vals <- seq(from = 0, to = 1, length = 1e2)
p_vals <- pi_vals*(delta + gamma -1) + (1-gamma)

grid_data <- data.frame(pi_vals = pi_vals, p_vals = p_vals)

grid_data <- grid_data %>% 
               mutate(prior = dbeta(p_vals, 1, 1),
                      likelihood = dbinom(x=y, size = 1e3, prob = p_vals),
                      unnorm_post = likelihood * prior,
                      posterior = unnorm_post / sum(unnorm_post),
                      theta = log(p_vals/(1-p_vals)),
                      e_theta = exp(theta))
  
ggplot(grid_data, aes(x= pi_vals, y=posterior)) + 
  geom_point() +
  geom_segment(aes(x = pi_vals, xend = pi_vals, y = 0, yend = posterior))

(theta_mean_grid <- mean(grid_data$theta))
(exp(theta_mean_grid))
(theta_var_grid <- var(grid_data$theta))
(const_pi_grid <- grid_data$unnorm_post |> sum())
quantile(prob = c(0.05, 0.95), grid_data$theta)




## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## (b) Gauss-Hermite rules
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
thetatilde <- theta_hat
VAR <- var_theta_hat
like_prior_ <- function(prob) {
  sum(rbinom(N, size=y, prob = prob))
}


# Constant 

myGH <- function(m=5, VAR, thetatilde){

   constant_integral <- 0
   quad <- statmod::gauss.quad(m,kind="hermite")
   for (i in 1:m){
     theta <- thetatilde + sqrt(2*VAR)*quad$nodes[i] 
     xi <- (exp(theta)/(1+exp(theta)))*(delta + gamma -1)+(1-gamma)
     print(xi)
     constant_integral <- constant_integral + 
                          quad$weights[i]*sqrt(2*VAR) *
                          like_prior_(prob=xi)
                         #sum(rbinom(N, size=y, prob = xi)) 
                         #like_prior_pi(theta)/sqrt(pi)
       #(1/dnorm(theta,thetatilde,sqrt(VAR))) * 
   }
   (constant_integral <- constant_integral)
   
   
   # Posterior Mean [first moment] EX
   EX <- 0
   quad <- statmod::gauss.quad(m,kind="hermite")
   for (i in 1:m){
     theta <- thetatilde + sqrt(2*var_pi_hat)*quad$nodes[i]
     xi <- (exp(theta)/(1+exp(theta)))*(delta + gamma -1)+(1-gamma)
     EX <- EX + 
           #quad$weights[i]*(1/dnorm(theta,thetatilde,sqrt(VAR))) * 
           quad$weights[i]*sqrt(2*VAR) *
           xi*like_prior_(xi)/sqrt(pi)
   }
   (p_EX <- EX/constant_integral)
   pi_EX <- p_EX*(delta+gamma-1) + (1-gamma)
   I_pi_EX <- (N*(delta + gamma -1)^2) / p_EX*(1-p_EX)
 
   theta_mean_GH <- log(pi_EX/(1-pi_EX))
   theta_var_GH <- 1/(I_pi_EX * pi_EX^2*(1-pi_EX^2))
   
   
   # Posterior Variance [use second moment] Var = EX^2 - (EX)^2
   # EX2 <- 0
   # quad <- statmod::gauss.quad(m,kind="hermite")
   # for (i in 1:m){
   #   theta <- thetatilde + sqrt(2*VAR)*quad$nodes[i]
   #   xi <- (exp(theta)/(1+exp(theta)))*(delta + gamma -1)+(1-gamma)
   #   EX2 <- EX2 + 
   #          #quad$weights[i]*(1/dnorm(theta,thetatilde,sqrt(VAR))) * 
   #          quad$weights[i]*sqrt(2*VAR) *
   #          xi^2*like_prior_(xi)/sqrt(pi)
   # }
   # 
   # theta_var_GH <- pi_EX^2*(1-pi_EX^2)
   # 
   # (p_EX2 <- EX2/constant_integral)
   # pi_EX2 <- p_EX2*(delta+gamma-1) + (1-gamma)
   # theta_var_GH <- log(pi_EX2/(1-pi_EX2)) - theta_mean_GH^2
   
   r <- list("constant" = constant_integral,
             "theta_mean_GH" = theta_mean_GH,
             "theta_var_GH" = theta_var_GH)
   return(r)
   
}

(GH_5 <- myGH(m=5, VAR=var_pi_hat, thetatilde = theta_hat))
(GH_4 <- myGH(m=4, VAR=var_pi_hat, thetatilde = theta_hat))







## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## (c)  Rejection Algo
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#pi_hat <- exp(theta_hat)/(1+exp(theta_hat))
mylike <- function(theta) {
  const <- choose(N, y)
  pi_ <- exp(theta)/(1+exp(theta))
  p_ <- pi_*(delta+gamma-1) + (1-gamma)
  #dbinom(length(theta), size=y, prob=p_)
  (p_)^y * (1-(p_))^(N-y)
}  

myloglike <- function(pi_) {
  #pi_ <- exp(theta)/(1+exp(theta))
  p_ <- pi_*(delta+gamma-1) + (1-gamma)
  theta <- log(pi_/(1-pi_))
  #log(dbinom(length(theta), size=y, prob=p_))
  log( (p_)^y * (1-p_)^(N-y) )
}  


logM <- myloglike(pi_hat)
accept <- 0
counter <- 0
n <- 5000 # No of samples to generate from the posterior
pi_samp <- rep(0,n)
unnorm <- rep(0,n)

while (accept < n) {
  counter <- counter + 1
  proposal_pi <- rbeta(1, 1, 1)
  u <- runif(1)
  test <- myloglike(proposal_pi) - logM
  if (log(u) < test) {
    accept <- accept + 1
    pi_samp[accept] <- proposal_pi
    pprop <- proposal_pi*(delta + gamma -1)+(1-gamma)
    unnorm[accept] <- sum(rbinom(N, size=y, prob=pprop) * pprop) 
  }
}
thetasamp <- log(pi_samp / (1-pi_samp))
(const_rej <- unnorm |> mean())
(theta_mean_rej <- mean(thetasamp))
(theta_var_rej <- var(thetasamp))
(CI_rej <- quantile(prob = c(0.05, 0.95), thetasamp))

par(mfrow=c(2,1))
hist(thetasamp, main = paste0("posteriror distribution of ", expression(theta)))
hist(pi_samp, main = paste0("posterior distribution of ", expression(pi)))





## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## (d) Importance Samp
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
n <- 5000 # pts to use to evaluate integrals
thsampn <- rt(n, df=4) 
thsampn <- thsampn / sd(thsampn)
thsampn <- thsampn * sqrt(var_pi_hat) + pi_hat
#thsampn <- thsampn * sqrt(var_theta_hat) + theta_hat

myloglike <- function(pi_) {
  #pi_ <- exp(theta)/(1+exp(theta))
  p_ <- pi_*(delta+gamma-1) + (1-gamma)
  #theta <- log(pi_/(1-pi_))
  log( choose(N,y) * (p_)^y * (1-p_)^(N-y) )
}  

mylognlike <- function(x) {
  # give me a theta, I give you pi:
  #pi_ <- exp(x)/(1+exp(x))
  
  # give me a pi, and I give you a p:
  #p_ <- x*(delta+gamma-1) + (1-gamma)
  
  #dmyt
  dt <- rep(0, length(x))
  for(i in 1:length(x)) {
    dt[i] <- log(mean(thsampn <= x[i] +.001 & thsampn >= x[i]-.001))
  }
  o <- myloglike(x) - dt
  return(o)
  
}  

# const
n0vals <- exp(mylognlike(thsampn))
(const_Imp     <- mean(n0vals))
(const_Imp_var <- var(n0vals)/n)
(const_Imp_CI  <- quantile(p=c(.05,.95), n0vals))

n0vals_pi <- n0vals*(delta+gamma-1) + (1-gamma)
n0vals_theta <- log(n0vals_pi/(1-n0vals_pi))


# EX
n1vals <- thsampn*exp(myloglike(thsampn)) / n0vals
n1vals_pi <- n1vals*(delta+gamma-1) + (1-gamma)
n1vals_theta <- log(n1vals_pi/(1-n1vals_pi))

# Posterior EX and Var
(p_EX_Imp     <- mean(n1vals_theta))
(p_EX_Imp_var <- var(n1vals_theta))
(p_EX_Imp_CI  <- quantile(p=c(.05,.95), n1vals_theta))
theta_mean_Imp <- p_EX_Imp
theta_var_Imp <- p_EX_Imp_var
CI_Imp <- p_EX_Imp_CI

# Var
#pi_EX_imp <- mean(n1vals_pi)
#I_pi_EX_Imp <- (N*(delta + gamma -1)^2) / pi_EX_Imp*(1-pi_EX_Imp)
#theta_var_Imp <- 1/(I_pi_EX_Imp * pi_EX_Imp^2*(1-pi_EX_Imp^2))




## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## (e) Metropolis Hastings
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
gfun <- function(theta, N=3330, y=50, delta=0.8, gamma=0.995){
  const <- choose(N, y)
  pi_ = exp(theta)/(1+exp(theta))
  y*log(pi_*(delta+gamma-1) + (1-gamma)) + (N-y)*log(1-(pi_*(delta+gamma-1) + (1-gamma)))
}

n <- 51000 # no if iterations
thetaMH <- NULL
pi_MH <- rep(0, 51000)
thetahat <- theta_hat #pi_hat = exp(theta_hat)/(1+exp(theta_hat))
sethetahat <- sqrt(var_theta_hat) #sqrt(  var_pi_hat )
thetaMH[1] <- thetahat
accept <- 0
for (i in 2:n){
  thetaMH[i] <- thetaMH[i-1]
  thetaprop <- rnorm(1,m=thetaMH[i-1],s=sqrt(3)*sethetahat)
  if (log(runif(1)) < gfun(theta=thetaprop) - gfun(theta=thetaMH[i-1])) {
    thetaMH[i] <- thetaprop 
    accept <- accept+1
  }
}
burnin <- 1000; indcalc <- seq(burnin+1,n)
# Estimate se's using batching
cat("MH:\n")
cat("Acceptance prob = ",accept/n,"\n")
B <- 1000 # length of each batch 
thetaMH <- thetaMH[indcalc]
(theta_mean_MH <- mean(thetaMH))
(theta_var_MH  <- var(thetaMH))
(CI_MH <- quantile(prob = c(0.05, 0.95), thetaMH))
pi_samp_MH <- exp(thetaMH)/(1+exp(thetaMH))

par(mfrow=c(2,1))
hist(thetaMH, main = paste0("posteriror distribution of ", expression(theta)))
hist(pi_samp_MH, main = paste0("posterior distribution of ", expression(pi)))



# TABLE OF RESULTS

tab_res <- data.frame("Mean"       = c(GH_4$theta_mean_GH, GH_5$theta_mean_GH, theta_mean_rej, theta_mean_Imp, theta_mean_MH),
                      "`Low 90% CI`"         = c(NA,                 NA,                 CI_rej[1], CI_Imp[1], CI_MH[1]),
                      "`Upp 90% CI`"         = c(NA,                 NA,                 CI_rej[2], CI_Imp[2], CI_MH[2]),
                      "Variance"        = c(GH_4$theta_var_GH,  GH_5$theta_var_GH,  theta_var_rej,  theta_var_Imp,  theta_var_MH),
                      "`Normilizing Constant`" = c(GH_4$constant,      GH_5$constant,      const_rej,      NA, NA))

colnames(tab_res) <- c("Mean", "Low 90% CI", "Upp 90% CI", "Variance", "Normilizing Constant")
row.names(tab_res) <- c("Gauss-Hermit m=4", "Gauss-Hermit m=5", "Rejection Algo", "Importance Samp", "Metropolis Algo")
tab_res


## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## (f)  Rejection Algo but with delta and gamma priors
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
myloglike <- function(delta, gamma, pi_) {
  #pi_ <- exp(theta)/(1+exp(theta))
  p_ <- pi_*(delta+gamma-1) + (1-gamma)
  theta <- log(pi_/(1-pi_))
  #log(dbinom(length(theta), size=y, prob=p_))
  log( (p_)^y * (1-p_)^(N-y) )
}  

x <- seq(from=0.0001, to=0.9999, length.out=5e3)

x <- pi_hat

# optim to find MLE for d and g
loglike <- function(b){
  log(choose(N,y)) + 
    y*log(x*(b[1]+b[2]-1)+(1-b[2])) + 
    (N-y)*log(1-(x*(b[1]+b[2]-1)+(1-b[2])))
}
my_optim <- optim(par = c(0.99,0.99), fn = loglike, control = list(fnscale = -1))

delta_hat <- my_optim$par[1]
gamma_hat <- my_optim$par[2]


logM <- myloglike(pi_hat, gamma_hat, pi_hat)
accept <- 0
counter <- 0
n <- 500 # No of samples to generate from the posterior
#pi_samp <- rep(0,n)
d_samp <- rep(0,n)
g_samp <- rep(0,n)
unnorm <- rep(0,n)

while (accept < n) {
  counter <- counter + 1
  #proposal_pi <- rbeta(1, 1, 1)
  proposal_delta <- rbeta(1, 160, 40)
  proposal_gamma <- rbeta(1, 1990, 10)
  u <- runif(1)
  test <- myloglike(pi_hat, proposal_delta, proposal_gamma) - logM
  if (log(u) < test) {
    accept <- accept + 1
    #pi_samp[accept] <- proposal_pi
     d_samp[accept] <- proposal_delta
     g_samp[accept] <- proposal_gamma
    #pprop <- proposal_pi*(delta + gamma -1)+(1-gamma)
    #unnorm[accept] <- sum(rbinom(N, size=y, prob=pprop) * pprop) 
  }
}

param_data_set <- expand.grid(pi_hat, d_samp, g_samp)
colnames(param_data_set) <- c("pi","d","g")

param_data_set$new_p <- param_data_set$pi*(param_data_set$d + param_data_set$g -1) + (1+param_data_set$g)
param_data_set$new_pi <-(param_data_set$p - (1-param_data_set$g))/(param_data_set$d + param_data_set$g -1) 
param_data_set$theta <- log(param_data_set$new_pi/(1-param_data_set$new_pi))

pi_samp_fun <- function(pi_hat, d_samp, g_samp) {
  c = 1
  sample <- rep(0, length(g_samp)^3)
  for(j in 1:length(d_samp)) {
    for(k in 1:g_samp) {
      sample[c] <- pi_hat*(d_samp[j] + g_samp[k] -1) + (1 - g_samp[k])
      c <- c + 1
    }
  }
  sample
} 

my_pi_samp <- pi_samp_fun(pi_hat, d_samp, g_samp)
  
my_thetasamp <- log(my_pi_samp / (1-my_pi_samp))
(const_rej <- unnorm |> mean())
(theta_mean_rej <- mean(my_thetasamp))
(theta_var_rej <- var(my_thetasamp))
(CI_rej <- quantile(prob = c(0.05, 0.95), my_thetasamp))

par(mfrow=c(2,1))
hist(param_data_set$theta, main = paste0("posteriror distribution of ", expression(theta)))
hist(param_data_set$new_pi, main = paste0("posterior distribution of ", expression(pi)))

## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++




save.image(file = "/Users/raulta/Library/CloudStorage/GoogleDrive-rauldavidt@gmail.com/My Drive/Stat 570/hw/midterm/midterm_objects.Rdata")
