#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# student: Raul
# course: STAT 570
# assignment: final
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(tidyverse)
library(dplyr)
rm(list=ls())

#~~~~~~~~~~~#
# Problem 1 #
#~~~~~~~~~~~#

smokers    <- c(29,16,17,4,3,9,4,5,1,1,1,3,7)
nonsmokers <- c(198,107,55,38,18,22,7,9,5,3,6,6,12)
stopifnot(length(smokers) == length(nonsmokers))

#~~~~~#
# 1.c #
#~~~~~#
get_phat <- function(cnts, t) {
  l <- length(cnts)
  stopifnot(l == length(t))
  num <- sum(cnts[1:l-1])
  denom <- cnts[l]*t[l] + sum(cnts[1:l-1]*(t[1:l-1]-1)) + num
  #num <- sum(cnts[l-1]*(t[l-1] - 1))
  #denom <- (2*num + cnts[l]*t[l])
  return(num/denom)
}

get_I <- function(cnts, phat) {
  #(1/phat^3 + 1/(phat*(1-phat)^2))*sum(cnts)
  #(1/(phat*(1-phat)^2))*sum(cnts[1:length(cnts)])
  sum(cnts[1:length(cnts)-1])/phat^2 + 
    length(cnts)*(cnts[length(cnts)])/(1-phat)^2 + 
      sum(cnts[1:length(cnts)-1])/(phat*(1-phat))
}

phat <- get_phat(cnts = smokers + nonsmokers, t = 1:13)
I <- get_I(smokers + nonsmokers, phat = phat)
var <- 1/I

#~~~~~#
# 1.d #
#~~~~~#
phat_smokers <- get_phat(cnts = smokers, t = 1:13)
var_smokers <- 1/get_I(smokers, phat_smokers)
phat_nonsmokers <- get_phat(cnts = nonsmokers, t = 1:13)
var_nonsmokers <- 1/get_I(nonsmokers, phat_nonsmokers)

get_ci <- function(theta, var, alpha, n) {
  t <- qnorm(1-alpha/2)
  return(c(theta - t*sqrt(var/n), theta + t*sqrt(var/n)))
}

ci_smokers    <- get_ci(phat_smokers,    var_smokers,    0.05, n=13)
ci_nonsmokers <- get_ci(phat_nonsmokers, var_nonsmokers, 0.05, n=13)


#~~~~~#
# 1.g #
#~~~~~#
mu = 0.2
sigma = 0.08
a <- mu*(mu*(1-mu)/sigma^2 - 1)
b <- a*(1/mu-1) 
test <- rbeta(1e4, shape1=a, shape2=b)
stopifnot(mean(test) - mu < 0.005)
stopifnot(sd(test) - sigma < 0.002)
rm(test)


get_Yt         <- function(x)   sum(x[1:length(x)-1])
get_Yt_tminus1 <- function(x,t) sum(x[1:length(x)-1]*t[1:length(t)-1])
NY_Nplus1      <- function(x)   13*x[length(x)]
get_post_mean <- function(a,b) a/(a+b)


post_mean_smokers <- get_post_mean(a=get_Yt(smokers) -a -1,
                                   b=get_Yt_tminus1(smokers, 1:13) -b -1)

post_mean_nonsmokers <- get_post_mean(a=get_Yt(nonsmokers) -a -1,
                                      b=get_Yt_tminus1(nonsmokers, 1:13) -b -1)

get_cred_int <- function(a, b, ndraws, pc) {
  post <- rbeta(ndraws, a, b)
  cred_int <- quantile(post, probs = c((1-pc)/2, (1+pc)/2))
  return(cred_int)
}

cred_int_smokers <- get_cred_int(a = get_Yt(smokers) -a -1,
                                 b = get_Yt_tminus1(smokers, 1:13) -b -1,
                                 ndraws = 1e4,
                                 pc = 0.95)

cred_int_nonsmokers <- get_cred_int(a = get_Yt(nonsmokers) -a -1,
                                    b = get_Yt_tminus1(nonsmokers, 1:13) -b -1,
                                    ndraws = 1e4,
                                    pc = 0.95)



q1g <- data.frame(cred_int_025 = c(cred_int_smokers[[1]], cred_int_nonsmokers[[1]]),
                  post_mean    = c(post_mean_smokers, post_mean_nonsmokers),
                  cred_int_975 = c(cred_int_smokers[[2]], cred_int_nonsmokers[[2]]))

row.names(q1g) <- c("smokers", "non-smokers")
colnames(q1g) <- c("2.5%","mean","97.5%")


#~~~~~#
# 2.d #
#~~~~~#

#mylikelihood <- function(theta, x=smokers) {
#  alpha <- theta[1]
#  beta <- theta[2]
#  p <- rbeta(1, alpha, beta)
#  (p*(1-p))^sum(x[1:12])*(1-p)^(13*x[13])
#}

mylikelihood <- function(theta, x) {
  alpha <- theta[1]; beta <- theta[2]; p<-theta[3]
  t <- 1:13
  sum(x[1:12])*lbeta(alpha, beta) + 
    sum(x[1:12]*(alpha-1))*log(p) +
    sum(x[1:12]*(beta -1))*log(1-p) + 
    sum(x[1:12]*(t[1:12]-1)*13*x[13]) * log(1-base::beta(alpha,beta)*(p)^(alpha-1)*(1-p)^(beta-1))
}

beta_params_smokers <- optim(par = c(2,5,0.5),
                             fn = mylikelihood, 
                             control = list(fnscale = -1),
                             x=smokers)$par

beta_params_nonsmokers <- optim(par = c(2,5,0.5),
                                fn = mylikelihood, 
                                control = list(fnscale = -1),
                                x=nonsmokers)$par

p_smokers_2d    <- mean(rbeta(1e5, beta_params_smokers[1]   , beta_params_smokers[2])) 
p_nonsmokers_2d <- mean(rbeta(1e5, beta_params_nonsmokers[1], beta_params_nonsmokers[2]))  


#~~~~~#
# 3.b #
#~~~~~#

get_phat_binoms <- function(cnts, t) {
  l <- length(cnts)
  stopifnot(l == length(t))
  num <- sum(cnts[1:l-1])
  denom <- num + sum(cnts[2:l]*(t[2:l]-1)) + cnts[l]
  return(num/denom)
}

phat_binoms_smokers <- get_phat_binoms(smokers, 1:13)
phat_binoms_nonsmokers <- get_phat_binoms(nonsmokers, 1:13)

get_I_binoms <- function(cnts,p) {
    I <- 1/(p*(1-p)) * cnts
    return(sum(I))
}
var_binoms_smokers    <- 1/get_I_binoms(smokers, phat_binoms_smokers)
var_binoms_nonsmokers <- 1/get_I_binoms(nonsmokers, phat_binoms_nonsmokers)

ci_binoms_smokers    <- get_ci(phat_binoms_smokers, var_binoms_smokers, 0.05, n=13)
ci_binoms_nonsmokers <- get_ci(phat_binoms_nonsmokers, var_binoms_nonsmokers, 0.05, n=13)



tab3b <- data.frame(geometric=c(phat_smokers, phat_nonsmokers),
                    binomial=c(phat_binoms_smokers, phat_binoms_nonsmokers),
                    geometric_c1 = c(ci_smokers[1], ci_nonsmokers[1]),
                    geometric_c2 = c(ci_smokers[2], ci_nonsmokers[2]),
                    binomial_c1 = c(ci_binoms_smokers[1], ci_binoms_nonsmokers[1]),
                    binomial_c2 = c(ci_binoms_smokers[2], ci_binoms_nonsmokers[2])
                    )
colnames(tab3b) <- c("Geom","Binom", "Geom 2.5%", "Geom 97.5%", "Binom 2.5%", "Binom 97.5%")
row.names(tab3b) <- c("smokers","non-smokers")


#~~~~~#
# 3.c #
#~~~~~#

# Grid method 
pi_vals <- seq(from = 0, to = 1, length = 5e2)

binom_likelihood <- function(x, pi_vals) {
  l <- length(x)
  N_sizes <- rep(0, l)
  for(i in 1:l) {
    N_sizes[i] <- sum(x[i:l])
  }
  like_val <- rep(0, length(pi_vals))
  i=0
  for(p in pi_vals){
    i <- i+1
    like_val[i] <- prod(dbinom(x, size=N_sizes, prob = p))
  }
  return(like_val)
}

get_posterior_df <- function(pi_vals, x) {
  grid_data <- data.frame(pi_vals = pi_vals)
  grid_data <- grid_data %>% 
    mutate(prior = dbeta(pi_vals, 1, 1),
           likelihood = binom_likelihood(x=x, pi_vals=pi_vals),
           unnorm_post = likelihood * prior,
           posterior = unnorm_post / sum(unnorm_post))
}

posterior_smokers <- get_posterior_df(pi_vals = pi_vals, x=smokers)

plot_p1 <- ggplot(posterior_smokers, aes(x= pi_vals, y=posterior)) + 
  geom_point() + 
  xlab("suggested p") + ylab("posterior probability") + 
  ggtitle("smokers") + 
  geom_segment(aes(x = pi_vals, xend = pi_vals, y = 0, yend = posterior)) + 
  theme(panel.grid=element_blank(),
        panel.background = element_blank())

posterior_nonsmokers <- get_posterior_df(pi_vals = pi_vals, x=nonsmokers)

plot_p2 <- ggplot(posterior_nonsmokers, aes(x= pi_vals, y=posterior)) + 
  geom_point() + 
  xlab("suggested p") + ylab("posterior probability") + 
  ggtitle("nonsmokers") + 
  geom_segment(aes(x = pi_vals, xend = pi_vals, y = 0, yend = posterior)) + 
  theme(panel.grid=element_blank(),
        panel.background = element_blank())

posterior_diff <- posterior_nonsmokers$pi_vals*(posterior_nonsmokers$posterior - 
                                                posterior_smokers$posterior)
hist(posterior_diff, 
     main="difference in probability between\n nonsmokers smokers",
     xlab = "posterior difference",
     ylab = "",
     freq=FALSE)

post_mean_nonsmokers_p3 <- posterior_nonsmokers$pi_vals[posterior_nonsmokers$posterior==max(posterior_nonsmokers$posterior)]

post_mean_smokers_p3 <- posterior_smokers$pi_vals[posterior_smokers$posterior==max(posterior_smokers$posterior)]

p3ctab <- data.frame(bayesian1 = c(post_mean_smokers, post_mean_nonsmokers),
                     bayesian2 = c(post_mean_smokers_p3, post_mean_nonsmokers_p3))

colnames(p3ctab) <- c("mean from problem 2","mode from problem 3")
row.names(p3ctab) <- c("smokers","non-smokers")


#~~~~~#
# 3.d #
#~~~~~#
prob_p2gtp1 <- mean(posterior_nonsmokers$pi_vals*posterior_nonsmokers$posterior > 
                    posterior_smokers$pi_vals*posterior_smokers$posterior)


#~~~~~#
# 3.e #
#~~~~~#

# simulate 13 periods with successes with prob post
get_sample <- function(N=100, t=1:13, post_p, probs=FALSE) {
  Y_t <- rep(0,length(t))
  Ns <- rep(0,length(t))
  N <- N
  for(i in t) {
    Ns[i] <- N
    Y_t[i] <- rbinom(1, size=N, prob=post_p)
    N <- N - Y_t[i]
  }
  if(probs==FALSE){
     return(Y_t)
  }else{
    Y_t_probs = Y_t/Ns
    return(Y_t_probs)
  }
}

#get_sample(N=100, post_p=post_mean_nonsmokers[1], probs=TRUE)

# weight each sample by posterior weight of suggested p
get_pred_dist <- function(Nsims=1, wts, pis) {
  predictive_dist_df <- data.frame(matrix(ncol = 15, nrow = 0))
  for(p in 1:length(pis)) {
    for(j in 1:Nsims) {
      row <- c(pis[p],
               wts[p],
               get_sample(N=100, t=1:13, pis[p])
      )
      predictive_dist_df <- rbind(predictive_dist_df, row)
    }
  }
  colnames(predictive_dist_df) <- c("Ps","pred_wt",paste("t", 1:13,sep=""))
  return(predictive_dist_df)
}

few <- seq(from=5, to=length(pi_vals), by=2)
pi_weights_smokers <- posterior_smokers$posterior
pi_weights_nonsmokers <- posterior_nonsmokers$posterior
pred_dist_smokers_df <- get_pred_dist(pis=pi_vals, wts=pi_weights_smokers)
pred_dist_nonsmokers_df <- get_pred_dist(pis=pi_vals, wts=pi_weights_nonsmokers)

par(mfrow = c(1,2))
plotrix::weighted.hist(pred_dist_smokers_df$t1,pred_dist_smokers_df$pred_wt, main="t1")
plotrix::weighted.hist(pred_dist_smokers_df$t2,pred_dist_smokers_df$pred_wt, main="t2")
plotrix::weighted.hist(pred_dist_smokers_df$t3,pred_dist_smokers_df$pred_wt, main="t3")
plotrix::weighted.hist(pred_dist_smokers_df$t4,pred_dist_smokers_df$pred_wt, main="t4")
par(mfrow = c(2,2))
plotrix::weighted.hist(pred_dist_smokers_df$t5,pred_dist_smokers_df$pred_wt, main="t5")
plotrix::weighted.hist(pred_dist_smokers_df$t8,pred_dist_smokers_df$pred_wt, main="t8")
plotrix::weighted.hist(pred_dist_smokers_df$t10,pred_dist_smokers_df$pred_wt, main="t10")
plotrix::weighted.hist(pred_dist_smokers_df$t13,pred_dist_smokers_df$pred_wt, main="t13")

par(mfrow = c(2,2))
plotrix::weighted.hist(pred_dist_nonsmokers_df$t1,pred_dist_nonsmokers_df$pred_wt, main="t1")
plotrix::weighted.hist(pred_dist_nonsmokers_df$t2,pred_dist_nonsmokers_df$pred_wt, main="t2")
plotrix::weighted.hist(pred_dist_nonsmokers_df$t3,pred_dist_nonsmokers_df$pred_wt, main="t3")
plotrix::weighted.hist(pred_dist_nonsmokers_df$t4,pred_dist_nonsmokers_df$pred_wt, main="t4")
par(mfrow = c(2,2))
plotrix::weighted.hist(pred_dist_nonsmokers_df$t5,pred_dist_nonsmokers_df$pred_wt, main="t5")
plotrix::weighted.hist(pred_dist_nonsmokers_df$t8,pred_dist_nonsmokers_df$pred_wt, main="t8")
plotrix::weighted.hist(pred_dist_nonsmokers_df$t10,pred_dist_nonsmokers_df$pred_wt, main="t10")
plotrix::weighted.hist(pred_dist_nonsmokers_df$t13,pred_dist_nonsmokers_df$pred_wt, main="t13")


get_pred_dist_probs <- function(Nsims=1, wts, pis) {
  predictive_dist_df <- data.frame(matrix(ncol = 15, nrow = 0))
  for(p in 1:length(pis)) {
    for(j in 1:Nsims) {
      row <- c(pis[p],
               wts[p],
               get_sample(N=100, t=1:13, pis[p], probs = TRUE)
      )
      predictive_dist_df <- rbind(predictive_dist_df, row)
    }
  }
  colnames(predictive_dist_df) <- c("Ps","pred_wt",paste("t", 1:13,sep=""))
  return(predictive_dist_df)
}

pred_dist_probs_smokers_df <- get_pred_dist_probs(pis=pi_vals, wts=pi_weights_smokers)
pred_dist_probs_nonsmokers_df <- get_pred_dist_probs(pis=pi_vals, wts=pi_weights_nonsmokers)

tab3e <- 
  purrr::map_dbl(pred_dist_probs_smokers_df, mean, na.rm=TRUE) %>% 
  rbind(purrr::map_dbl(pred_dist_probs_nonsmokers_df, mean, na.rm=TRUE)) %>% 
  as.data.frame()

tab3e <- tab3e[,3:ncol(tab3e)]
row.names(tab3e) <- c("smokers","non-smokers")


save.image(file = "/Users/raulta/Library/CloudStorage/GoogleDrive-rauldavidt@gmail.com/My Drive/Stat 570/hw/zfinal/final_objects.Rdata")










