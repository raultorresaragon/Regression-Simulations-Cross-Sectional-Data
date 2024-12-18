# Raul Torres Aragon
# Stat 570
# Date: 2020/11/01
# Assignment 5


# ~~~~~~~~~~ #
# QUESTION 1 #
# ~~~~~~~~~~ #

# (b)
p1 <- (96/(96+109))
p2 <- (104/(104+666))
n1 <- 200
n2 <- 775
num <- (96/(96+109)) / (109/(96+109))
den <- (104/(104+666)) / (666/(104+666))
MLE <- num/den

FIM = matrix(c( n1/(p1*(1-p1)), 0, 0, n2/(p2*(1-p2))), nrow=2, ncol=2)

Grad <- matrix(c(1/(p1*(1-p1)),  1/(p2*(1-p2))), nrow=1, ncol=2)

VarLogTheta <- Grad %*% solve(FIM) %*% t(Grad)

CI_MLE <- c(exp(log(MLE)-qnorm(0.95,mean=0,sd=1)*sqrt(VarLogTheta)),
            exp(log(MLE)+qnorm(0.95,mean=0,sd=1)*sqrt(VarLogTheta)))

# (d)

a1 <- 97; b1 <- 105
a2 <- 109; b2 <- 667

get_beta_stuff <- function(a, b) {
  EX <- a/(a+b)
  Var <- (a*b) / ((a+b+1)*(a+b)^2)
  Mode <- (a-1)/(a+b-2)
  
  r <- list("EX" = EX, "Var" = Var, "Mode" = Mode)
  return(r)
}

p1_stuff <- get_beta_stuff(a1, b1)
p2_stuff <- get_beta_stuff(a2, b2)


CI_1 <- c(qbeta(p = .05, shape1 = 97, shape2 = 105),
          qbeta(p = .95, shape1 = 97, shape2 = 105))
CI_2 <- c(qbeta(p = .05, shape1 = 109, shape2 = 667),
          qbeta(p = .95, shape1 = 109, shape2 = 667))

# (e)
# 
p1_dist <- rbeta(1e3, shape1 = 1, shape2 = 1)
p2_dist <- rbeta(1e3, shape1 = 1, shape2 = 1)
prior_dist <- p1_dist + p2_dist
C_I_prior <- quantile(prior_dist, probs=c(0.05,0.95))

hist(prior_dist, main = "A thousand simulated draws \n from prior(p1,p2)",
     xlab = "prior distribution")
abline(v = C_I_prior[[1]], col="black", lwd=3, lty=1)
abline(v = C_I_prior[[2]], col="black", lwd=3, lty=1)




# (f)
n <- 1e4
p1s <- rbeta(n=n, a1, b1)
p2s <- rbeta(n=n, a2, b2)

CI_1_from_hist <- c(quantile(p1s, probs = 0.05), quantile(p1s, probs = 0.95))
CI_2_from_hist <- c(quantile(p2s, probs = 0.05), quantile(p2s, probs = 0.95))

par(mfrow=c(1,2))

hist(p1s, main = "1000 draws from P(p1|X1)") 
abline(v = CI_1_from_hist[[1]], col="black", lwd=3, lty=1)
abline(v = CI_1_from_hist[[2]], col="black", lwd=3, lty=1)

hist(p2s, main = "1000 draws from P(p2|X2)")
abline(v = CI_2_from_hist[[1]], col="black", lwd=3, lty=1)
abline(v = CI_2_from_hist[[2]], col="black", lwd=3, lty=1)




# (g)

theta <- function(p1, p2) {
  num <- p1/(1-p1)
  denom <- p2/(1-p2)
  return(num/denom)
}

thetas <- theta(p1s, p2s)
CI_theta <- c(quantile(thetas, 0.05), quantile(thetas, 0.95))

hist(thetas, main = "histogram of 1,000 simulated thetas")
abline(v = CI_theta[[1]], col="black", lwd=3, lty=1)
abline(v = CI_theta[[2]], col="black", lwd=3, lty=1)



# (h)
PY1 = 18/100000
PY0 = 1-PY1
PY1gX1 = ((96/200)*PY1)/(((96/200)*PY1)+((109/775)*PY0))
PY1gX0 = ((104/200)*PY1) / ((104/200)*PY1 + (666/775)*PY0)


# (i)
# numerically find a and b
myfun <- function(theta, low=16/1e5, upp=20/1e5) {
    (0.05 - pbeta(low, theta[1], theta[2]))^2 + 
    (0.95 - pbeta(upp, theta[1], theta[2]))^2
}

beta_params <- optim(par = c(1,1),
                      fn = myfun,
                      control = list(abstol = 1/1e7))$par

qbeta(0.95, shape1 = beta_params[1], shape2 = beta_params[2])
qbeta(0.05, shape1 = beta_params[1], shape2 = beta_params[2])

alpha <- beta_params[1]
beta <- beta_params[2]

py1s <- rbeta(n=1000, alpha, beta)
py0s <- 1-py1s

PY1gX1s <- ((96/200)*py1s)/(((96/200)*py1s)+((109/775)*py0s))
median_q1 <- quantile(probs= 0.5, PY1gX1s)

PY1gX0s <- ((104/200)*py1s) / ((104/200)*py1s + (666/775)*py0s)
median_q0 <- quantile(probs= 0.5, PY1gX0s)


theta_q <- PY1gX1s/PY1gX0s


# ~~~~~~~~~~ #
# QUESTION 2 #
# ~~~~~~~~~~ #

# ~~~~~~~~~~ #
# QUESTION 3 #
# ~~~~~~~~~~ #
library("INLA")
library("dplyr")

dlung <- read.table(url("http://faculty.washington.edu/jonno/book/MNlung.txt"),
                    sep = "\t", header = TRUE) |> as_tibble()
drad  <- read.table(url("http://faculty.washington.edu/jonno/book/MNradon.txt"),
                    sep = " ", header = TRUE) |> as_tibble()

colnames(dlung) <- c("county", "county_n", 
                     "obs_male",  "exposure_male", 
                     "obs_female","exposure_female", 
                     "avg_male", "avg_female")

colnames(drad) <- c("county_n", "radon")

df <- dlung |> 
      left_join(drad |> 
                group_by(county_n) |> 
                summarise(mean_rad = mean(radon, na.rm = TRUE)),
      by = "county_n") |>
      mutate(obs_all = obs_male + obs_female)


## (a)


inla_fit <- inla(formula = obs_all ~ mean_rad, family = "poisson", data = df) #, 
#                 control.compute=list(config = TRUE))

inla.posterior.sample(n=100, result = inla_fit)

beta_0_marg <- inla_fit$marginals.fixed$`(Intercept)`[,1]
beta_1_marg <- inla_fit$marginals.fixed$mean_rad[,1]
inla_beta_results <- inla_fit$summary.fixed |> as_tibble() |> dplyr::select(-mode, -kld)
inla_post_results <- inla_fit$summary.fitted.values |> as_tibble()
#inla_post_samples <- inla.posterior.sample(n = 5, inla_fit2)


## (b)

# demean X
df$demeaned_rad <- df$mean_rad - mean(df$mean_rad)

# recall INLA takes in precision not sigma
inla_fit2 <- inla(formula = obs_all ~ demeaned_rad, 
                  family = "poisson", 
                  data = df,
                  control.fixed = list(mean.intercept=0, 
                                       prec.intercept=1/0.21^2,
                                       mean          = -0.02,
                                       prec          =list(param=1/0.10^2, prior = "lognormal"))
              )

beta_0_marg2 <- inla_fit2$marginals.fixed$`(Intercept)`[,1]
beta_1_marg2 <- inla_fit2$marginals.fixed$demeaned_rad[,1]
inla_beta_results2 <- inla_fit2$summary.fixed |> as_tibble() |> dplyr::select(-mode, -kld)
inla_post_results2 <- inla_fit2$summary.fitted.values |> as_tibble()

save.image(file = "/Users/raulta/Library/CloudStorage/GoogleDrive-rauldavidt@gmail.com/My Drive/Stat 570/hw/hw5/hw5_objects.Rdata")








