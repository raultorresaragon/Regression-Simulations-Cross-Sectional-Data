#~~~~~~~~~~~~~~~~~~~~~~~
# student: Raul Aragon
# class: STAT 570
# assignment: hw3
# notes:
#~~~~~~~~~~~~~~~~~~~~~~~

library(MASS)
library(stringr)
library(tidyverse)
library(tibble)
library(dplyr)
rm(list=ls())
set.seed(570)

#~~~~~~~~~~~~#
# Question 1 #
#~~~~~~~~~~~~#


# HELPER FUNCTIONS
# ----------------

# negative binomial estimator
my_nbin <- function(X, Y, b0, b1) {
  mod <- try(glm.nb(Y~X))
  if(class(mod) == "try-error") { 
    b0_inside <- NA
    b1_inside <- NA
  } else {
    mat <- coef(summary(mod))
    CIs <- confint.default(mod)
    b0_inside <- (CIs[1,1] < b0 & b0 < CIs[1,2]) |> as.numeric()
    b1_inside <- (CIs[2,1] < b1 & b1 < CIs[2,2]) |> as.numeric()
  }
  return(list("b0_in_nbin" = b0_inside,
              "b1_in_nbin" = b1_inside))

}

# poisson estimator using MASS::glm
my_pois <- function(X, Y, b0, b1) {
  mod <- try(glm(Y~X, family = poisson))
  if(class(mod) == "try-error") { 
    b0_inside <- NA
    b1_inside <- NA
  } else {
  mat <- coef(summary(mod))
  CIs <- confint.default(mod)
  b0_inside <- (CIs[1,1] < b0 & b0 < CIs[1,2]) |> as.numeric()
  b1_inside <- (CIs[2,1] < b1 & b1 < CIs[2,2]) |> as.numeric()
  }
  return(list("b0_in_pois" = b0_inside,
              "b1_in_pois" = b1_inside))  
}

# poisson quasi-likelihood using MASS::glm
my_quas <- function(X, Y, b0, b1) {
  mod <- glm(Y~X, family = quasipoisson) |> try()
  if(class(mod) == "try-error") { 
    b0_inside <- NA
    b1_inside <- NA
  } else {  
  mat <- coef(summary(mod))
  CIs <- confint.default(mod)
  b0_inside <- (CIs[1,1] < b0 & b0 < CIs[1,2]) |> as.numeric()
  b1_inside <- (CIs[2,1] < b1 & b1 < CIs[2,2]) |> as.numeric()
  }
  return(list("b0_in_quas" = b0_inside,
              "b1_in_quas" = b1_inside))  
}

# sandwich estimator using sandwich::sandwich
my_sand <- function(X, Y, b0, b1) {
  mod <- glm(Y~X, family = poisson) |> try()
  if(class(mod) == "try-error") { 
    b0_inside <- NA
    b1_inside <- NA
  } else {  
  sand <- sandwich::sandwich(mod)
  
  b0_inside <- ((mod$coefficients[[1]] - sqrt(sand[1,1])) < b0 &
                 b0 < mod$coefficients[[1]] + sqrt(sand[1,1])) |> as.numeric()
  b1_inside <- ((mod$coefficients[[2]] - sqrt(sand[2,2])) < b1 &
                 b1 < mod$coefficients[[2]] + sqrt(sand[2,2])) |> as.numeric()
  }
  return(list("b0_in_sand" = b0_inside,
              "b1_in_sand" = b1_inside))  
  
}
  
  

# simulations function
# --------------------

run_sims <- function(n, b, b0=0.5, b1=log(2.5), N = 1000) {
  
  b0_nbin_inside <- rep(NA, N) #vector(mode = "numeric", length = N)
  b1_nbin_inside <- rep(NA, N) #vector(mode = "numeric", length = N)
  b0_pois_inside <- rep(NA, N) #vector(mode = "numeric", length = N)
  b1_pois_inside <- rep(NA, N) #vector(mode = "numeric", length = N)
  b0_quas_inside <- rep(NA, N) #vector(mode = "numeric", length = N)
  b1_quas_inside <- rep(NA, N) #vector(mode = "numeric", length = N)
  b0_sand_inside <- rep(NA, N) #vector(mode = "numeric", length = N)
  b1_sand_inside <- rep(NA, N) #vector(mode = "numeric", length = N)
  
  j = 0
  for(i in 1:N) {
    j = j+1
    
    X <- rnorm(mean=0, sd=1, n=n)
    theta <- rgamma(shape=b, rate=b, n=n)
    mu <- exp(b0 + b1*X)
    Y <- rpois(lambda = theta*mu, n=n) 
    
    # negative binomial
    mod_nbin <- my_nbin(Y=Y, X=X, b0=b0, b1=b1) 
    b0_nbin_inside[j] <- mod_nbin[[1]]
    b1_nbin_inside[j] <- mod_nbin[[2]]
  
    # poisson
    mod_pois <- my_pois(Y=Y, X=X, b0=b0, b1=b1)
    b0_pois_inside[j] <- mod_pois[[1]]
    b1_pois_inside[j] <- mod_pois[[2]]
      
    # quasi-likelihood  
    mod_quas <- my_quas(Y=Y, X=X, b0=b0, b1=b1)
    b0_quas_inside[j] <- mod_quas[[1]]
    b1_quas_inside[j] <- mod_quas[[2]]
      
    # sandwich estimator  
    mod_sand <- my_sand(Y=Y, X=X, b0=b0, b1=b1)
    b0_sand_inside[j] <- mod_sand[[1]]
    b1_sand_inside[j] <- mod_sand[[2]]
  }
  
  o <- list("b0_nbin_inside" = mean(b0_nbin_inside),
            "b1_nbin_inside" = mean(b1_nbin_inside),
            "b0_pois_inside" = mean(b0_pois_inside),
            "b1_pois_inside" = mean(b1_pois_inside),
            "b0_quas_inside" = mean(b0_quas_inside),
            "b1_quas_inside" = mean(b1_quas_inside),
            "b0_sand_inside" = mean(b0_sand_inside),
            "b1_sand_inside" = mean(b1_sand_inside),
            "b0_vector_nbin" = b0_nbin_inside,
            "b1_vector_nbin" = b1_nbin_inside,
            "b0_vector_pois" = b0_pois_inside,
            "b1_vector_nbin" = b1_pois_inside,
            "b0_vector_quas" = b0_quas_inside,
            "b1_vector_quas" = b1_quas_inside,
            "b0_vector_sand" = b0_sand_inside,
            "b1_vector_sand" = b1_sand_inside)
  return(o)
}

ns <- c(10,20,50,100,250)
bs <- c(1, 10, 1000, 0.1)

# invoke simulations function for each combination of ns and bs

b0_table <- tibble("n" = integer(), 
                   "b" = integer(), 
                   "nbinomial" = double(),
                   "poisson" = double(),
                   "quasi" = double(),
                   "sandwich" = double()) 

b1_table <- tibble("n" = integer(), 
                   "b" = integer(), 
                   "nbinomial" = double(),
                   "poisson" = double(),
                   "quasi" = double(),
                   "sandwich" = double())
for(b in bs) {
  for(n in ns) {
    print(paste0("for n=", n, " and b=", b))
    vals <- run_sims(n=n, b=b, b0=0.5, b1=log(2.5), N=1e4)
    b0_table <- b0_table |> 
                add_row("n" = n, "b" = b, 
                        "nbinomial" = vals[[1]], 
                        "poisson" = vals[[3]], 
                        "quasi" = vals[[5]], 
                        "sandwich" = vals[[7]])
    
    b1_table <- b1_table |> 
                add_row("n" = n, "b" = b, 
                        "nbinomial" = vals[[2]], 
                        "poisson" = vals[[4]], 
                        "quasi" = vals[[6]], 
                        "sandwich" = vals[[8]])
  }
}

b1_tab <- pivot_wider(b1_table, id_cols = n, 
                      names_from = b, 
                      values_from = nbinomial:sandwich
          ) |> dplyr::select(n, 
                       `nbinomial_0.1`, `poisson_0.1`, `quasi_0.1`, `sandwich_0.1`,
                        nbinomial_1,    poisson_1,      quasi_1,     sandwich_1,
                        nbinomial_10,   poisson_10,     quasi_10,    sandwich_10,
                        nbinomial_1000, poisson_1000,   quasi_1000,  sandwich_1000)

b0_tab <- pivot_wider(b0_table, id_cols = n, 
                      names_from = b, 
                      values_from = nbinomial:sandwich
          ) |> dplyr::select(n, 
                   `nbinomial_0.1`, `poisson_0.1`, `quasi_0.1`, `sandwich_0.1`,
                   nbinomial_1,    poisson_1,      quasi_1,     sandwich_1,
                   nbinomial_10,   poisson_10,     quasi_10,    sandwich_10,
                   nbinomial_1000, poisson_1000,   quasi_1000,  sandwich_1000)



#~~~~~~~~~~~~#
# Question 2 #
#~~~~~~~~~~~~#


# (b)

g1<-c(2.247,2.640,2.842,2.908,3.099,3.126,3.245,3.328,3.355,3.383,3.572,3.581,3.681)
g2<-c(1.901,2.132,2.203,2.228,2.257,2.350,2.361,2.396,2.397,2.445,2.454,2.454,2.474)
g3<-c(1.312,1.314,1.479,1.552,1.700,1.803,1.861,1.865,1.944,1.958,1.966,1.997,2.006)
g4<-c(1.339,1.434,1.549,1.574,1.589,1.613,1.746,1.753,1.764,1.807,1.812,1.840,1.852)

get_lambda <- function(g) {
  lambda <- length(g)/sum(g)
  var <- (1/length(g)) * (lambda^2)
  s <- sqrt(var)
  return(list("lambda" = lambda, "var" = var, "s" = s))
}

g1_stuff <- get_lambda(g1)
g2_stuff <- get_lambda(g2)
g3_stuff <- get_lambda(g3)
g4_stuff <- get_lambda(g4)

var_lamb1 <- lambda_1^2*(1/13)

# exponential qq-plot 
gps <- list(g1, g2, g3 ,g4)
lambdas <- list(g1_stuff$lambda, g2_stuff$lambda, g3_stuff$lambda, g4_stuff$lambda)

par(mfrow=c(2,2))
for(i in 1:4) {
qqnorm(gps[[i]] - lambdas[[i]], pch = 1, 
       ylab = paste0("Group ", i, " - MLE(", i,")"), 
       main = paste0("Group ", i))
qqline(gps[[i]] - lambdas[[i]], col = "darkblue", lwd =2)
}

# exponential QQ plots
### par(mfrow=c(2,2))
### pts <- ppoints(13)
### q <- quantile(g1, p = pts)
### plot(qexp(pts), q, 
###      main="Exponential Q-Q Plot",
###      xlab="Theoretical Quantiles",
###      ylab="Group 1 Quantiles")
### qqline(q, distribution=qexp,col="blue", lty=2)


# (c)
# estimator for \alpha

alphas <- vector(mode = "numeric", length = 4)
for(i in 1:4) {  
 alphas[i] <- (lambdas[[i]]^2/(length(gps[[i]])-2)) * sum((gps[[i]] - (lambdas[[i]])^-1)^2)
}
alphas

g1_stuff$var_quas <- alphas[[1]]*(lambdas[[1]])^2
g1_stuff$s_quas <- sqrt(alphas[[1]]*(lambdas[[1]])^2)

g2_stuff$var_quas <- alphas[[2]]*(lambdas[[2]])^2
g2_stuff$s_quas <- sqrt(alphas[[2]]*(lambdas[[2]])^2)

g3_stuff$var_quas <- alphas[[3]]*(lambdas[[3]])^2
g3_stuff$s_quas <- sqrt(alphas[[3]]*(lambdas[[3]])^2)

g4_stuff$var_quas <- alphas[[4]]*(lambdas[[4]])^2
g4_stuff$s_quas <- sqrt(alphas[[4]]*(lambdas[[4]])^2)

# (d)
# 

G <- function(Y, lambda) {
  length(Y)/lambda - sum(Y)
}

A <- function(Y, lambda) {
   -length(Y) / lambda^2
}

B <- function(Y, lambda) {
  (1/length(Y)) * sum((length(Y)/lambda - sum(Y))^2)
}


Var_sand <- function(Y, lambda) {
  lambda^4/length(Y)^2 * sum(((1/lambda) - Y)^2)
}

g1_stuff[["s_sand"]] <- sqrt(Var_sand(gps[[1]], lambdas[[1]]))
g2_stuff[["s_sand"]] <- sqrt(Var_sand(gps[[2]], lambdas[[2]]))
g3_stuff[["s_sand"]] <- sqrt(Var_sand(gps[[3]], lambdas[[3]]))
g4_stuff[["s_sand"]] <- sqrt(Var_sand(gps[[4]], lambdas[[4]]))

glimpse(g1_stuff)
glimpse(g2_stuff)
glimpse(g3_stuff)
glimpse(g4_stuff)


# (h)
alpha_fun <- function(Y, eta) {
  ((1/length(Y)) * sum(Y^eta))^(1/eta)
}

my_weibul_fun <- function(Y) {

  Y <- Y
  eta_fun <- function(eta) {
    (1/eta[1]) - (sum(Y^eta[1] * log(Y)) / sum(Y^eta[1])) + (1/length(Y))*sum(log(Y))
  }
  eta <- optim(par = 2,
               fn = eta_fun, 
               method = "CG")$par
  
  alpha <- alpha_fun(Y, eta=eta)
  
  my_weibul_params <- list("alpha" = alpha, "eta" = eta)
  my_weibul_params
  
  # for comparison
  weibul_fit <- EnvStats::eweibull(Y, method = "mle")
  weibul_fit$parameters
  
  return(list("alpha" = alpha, "eta" = eta, 
              "eweibul_alpha" = weibul_fit$parameters[[2]],
              "eweibul_eta" = weibul_fit$parameters[[1]]))
}
g1_stuff$weibul <- my_weibul_fun(Y=gps[[1]])
g2_stuff$weibul <- my_weibul_fun(Y=gps[[2]])
g3_stuff$weibul <- my_weibul_fun(Y=gps[[3]])
g4_stuff$weibul <- my_weibul_fun(Y=gps[[4]])
