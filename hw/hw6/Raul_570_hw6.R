#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# student: Raul
# course: STAT 570
# assignment: hw6
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list=ls())
library(tidyverse)
library(dplyr)


nh <-c(375,	375,	752,	208,	151,	116,	736,	192,	315,	1252,	675,	700,	440,	771,	
       688,	426,	410,	979,	377,	503)
ho <-c(396,	568,	1212,	171,	554,	1104,	257,	435,	295,	397,	288,	1004,	431,	795,	
       1621,	1378,	902,	958,	1283,	2415)


#~~~~~~~~~~~#
# Problem 1 # 
#~~~~~~~~~~~#

## (a)
means <- c(mean(nh), mean(ho))
boxplot(data.frame("Non_Hodgkins" = nh, "Hodgkins" = ho),
        main = "Cell counts\n between two diseases",
        horizontal = TRUE)
text(c(mean(nh), min(nh), max(nh), min(ho), max(ho), mean(ho)),  
     c(rep(0.8,3), rep(1.8,3)),
     c(mean(nh), min(nh), max(nh), min(ho), max(ho), mean(ho)), 
     pos = 1, cex = 0.5)  
points(means, c(1,2), col="red",pch=1, cex=2)


## (b)  
get_bootstrap_var <- function(x, N) {
    B <- rep(0, N) 
    for(i in 1:N) {
        B[i] <- sample(x, size = length(x), replace = TRUE) |> var()
    }
    var <- mean(B)
    var
}
N <- 1e3
get_CI <- function(x1, x2, level = .90, scale = "orig") {
  if(scale == "log") {
    x1 <- log(x1); x2 <- log(x2)
  } else if(scale == "sqrt") {
    x1 <- sqrt(x1); x2 <- sqrt(x2)
  } else {
    TRUE
  }
  var1 <- get_bootstrap_var(x1, N=1e3)
  var2 <- get_bootstrap_var(x2, N=1e3)
  Phi <- abs(qnorm((1-level)/2, mean = 0, sd = 1))
  ci_lo <- mean(x1)-mean(x2) - Phi*sqrt(var1/length(x1) + var2/length(x2))
  ci_hi <- mean(x1)-mean(x2) + Phi*sqrt(var1/length(x1) + var2/length(x2))
  return(c(ci_lo, ci_hi))
}

CI_orig <- get_CI(ho, nh, level = 0.90, scale = "orig")
CI_log  <- get_CI(ho, nh, level = 0.90, scale = "log")
CI_sqrt <- get_CI(ho, nh, level = 0.90, scale = "sqrt")


## (c)  
data <- data.frame("count" = c(nh, ho), 
                   "hodgkins" = c(rep(0,length(nh)), rep(1,length(ho))))

fit_pois <- glm(count~hodgkins, family=poisson(link="log"), data=data)
fit_gamm <- glm(count~hodgkins, family=Gamma(link="inverse"), data=data) 
fit_invg <- glm(count~hodgkins, family=inverse.gaussian(link="1/mu^2"), data=data)

coef_pois <- summary(fit_pois)$coefficients |> as.data.frame()
coef_gamm <- summary(fit_gamm)$coefficients |> as.data.frame()
coef_invg <- summary(fit_invg)$coefficients |> as.data.frame()

res_tab <- data.frame(
  "Intercept"    = c(coef_pois$Estimate[1], coef_gamm$Estimate[1], coef_invg$Estimate[1]),
  "SE_Intercept" = c(coef_pois$`Std. Error`[1], coef_gamm$`Std. Error`[1], coef_invg$`Std. Error`[1]),
  "Beta1"        = c(coef_pois$Estimate[2], coef_gamm$Estimate[2], coef_invg$Estimate[2]), 
  "SE_Beta1"     = c(coef_pois$`Std. Error`[2], coef_gamm$`Std. Error`[2], coef_invg$`Std. Error`[2])
  )
row.names(res_tab) <- c("Poisson","Gamma","InverseGaussian")

undo_gamma <- function(x) { -1/x }
undo_invg  <- function(x) { 1/sqrt(x) }

coef_pois_undone <- fit_pois$coefficients[,1] |> exp() |> as.data.frame()
coef_gamm_undone <- fit_gamm$coefficients[,1] |> undo_gamma() |> as.data.frame()
coef_invg_undone <- fit_invg$coefficients[,1] |> abs() |> undo_invg() |> as.data.frame()

res_tab_undone <- data.frame(
  "Intercept"    = c(coef_pois_undone[1,1], coef_gamm_undone[1,1], coef_invg_undone[1,1]),
  "Beta1"        = c(coef_pois_undone[2,1], coef_gamm_undone[2,1], coef_invg_undone[2,1]) 
)
row.names(res_tab_undone) <- c("Poisson","Gamma","InverseGaussian")



## (d)  
I_pois <- solve(vcov(fit_pois))
I_gamm <- solve(vcov(fit_gamm))
I_invg <- solve(vcov(fit_invg))

get_CI_model <- function(theta, I_jj) {
  r <- c(theta - qnorm(0.95)*sqrt(1/I_jj), theta + qnorm(0.95)*sqrt(1/I_jj))
  r
}

CI_b0_pois <- get_CI_model(coef_pois[1,1], I_pois[1,1])
CI_b1_pois <- get_CI_model(coef_pois[2,1], I_pois[2,2])

CI_b0_gamm <- get_CI_model(coef_gamm[1,1], I_gamm[1,1])
CI_b1_gamm <- get_CI_model(coef_gamm[2,1], I_gamm[2,2])

CI_b0_invg <- get_CI_model(coef_invg[1,1], I_pois[1,1])
CI_b1_invg <- get_CI_model(coef_invg[2,1], I_pois[2,2])


CI_tab <- data.frame("intercept" = c(paste0(round(CI_b0_pois, 3), collapse = " to "), 
                                     paste0(round(CI_b0_gamm, 3), collapse = " to "),
                                     paste0(round(CI_b0_invg, 3), collapse = " to ")),
                     "b1" = c(paste0(round(CI_b1_pois, 3), collapse = " to "), 
                              paste0(round(CI_b1_gamm, 3), collapse = " to "),
                              paste0(round(CI_b1_invg, 3), collapse = " to "))
           )

row.names(CI_tab) <- c("Poisson","Gamma","InverseGaussian")




#~~~~~~~~~~~#
# Problem 2 # 
#~~~~~~~~~~~#

x <- c(seq(from = 2, to = 10, by = 2), 24, 28, 32)

y <- c(1.63, 1.01, 0.73, 0.55, 0.41, 0.01, 0.06, 0.02)

## a
## b

my_likelihood <- function(theta) {
  b0 <- log(30/theta[1])
  b1 <- theta[2]
  n <- length(x)
  sigma2 <- theta[3]
  -n/2 * log(2*pi*sigma2) -sum(log(y)) -(1/(2*sigma2))*sum( (log(y) -b0+b1*x)^2 ) 
}
params <- optim(par = c(1,1,1), fn = my_likelihood, control = list(fnscale = -1))$par

V <-params[1]; Ke <- params[2]; sigma2 <- params[3]


# TO DO
I <- matrix(1, ncol=3, nrow=3)

q975 <- qnorm(.975, 0, 1)
ci_V <- c(V-q975*sqrt(I[1,1]), V+q975*sqrt(I[1,1])) |> round(2)
ci_Ke <- c(Ke-q975*sqrt(I[2,2]), Ke+q975*sqrt(I[2,2]))  |> round(2)
ci_sigma2 <- c(sigma2-q975*sqrt(I[3,3]), sigma2+q975*sqrt(I[3,3]))  |> round(2)

confidence_intervals <- c(paste0(ci_V, collapse = " to "),
                          paste0(ci_Ke,collapse = " to "),
                          paste0(ci_sigma2, collapse = " to "))

Q2b_tab <- data.frame(round(params,2), confidence_intervals)
row.names(Q2b_tab) <- c("V", "-Ke", "sigma2")
colnames(Q2b_tab) <- c("parameter", "95% conf interval")


## (c)  
#We now plot the data, along with the fitted curve.  
b0 <- log(30/V)
b1 <- Ke
x_grid <- seq(from = 2, to = 33, by = 1)
plot(y~x, main = "Cadralazine concentration at 8 given times", xaxt='n')
lines(exp(b0-b1*x_grid), lty=1, lwd=3)
axis(1,at=x,labels=x)

## (d)  
#Using residuals, we examine the appropriateness of the assumptions of the above model. 
#And answer the question: Does the model seem reasonable for these data?  
yhat <- exp(b0-b1*x)
r <- y - exp(b0-b1*x)
std_r <- (y - yhat) / sqrt(sigma2)
std_re <- (y - yhat) / sd(y - yhat)
#hist(r, breaks = 13)


## (e)  
#The clearance Cl = V * K_e and elimination half-life x_{1/2} = log(2/K_e) are 
#parameters of interest in this experiment.  
#We find the MLEs of these parameters along with asymptotic 95% confidence intervals. 

Cl_mle <- V*Ke
xhalf_mle <- log(2)/Ke 
n <- length(x)

dKe <- -Ke/sigma2^2 * sum(x^2) + 1/sigma2 * sum(log(y)*x) + 1/sigma2 * log(30/V)*sum(x)
dV <- n*log(30)*1/(sigma2*V^2) + n*log(V)/(sigma2*V^2) + sum(log(y))/(sigma2*V^2) + Ke*sum(x)/(sigma2*V^2)
I_11 <- 1/(sigma2*V^2) * sum(x)
I_22 <- 1/(sigma2*V^2) * sum(x)
I_21 <- -sum(x^2)/(sigma2*V^2)

grad_loglike <- matrix(c(-dKe,dV), ncol=2, nrow=1)
FIM <- matrix(c(I_11, I_21, I_21, I_22), nrow=2, ncol=2)

var_Cl <- grad_loglike %*% solve(FIM) %*% t(grad_loglike)
var_xhalf <- -dKe^2 * 1/I_22

CI_cl <- paste0(c(Cl_mle - qnorm(0.975)*sqrt(var_Cl), 
                  Cl_mle + qnorm(0.975)*sqrt(var_Cl)) |> round(3),
                collapse = " to ")

CI_xhalf <- paste0(c(xhalf_mle - qnorm(0.975)*sqrt(var_xhalf), 
                     xhalf_mle + qnorm(0.975)*sqrt(var_xhalf)) |> round(3),
                     collapse = " to ")









save.image(file = "/Users/raulta/Library/CloudStorage/GoogleDrive-rauldavidt@gmail.com/My Drive/Stat 570/hw/hw6/hw6_objects.Rdata")









