#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Author: Raul Torres Aragon
# Date: 2022-09-06
# Assignment: Stat 570 hw1
# Notes:
#
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls())
library(tidyverse)
library(gamlss.data)
data(usair)

# Computation: Question 1
# ~~~~~~~~~~~~~~~~~~~~~~~
data <- usair
data$intercept <- 1
X <- data[, c(8, 2:7)] |> as.matrix()
y <- as.matrix(data[,1])

XX_inv <- solve(t(X)%*%X)
beta <- XX_inv %*% t(X) %*% y

res <- y -(X %*% beta)
rss <- sum(res^2) 
df <- nrow(X) - ncol(X)
sigma <- rss/df
vcovar <- sigma * XX_inv
beta_SE <- diag(vcovar) |> sqrt()

tvals <- (beta-0)/beta_SE
alpha <- 0.05
pvals <- 2*pt(-abs(tvals),df=df)

rse <- sqrt(rss/df)


# Interpretation: Question 1
# ~~~~~~~~~~~~~~~~~~~~~~~~~~
par(mfrow = c(2,3))
plot(X[,2], y, 
     xlab = "avg. annual temperature (F)", 
     ylab = "observed SO2 values")
abline(beta[1], beta[2])

plot(X[,3], y, 
     xlab = "# manufacturers employing > 20 workers", 
     ylab = "observed SO2 values")
abline(beta[1], beta[3])

plot(X[,4], y, 
     xlab = "population size in 1000s", 
     ylab = "observed SO2 values")
abline(beta[1], beta[4])

plot(X[,5], y, 
     xlab = "avg. annual wind speed in miles per hour", 
     ylab = "observed SO2 values")
abline(beta[1], beta[5])

plot(X[,6], y, 
     xlab = "avg. annual rainfall in inches", 
     ylab = "observed SO2 values")
abline(beta[1], beta[6])

plot(X[,7], y, 
     xlab = "avg. number of days of rainfall per year", 
     ylab = "observed SO2 values")
abline(beta[1], beta[7])

mtext("Linear association between SO2 level and metric", side = 3, line = -2, outer = TRUE)

