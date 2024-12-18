# student: Raul Torres Aragon
# date: 2022-10-13
# assignment: 570 hw 2
# notes:



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# QUESTION 1
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list = ls())
set.seed(570)
b0 <- 3
b1 <- -3
Nsim <- 1000

# create function to do simulations 

run_sims <- function(Nsim = Nsim, b0=3, b1=-3, n=15) {
  
 x <- rnorm(mean=20, sd=2, n=n)

 ModName      <- vector(mode = "character", length = n*3)
 Biases       <- vector(mode = "numeric", length = n*3)
 Var_beta0    <- vector(mode = "numeric", length = n*3)
 Var_beta1    <- vector(mode = "numeric", length = n*3)
 Beta0s       <- vector(mode = "numeric", length = n*3)
 Beta1s       <- vector(mode = "numeric", length = n*3)
 CI_beta0_low <- vector(mode = "numeric", length = n*3) 
 CI_beta0_hig <- vector(mode = "numeric", length = n*3)
 CI_beta1_low <- vector(mode = "numeric", length = n*3)
 CI_beta1_hig <- vector(mode = "numeric", length = n*3)
 
 ynames <- c("Y_norm","Y_unif","Y_skew")
 j = 0
 for(i in 1:Nsim) {
   
   y_norm <- b0 + b1*x + rnorm(mean=0, sd=2, n=n)
   y_unif <- b0 + b1*x + runif(min=-5, max=5, n=n)
   y_skew <- b0 + b1*x + sn::rsn(alpha=5, omega=1, xi=-5/(sqrt(13*pi)), n=n)[1:n]
   
   ys <- list(y_norm, y_unif, y_skew)
   
   k = 0
   for(y in ys){
     k = k+1
     j = j+1
     mymod <- lm(y~x)
     yhat  <- predict(mymod)
     bias  <- mean(yhat - y)
     CIs   <- confint(mymod, param = c(0,1), level = 0.80)
     
     ModName[j]   <- paste(ynames[k], i)
     Biases[j]    <- bias
     Var_beta0[j] <- vcov(mymod)[1,1]
     Var_beta1[j] <- vcov(mymod)[2,2]
     Beta0s[j]    <- mymod$coefficients[[1]]
     Beta1s[j]    <- mymod$coefficients[[2]]
     CI_beta0_low[j] <- CIs[1]
     CI_beta1_low[j] <- CIs[2]
     CI_beta0_hig[j] <- CIs[3]
     CI_beta1_hig[j] <- CIs[4]    
     
   }
 }
 
 results <- dplyr::tibble(model = ModName, 
                          bias = Biases,
                          var_beta0 = Var_beta0,
                          var_beta1 = Var_beta1,
                          b0 = Beta0s,
                          b1 = Beta1s,
                          b0_ci_low = CI_beta0_low,
                          b0_ci_hig = CI_beta0_hig,
                          b1_ci_low = CI_beta1_low,
                          b1_ci_hig = CI_beta1_hig                         
 )
 results$type <- stringr::str_sub(results$model, 3, 6)
 results$b0_inside <- 0
 results$b0_inside[results$b0_ci_low < b0 & b0 < results$b0_ci_hig] <- 1
 results$b1_inside <- 0
 results$b1_inside[results$b1_ci_low < b1 & b1 < results$b1_ci_hig] <- 1
 
 
 # confirm numerically that the bias is zero
 bias_results <- sum(results$bias != 0)

 # compare variance of betas with the sampling distribution of betas_hat
 sampl_var_b0_norm <- results$b0[results$type == "norm"] |> var() 
 sampl_var_b1_norm <- results$b1[results$type == "norm"] |> var() 
 
 sampl_var_b0_unif <- results$b0[results$type == "unif"] |> var() 
 sampl_var_b1_unif <- results$b1[results$type == "unif"] |> var() 
 
 sampl_var_b0_skew <- results$b0[results$type == "skew"] |> var() 
 sampl_var_b1_skew <- results$b1[results$type == "skew"] |> var() 
 
 mean_varb0_norm <- results$var_beta0[results$type == "norm"] |> mean() 
 mean_varb0_unif <- results$var_beta0[results$type == "unif"] |> mean() 
 mean_varb0_skew <- results$var_beta0[results$type == "skew"] |> mean() 
 
 mean_varb1_norm <- results$var_beta1[results$type == "norm"] |> mean() 
 mean_varb1_unif <- results$var_beta1[results$type == "unif"] |> mean() 
 mean_varb1_skew <- results$var_beta1[results$type == "skew"] |> mean() 
 
 # examine the distribution of the resultant estimators (across sims) of betas...
 ci_b0_norm_inside <- mean(results[results$type == "norm", ]$b0_inside)
 ci_b0_unif_inside <- mean(results[results$type == "unif", ]$b0_inside)
 ci_b0_skew_inside <- mean(results[results$type == "skew", ]$b0_inside)
 ci_b1_norm_inside <- mean(results[results$type == "norm", ]$b1_inside)
 ci_b1_unif_inside <- mean(results[results$type == "unif", ]$b1_inside)
 ci_b1_skew_inside <- mean(results[results$type == "skew", ]$b1_inside)
 
 
 o <- list("bias_results"      = bias_results, 
           "ci_b0_norm_inside" = ci_b0_norm_inside,
           "ci_b0_unif_inside" = ci_b0_unif_inside,
           "ci_b0_skew_inside" = ci_b0_skew_inside,
           "ci_b1_norm_inside" = ci_b1_norm_inside,
           "ci_b1_unif_inside" = ci_b1_unif_inside,
           "ci_b1_skew_inside" = ci_b1_skew_inside,
           "sampl_var_b0_norm" = sampl_var_b0_norm,
           "sampl_var_b1_norm" = sampl_var_b1_norm,
           "sampl_var_b0_unif" = sampl_var_b0_unif,
           "sampl_var_b1_unif" = sampl_var_b1_unif,
           "sampl_var_b0_skew" = sampl_var_b0_skew,
           "sampl_var_b1_skew" = sampl_var_b1_skew,
           "mean_varb0_norm"   = mean_varb0_norm,
           "mean_varb0_unif"   = mean_varb0_unif,
           "mean_varb0_skew"   = mean_varb0_skew,
           "mean_varb1_norm"   = mean_varb1_norm,
           "mean_varb1_unif"   = mean_varb1_unif,
           "mean_varb1_skew"   = mean_varb1_skew,           
           "results_df"        = results
           )
 
 return(o)
 
} # end function

o15 <- run_sims(Nsim = Nsim, n = 15)
o40 <- run_sims(Nsim = Nsim, n = 40)


### (b)
par(mfrow = c(1,3))

hist(o15$results_df$b1[o15$results_df$type == "norm"], main = "b1 normal epsilon", 
     xlab = paste0("mean(lm var) = ", 
                   round(o15$mean_varb1_norm,2), "\n sample var =", 
                   round(o15$sampl_var_b1_norm,2)))

hist(o15$results_df$b1[o15$results_df$type == "unif"], main = "b1 normal epsilon", 
     xlab = paste0("mean(lm var) = ", 
                   round(o15$mean_varb1_unif,2), "\n sample var =", 
                   round(o15$sampl_var_b1_unif,2)))

hist(o15$results_df$b1[o15$results_df$type == "skew"], main = "b1 normal epsilon", 
     xlab = paste0("mean(lm var) = ", 
                   round(o15$mean_varb1_skew,2), "\n sample var =", 
                   round(o15$sampl_var_b1_skew,2)))


### (c)
library(plotrix)

sn <- 10

samp <- sample(1:1000, size = sn)

plotrix::plotCI(x=1:sn, rep(-3,sn), 
                abs(o15$results_df$b0_ci_low[o15$results_df$type == "norm"][samp]), 
                abs(o15$results_df$b0_ci_hig[o15$results_df$type == "norm"][samp]), 
                main = "100 randomly selected 10% CI for b1", ylab = "", xlab = "")

plotrix::plotCI(x=1:sn, rep(-3,sn), 
                abs(o15$results_df$b1_ci_low[o15$results_df$type == "norm"][samp]), 
                abs(o15$results_df$b1_ci_hig[o15$results_df$type == "norm"][samp]), 
                main = "100 randomly selected 10% CI for b1", ylab = "", xlab = "")


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# QUESTION 2
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls())

x <- c(6.2, 4.2, 0.5,  8.8, 1.5, 9.2, 8.5, 8.7, 
       6.7, 6.5, 6.3, 6.7,  0.2, 8.7, 7.5)

y <- c(0.8, 3.5, 12.4, 1.1, 8.9, 2.4, 0.1, 0.5, 
       3.5, 8.3, 2.6, 1.5, 16.6, 0.1, 1.4)

loglike <- function(betas){
  sum(betas[1] + betas[2] * x) - sum(y*exp((betas[1] + betas[2] * x)))
}
my_optim <- optim(par = c(1,1), fn = loglike, control = list(fnscale = -1))

b0_hat <- my_optim$par[1]
b1_hat <- my_optim$par[2]

FIM <- matrix(c(length(x), sum(x), sum(x), sum(x^2)), ncol = 2, nrow = 2)
varcovar <- solve(FIM)

b0_ci <- matrix(c(b0_hat - 1.96 * sqrt(varcovar[1,1]), 
                  b0_hat + 1.96 * sqrt(varcovar[1,1])), nrow=1)

b1_ci <- matrix(c(b1_hat - 1.96 * sqrt(varcovar[2,2]), 
                  b1_hat + 1.96 * sqrt(varcovar[2,2])), nrow=1)



# (d) Plot
library(plotly)

b0 <- seq(-5, -1, length.out = 10)
b1 <- seq(0.1, 0.6, length.out = 10)

grid <- expand.grid(b0, b1)
z <- NULL

for(i in 1:length(grid[,1])) {
  z <- cbind(z, loglike(c(grid[i, 1], grid[1, 2])))
}

z_mat <- matrix(data = z, nrow = length(b0), ncol =length(b1))
fig <- plot_ly() %>% add_surface(x = b1, y = b0, z = z_mat, type = "mesh3d")
fig


# (f)
b0_null <- log(15/63.7)
I <- FIM
S <- matrix(c(15 - sum(y * exp(b0_null)), 
              sum(x) - sum(y * x * exp(b0_null))), ncol = 2, nrow = 1)
score <- S %*% solve(I) %*% t(S)
pval_score <- pchisq(score, 1, lower.tail = F)

# wald
I_112 <- I[1,1] - I[1,2]*(I[2,2])^(-1)*I[2,1]
wald_stat <- b1_hat^2 * I_112
pval_wald <- pchisq(wald_stat, 1, lower.tail = FALSE)

# likelihood ratio 
likeli_stat <- 2 * (loglike(c(b0_hat,b1_hat)) - loglike(c(b0_null, 0)))
pval_likeratio <- pchisq(likeli_stat, 1, lower.tail = FALSE)


