n <- 100
X <- rnorm(mean=0, sd=1, n=n)
theta <- rgamma(shape=b, scale=b, n=n)
mu <- exp(b0 + b1*X)
Y <- rpois(lambda = theta*mu, n=n) 

by_hand <- function(X, Y, b) {
  
  # do we optimize loglikelihood 
  ##############################
  
  n <- length(Y)
  
  my_loglike <- function(betas) {
    sum(Y*log(exp(betas[1] + betas[2]*X))) 
    - n*exp(betas[1] + betas[2]*X) 
    -sum(log(factorial(Y)))
  }  
  
  my_optim <- optim(par = c(1,1),
                    fn = my_loglike,
                    control = list(fnscale = -1))
  
  yhat <- my_optim$par[1] + my_optim$par[2]*X 
  alpha <- (n-1)^-1 * sum((Y-yhat)^2/yhat)
  s2 <- (n-1)^-1 * sum(Y-yhat)^2
  
  Xmat <- matrix(c(rep(1,length(X)),X), nrow = length(X), ncol = 2)
  varcovar <- s2 * solve(t(Xmat)%*%Xmat)
  
  var_b0 <- varcovar[1,1] 
  var_b1 <- varcovar[2,2]
  
  CI_b0 <- c(my_optim$par[1] - 1.96*var_b0, my_optim$par[1] + 1.96*var_b0)
  CI_b1 <- c(my_optim$par[2] - 1.96*var_b0, my_optim$par[2] + 1.96*var_b0)
  
  b0_inside <- (CI_b0[1] < b0 & b0 < CI_b0[2]) |> as.numeric()
  b1_inside <- (CI_b1[1] < b0 & b0 < CI_b1[2]) |> as.numeric() 
  
  return(list("b0_inside" = b0_inside, 
              "b1_inside" = b1_inside))
}