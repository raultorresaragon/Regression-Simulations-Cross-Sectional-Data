# Student: Raul
# Course: STAT 570
# Assignment: hw4
# Date: 2022-10-25
#------------------------------------------# 

# QUESTION 1 #
# -----------#
rm(list = ls())
set.seed(570)
ns <- 10 # c(5,10,20,30,40,50)
rhos <- 0.5 # c(0.1, 0.5, 0.9)
b0 <- 4
b1 <- 1.75
sigma2 <- 1

# proof of concept
t <- letters[1:5]
expand.grid(t(t), t) -> d
tjtk <- stringr::str_c(d$Var1, d$Var2)
V <- matrix(tjtk, nrow=5,ncol=5)
for(row in 1:5) {
  for(col in 1:5){
    if(row==col) { V[row,col] <- paste0(V[row,col],"^",0) }
    if(row!=col) { V[row,col] <- paste0(V[row,col],"^",max(row,col)-1)}
  }
}
V
rm(list = c("t","d","V","tjtk"))


# run model
run_model <- function(b0=4, b1=1.75, sigma2=1, n, rho, use_sigma2_hat=FALSE) {
  
  t <- seq(-2, 2, length.out = n)
  t_star <- t - mean(t)
  X <- cbind(1, t_star)
  Mu<- mean(b0 + b1*(t_star)) #<-mean
  
  results <- list("m1" = NA, "m2" = NA, "m3" = NA)
  for(type in 1:3) {
    if(type==1) { e <- rnorm(n=n, mean=0, sqrt(Mu*sigma2)) }
    if(type==2) { e <- rnorm(n=n, mean=0, sqrt(Mu^2*sigma2)) }
    if(type==3) { 
      expand.grid(t(t), t) -> d; rho^abs(d$Var1 - d$Var2) -> deltas
      V <- matrix(deltas, nrow = n, ncol = n)
      for(row in 1:n) {
        for(col in 1:n){
          if(row==col) {
            V[row,col] <- V[row,col]^0 
          } else{ 
            V[row,col] <- V[row,col]^(max(row,col)-1) 
          }
        }
      }
      e <- mgcv::rmvn(n=1, mu = rep(0,n), V=V)
    }
    
    Y <- b0 + b1*(t_star) + e
     
    if(type==1) { vcovar <- Mu*sigma2 * solve(t(X)%*%X) }
    if(type==2) { vcovar <- Mu^2*sigma2 * solve(t(X)%*%X) }
    if(type==3) {
      vcovar <- sigma2 * solve(t(X)%*%X)
      vcovar[1,2] <- vcovar[2,1] <- sigma2*rho^sum(abs(1-t)) #<-sum???
    }
    betas <- solve(t(X)%*%X)%*%t(X)%*%Y
    
    # rss/dfreed is an unbiased estimator for var(e)
    # do we compute C.I. with rss/dfreed or with the known var(e)?
    if(use_sigma2_hat == TRUE) {
      #print("using rss/dfreed to estimate the error variance")
      sigma2_hat <- sum((Y -(X %*% betas))^2) / (nrow(X) - ncol(X))
      vcov_hat <- sigma2_hat * solve(t(X)%*%X)
      vcovar <- vcov_hat
    }
    
    c_i_b0 <- c(betas[1] - 1.96 * sqrt(vcovar[1,1]),
                betas[1] + 1.96 * sqrt(vcovar[1,1]))
    c_i_b1 <- c(betas[2] - 1.96 * sqrt(vcovar[2,2]),
                betas[2] + 1.96 * sqrt(vcovar[2,2]))    
    results[[type]] <- list("c_i_b0" = c_i_b0, 
                            "c_i_b1" = c_i_b1,
                            "betas" = betas,
                            "vcovar" = vcovar)
  }
  return(results)
}


simulate <- function(N, n, rho) {

   coverage_m1_b0 <- vector(mode = "numeric", length = N)
   coverage_m1_b1 <- vector(mode = "numeric", length = N)
   coverage_m2_b0 <- vector(mode = "numeric", length = N)
   coverage_m2_b1 <- vector(mode = "numeric", length = N)
   coverage_m3_b0 <- vector(mode = "numeric", length = N)
   coverage_m3_b1 <- vector(mode = "numeric", length = N)
   
   for(s in 1:N) {
     
     r <- run_model(n=n, rho=rho, use_sigma2_hat = TRUE)
     
     if(r$m1$c_i_b0[1] < b0 & b0 < r$m1$c_i_b0[2]) {
       coverage_m1_b0[s] <- 1
     }  else { coverage_m1_b0[s] <- 0}
     if(r$m1$c_i_b1[1] < b1 & b1 < r$m1$c_i_b1[2]) {
       coverage_m1_b1[s] <- 1
     }  else { coverage_m1_b1[s] <- 0}
     
     
     if(r$m2$c_i_b0[1] < b0 & b0 < r$m2$c_i_b0[2]) {
       coverage_m2_b0[s] <- 1
     }  else { coverage_m2_b0[s] }
     if(r$m2$c_i_b1[1] < b1 & b1 < r$m2$c_i_b1[2]) {
       coverage_m2_b1[s] <- 1
     }  else { coverage_m2_b1[s] <- 0}
     
     
     if(r$m3$c_i_b0[1] < b0 & b0 < r$m3$c_i_b0[2]) {
       coverage_m3_b0[s] <- 1
     }  else { coverage_m3_b0[s] <- 0}
     if(r$m3$c_i_b1[1] < b1 & b1 < r$m3$c_i_b1[2]) {
       coverage_m3_b1[s] <- 1
     }  else { coverage_m3_b1[s] <- 0}
   }
  
  return(list(
    "coverage_m1_b0" = coverage_m1_b0 |> mean(),
    "coverage_m1_b1" = coverage_m1_b1 |> mean(),
  
    "coverage_m2_b0" = coverage_m2_b0 |> mean(),
    "coverage_m2_b1" = coverage_m2_b1 |> mean(),
  
    "coverage_m3_b0" = coverage_m3_b0 |> mean(),
    "coverage_m3_b1" = coverage_m3_b1 |> mean()  
  ))
}
mytab <- tibble::tibble("rho" = numeric(),
                        "n"   = numeric(),
                        "beta" = character(),
                        "coverage_m1" = numeric(),
                        "coverage_m2" = numeric(),
                        "coverage_m3" = numeric())
for(r in rhos) {
  for(n in ns) {
    s <- simulate(N=1e2, n=n, rho=r)
    mytab <- mytab |>
      dplyr::add_row("rho"=r, "n"=n, "beta"="b0", 
            "coverage_m1"=s$coverage_m1_b0,
            "coverage_m2"=s$coverage_m2_b0,
            "coverage_m3"=s$coverage_m3_b0) |>
      dplyr::add_row("n"=n, "rho"=r, "beta"="b1", 
            "coverage_m1"=s$coverage_m1_b1,
            "coverage_m2"=s$coverage_m2_b1,
            "coverage_m3"=s$coverage_m3_b1)
  }
}














