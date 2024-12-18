#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# student: Raul
# course: STAT 570
# assignment: hw7
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(tidyverse)
library(dplyr)
rm(list=ls())

#~~~~~~~~~~~#
# Problem 1 #
#~~~~~~~~~~~#
n <- 1e3
b0 <- -2.5

mysim <- function(b1, b2, N) {
  betas_1 <- rep(0,N)
  betas_1_star <- rep(0,N)
  inside <- rep(0,N)
  inside_star <- rep(0,N)
  signif_b1 <- rep(0,N)
  signif_b1_star <- rep(0,N)
  cz <- qnorm(0.975, 0, 1)
  
  for(i in 1:N){
    x <- rbinom(n=n, size=1, prob=0.5)
    z <- rnorm(n=n, mean=0, sd=1)
    p <- exp(b0+b1*x+b2*z)/(1 + exp(b0+b1*x+b2*z))
    y <- rbinom(n=n, size=1, prob=p)
    
    fit_tru <- glm(y~x+z, family = binomial(link="logit")) |> summary() %>% .$coefficients
    fit_noz <- glm(y~x  , family = binomial(link="logit")) |> summary() %>% .$coefficients
    
    betas_1[i] <- fit_tru[2,1]
    betas_1_star[i] <- fit_noz[2,1]
    inside[i]      <- as.numeric(fit_tru[2,1]-cz*fit_tru[2,2] < b1 & b1 < fit_tru[2,1]+cz*fit_tru[2,2])
    inside_star[i] <- as.numeric(fit_noz[2,1]-cz*fit_noz[2,2] < b1 & b1 < fit_noz[2,1]+cz*fit_noz[2,2])
    signif_b1[i]      <- as.numeric(fit_tru[2,4] < 0.05)
    signif_b1_star[i] <- as.numeric(fit_noz[2,4] < 0.05)
  }
  return(data.frame("betas_1" = betas_1,
                    "betas_1_star" = betas_1_star,
                    "inside" = inside,
                    "inside_star" = inside_star,
                    "signif_b1" = signif_b1,
                    "signif_b1_star" = signif_b1_star))
}

N <- 1e3
sim_04_06 <- mysim(0.4, 0.6, N=N)
sim_04_12 <- mysim(0.4, 1.2, N=N)
sim_04_25 <- mysim(0.4, 2.5, N=N)
sim_11_06 <- mysim(1.1, 0.6, N=N)
sim_11_12 <- mysim(1.1, 1.2, N=N)
sim_11_25 <- mysim(1.1, 2.5, N=N)

get_SE <- function(x) sqrt(var(x)/length(x))

row1 <-c(sapply(sim_04_06, mean), "SE_b1"=get_SE(sim_04_06$betas_1), "SE_b1_star"=get_SE(sim_04_06$betas_1_star))
row2 <-c(sapply(sim_04_12, mean), "SE_b1"=get_SE(sim_04_12$betas_1), "SE_b1_star"=get_SE(sim_04_12$betas_1_star))
row3 <-c(sapply(sim_04_25, mean), "SE_b1"=get_SE(sim_04_25$betas_1), "SE_b1_star"=get_SE(sim_04_25$betas_1_star))
row4 <-c(sapply(sim_11_06, mean), "SE_b1"=get_SE(sim_11_06$betas_1), "SE_b1_star"=get_SE(sim_11_06$betas_1_star))
row5 <-c(sapply(sim_11_12, mean), "SE_b1"=get_SE(sim_11_12$betas_1), "SE_b1_star"=get_SE(sim_11_12$betas_1_star))
row6 <-c(sapply(sim_11_25, mean), "SE_b1"=get_SE(sim_11_25$betas_1), "SE_b1_star"=get_SE(sim_11_25$betas_1_star))

res_tab <- rbind(row1,row2,row3,row4,row5,row6) |> as.data.frame()

res_tab <- res_tab |> 
           dplyr::select(betas_1, betas_1_star, SE_b1, SE_b1_star, 
                         inside, inside_star, signif_b1, signif_b1_star)

row.names(res_tab) <- c("b1=0.4, b2=0.6","b1=0.4, b2=1.2","b1=0.4, b2=2.5",
                        "b1=1.1, b2=0.6","b1=1.1, b2=1.2","b1=1.1, b2=2.5")
colnames(res_tab) <- c("b1","b1*","SE b1","SE b1*",
                       "95% covrge","95% covrge*",
                       "Power b1", "Power b1*")

# compare with linear model
mysim_lm <- function(b1, b2, N, star = FALSE) {
  betas_1_lm   <- rep(0,N)
  inside_lm    <- rep(0,N)
  sig_b1_lm    <- rep(0,N)
  cz <- qnorm(0.975, 0, 1)
  
  for(i in 1:N){
    x <- rbinom(n=n, size=1, prob=0.5)
    z <- rnorm(n=n, mean=0, sd=1)
    p <- exp(b0+b1*x+b2*z)/(1 + exp(b0+b1*x+b2*z))
    y <- rbinom(n=n, size=1, prob=p)
    if(star == FALSE) {
      fit_lm <-  lm(y~x+z) |> summary() %>% .$coefficients
      if(i<10) print("lm(y~x+z)")
    }
    else {
      fit_lm <-  lm(y~x) |> summary() %>% .$coefficients
    }
    betas_1_lm[i] <- fit_lm[2,1]
    inside_lm[i]  <- as.numeric( fit_lm[2,1]-cz*fit_lm[2,2]  < b1 & b1 <  fit_lm[2,1]+cz*fit_lm[2,2])
    sig_b1_lm[i] <- as.numeric(fit_lm[2,4] < 0.05)
  }
  return(data.frame("betas_1_lm" = betas_1_lm,
                    "inside_lm" = inside_lm,
                    "signif_b1_lm" = sig_b1_lm))
}



get_results_tab <- function(sim1, sim2, sim3, sim4, sim5, sim6, star = FALSE){ 
  
  row1 <-c(sapply(sim1, mean), "SE_b1"=get_SE(sim1[,3]))
  row2 <-c(sapply(sim2, mean), "SE_b1"=get_SE(sim2[,3]))
  row3 <-c(sapply(sim3, mean), "SE_b1"=get_SE(sim3[,3]))
  row4 <-c(sapply(sim4, mean), "SE_b1"=get_SE(sim4[,3]))
  row5 <-c(sapply(sim5, mean), "SE_b1"=get_SE(sim5[,3]))
  row6 <-c(sapply(sim6, mean), "SE_b1"=get_SE(sim6[,3]))
  
  rt <- rbind(row1,row2,row3,row4,row5,row6) |> as.data.frame()
  
  rt <- rt[,c(1,4,2,3)] #|> dplyr::select(betas_1_lm, SE_b1_lm, inside_lm, signif_b1_lm)
  
  row.names(rt) <- c("b1=0.4, b2=0.6","b1=0.4, b2=1.2","b1=0.4, b2=2.5",
                     "b1=1.1, b2=0.6","b1=1.1, b2=1.2","b1=1.1, b2=2.5")
  
  if(star==TRUE) {
    colnames(rt) <- c("b1*","SE b1*","95% covrge*","Power b1*")
  } else {
    colnames(rt) <- c("b1","SE b1","95% covrge","Power b1")
  }
  rt
}

N <- 1e3
sim_lm_04_06 <- mysim_lm(0.4, 0.6, N=N, star = FALSE)
sim_lm_04_12 <- mysim_lm(0.4, 1.2, N=N, star = FALSE)
sim_lm_04_25 <- mysim_lm(0.4, 2.5, N=N, star = FALSE)
sim_lm_11_06 <- mysim_lm(1.1, 0.6, N=N, star = FALSE)
sim_lm_11_12 <- mysim_lm(1.1, 1.2, N=N, star = FALSE)
sim_lm_11_25 <- mysim_lm(1.1, 2.5, N=N, star = FALSE)
res_tab_lm_z <- get_results_tab(sim_lm_04_06, sim_lm_04_12, sim_lm_04_25, 
                                sim_lm_11_06, sim_lm_11_12, sim_lm_11_25, star = FALSE)

N <- 1e3
sim_lm_04_06 <- mysim_lm(0.4, 0.6, N=N, star = TRUE)
sim_lm_04_12 <- mysim_lm(0.4, 1.2, N=N, star = TRUE)
sim_lm_04_25 <- mysim_lm(0.4, 2.5, N=N, star = TRUE)
sim_lm_11_06 <- mysim_lm(1.1, 0.6, N=N, star = TRUE)
sim_lm_11_12 <- mysim_lm(1.1, 1.2, N=N, star = TRUE)
sim_lm_11_25 <- mysim_lm(1.1, 2.5, N=N, star = TRUE)
res_tab_lm <- get_results_tab(sim_lm_04_06, sim_lm_04_12, sim_lm_04_25, 
                              sim_lm_11_06, sim_lm_11_12, sim_lm_11_25, star = FALSE)

res_tab_lm
res_tab_lm_z


#~~~~~~~~~~~#
# Problem 2 #
#~~~~~~~~~~~#

library(mlbench)
data(PimaIndiansDiabetes2)


# (a)
# remove records with missing values
mydata <- PimaIndiansDiabetes2[, c("diabetes","pregnant","glucose","mass","pedigree","age")]
mydata <- mydata[complete.cases(mydata), ]
mydata$diabetes <- as.numeric(mydata$diabetes == "pos")
glimpse(mydata)

# (b)
# fit the three models
myformula <- as.formula("diabetes ~ pregnant + glucose + mass + pedigree + age")
fit_logit   <- glm(myformula, data=mydata, binomial(link="logit"))
fit_probit  <- glm(myformula, data=mydata, binomial(link="probit"))
fit_cloglog <- glm(myformula, data=mydata, binomial(link="cloglog"))

get_model_stuff <- function(model) {
  l <- model$rank
  even_i <- seq(from=0, by=2, len=l)+2
  odd_i <- seq(from=1, by=2, len=l)
  v <- rep(0, 2*l)
  v[odd_i] <- summary(model)$coefficients[,1] |> as.vector() |> formatC(format="f",digits=2)
  v[even_i] <- summary(model)$coefficients[,2] |> as.vector() |> formatC(format="f",digits=4)
  return(v)
}

tab_b <- data.frame("logit"  = get_model_stuff(fit_logit),
                    "probit" = get_model_stuff(fit_probit),
                    "cloglog"= get_model_stuff(fit_cloglog))
mynames <- rep("",12)
mynames[seq(from=0, by=2, len=6)+2] <- paste("SE", names(coef(fit_logit)))
mynames[seq(from=1, by=2, len=6)] <-  names(coef(fit_logit))
row.names(tab_b) <- mynames


## (c)
tab_c <- tab_b[seq(from=1, by=2, len=6), 1] |> 
               cbind(as.numeric(tab_b[seq(from=1, by=2, len=6),1]) |> 
                     exp() |> 
                     formatC(format="f", digits=2))
row.names(tab_c) <- names(coef(fit_logit))
colnames(tab_c) <- c("logit  b", "logit exp(b)")


## (d)  
xseq <- seq(-5, 5, length.out=200)
X_means_df <- data.frame("pregnant" = median(mydata$pregnant),
                         "glucose"  = median(mydata$glucose),
                         "mass"     = median(mydata$mass),
                         "pedigree" = xseq,
                         "age"      = median(mydata$age))

fitted_log <- predict(fit_logit,   newdata = X_means_df, se.fit=TRUE, type="response") 
fitted_cll <- predict(fit_cloglog, newdata = X_means_df, se.fit=TRUE, type="response")
fitted_pro <- predict(fit_probit,  newdata = X_means_df, se.fit=TRUE, type="response")


# Three in same plot
plot(diabetes~pedigree, data = mydata,
     ylim=c(0,1),
     ylab="Probability of Diabetes",
     xlab="Diabetes Pedigree",
     xlim=c(-1.5,3.5),
     main = "Association between pedigree and prevalence of diabetes\nfor three link functions in a binomial model")
lines(xseq, fitted_log$fit, col="darkgreen")
lines(xseq, fitted_pro$fit, col="darkorange")
lines(xseq, fitted_cll$fit, col="darkblue")
legend("topleft",
       fill=c("darkgreen","darkblue","darkorange"),
       col=c("darkgreen","darkblue","darkorange"),
       bty="n",
       legend=c("Logistic","Comp log-log","Probit"))


#Three plots
z <- qnorm(0.96) # for 92% confidence intervals

par(mfrow=c(1,3))

# logistic
plot(diabetes~pedigree, data = mydata,
     ylim=c(0,1),
     ylab="Probability of Diabetes",
     xlab="Diabetes Pedigree",
     xlim=c(-1,3), main="link=logit")
lines(xseq, fitted_log$fit, col="darkgreen")
lines(xseq, fitted_log$fit + z*fitted_log$se.fit, col="darkgreen", lty=2)
lines(xseq, fitted_log$fit - z*fitted_log$se.fit, col="darkgreen", lty=2)

# cloglog
plot(diabetes~pedigree, data = mydata,
     ylim=c(0,1),
     ylab="Probability of Diabetes",
     xlab="Diabetes Pedigree",
     xlim=c(-1,3), main = "link=cloglog")
lines(xseq, fitted_cll$fit, col="darkblue")
lines(xseq, fitted_cll$fit + z*fitted_cll$se.fit, col="darkblue", lty=2)
lines(xseq, fitted_cll$fit - z*fitted_cll$se.fit, col="darkblue", lty=2)

# probit
plot(diabetes~pedigree, data = mydata,
     ylim=c(0,1),
     ylab="Probability of Diabetes",
     xlab="Diabetes Pedigree",
     xlim=c(-1,3), main = "link=probit")
lines(xseq, fitted_pro$fit, col="darkorange")
lines(xseq, fitted_pro$fit + z*fitted_pro$se.fit, col="darkorange", lty=2)
lines(xseq, fitted_pro$fit - z*fitted_pro$se.fit, col="darkorange", lty=2)
mtext("Predicted probabilities and 92% confidence interval", side = 3, line = -1.25, outer = TRUE)


# (d) DEPRECATED
### xseq <- seq(0, 2.5, length.out=100)

#fitted_log <- predict(fit_logit,   newdata = X_means_df, se.fit=TRUE, type="link") 
#fitted_cll <- predict(fit_cloglog, newdata = X_means_df, se.fit=TRUE, type="link")
#fitted_pro <- predict(fit_probit,  newdata = X_means_df, se.fit=TRUE, type="link")

# logit
#lines(xseq, exp(fitted_log$fit)/(1+exp(fitted_log$fit)), col="darkgreen")
#lines(xseq, exp(fitted_log$fit+z*fitted_log$se.fit)/(1+exp(fitted_log$fit+z*fitted_log$se.fit)), col="darkgreen", lty=2)
#lines(xseq, exp(fitted_log$fit-z*fitted_log$se.fit)/(1+exp(fitted_log$fit-z*fitted_log$se.fit)), col="darkgreen", lty=2)

# cloglog
#lines(xseq, 1-exp(-exp(fitted_cll$fit)), col = "darkblue")
#lines(xseq, 1-exp(-exp(fitted_cll$fit +z*fitted_cll$se.fit)), col = "darkblue", lty=2)
#lines(xseq, 1-exp(-exp(fitted_cll$fit -z*fitted_cll$se.fit)), col = "darkblue", lty=2)



save.image(file = "/Users/raulta/Library/CloudStorage/GoogleDrive-rauldavidt@gmail.com/My Drive/Stat 570/hw/hw7/hw7_objects.Rdata")







