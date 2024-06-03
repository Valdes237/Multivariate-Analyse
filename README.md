# Multivariate-Analyse
Mietspiegel Analyse

rm(list=ls())
if (!require("VIM")) install.packages("VIM")
if (!require("mice")) install.packages("mice")
if (!require("dplyr")) install.packages("dplyr")
if(!require("normtest")) install.packages("normtest")
library(mice)
library(VIM)
library(mvtnorm)
library(dplyr)
library(lattice)
library(normtest)

set.seed(123)

  
## Data

mietspiegel # Data set

mietspiegel <- mietspiegel %>% dplyr::mutate(miete = as.numeric(miete),
                                             mieteqm = as.numeric(mieteqm),
                                             bjahr = as.numeric(bjahr),
                                             lage = as.factor(lage),
                                             bad = as.factor(bad),
                                             kueche = as.factor(kueche),
                                             bezv = as.factor(bezv),
                                             zh = as.factor(zh))


# Test regression to get an overview of the data
model_1_mietspiegel <- lm(miete~., data = mietspiegel)
summary(model_1_mietspiegel)

# Test of normal distribution
# Shapiro-wilk-test
attach(mietspiegel)
num.Mat <- apply(mietspiegel,2,as.numeric)
Shap.test <- apply(num.Mat,2,shapiro.test) # the P-value of each variable is lower than the significant levels
detach(mietspiegel)

# remove: mieteqm,bezv
dat_2 <- mietspiegel[,-7]

# recode: bjahr
dat_2$bjahr <- 2019-dat_2$bjahr

# recode: lage as binary variable 0/1
dat_2 <- dat_2 %>% dplyr::mutate(lage =  ifelse(as.numeric(dat_2$lage)>1, 1, 0))

head(dat_2)


model_2_mietspiegel <- lm(miete~., data = dat_2)
summary(model_2_mietspiegel) # all variable are significant

sapply(dat_2, class)

dim(dat_2)





### MC study ###

R <- 1000             # number of simulation cycles
M <- 10             # number of multiple imputations
n <- 500             # number of observations
runs <- 1           # number of iterations


# True Value


true_val_miete <- mean(dat_2$miete)
true_val_flaeche <- mean(dat_2$flaeche)



# Create functions for bias, MSE and coverage

bias <- function(true_val,value) {
  bias <- (true_val-value)/value
  return(bias)
}

MSE <- function(true_val, value) {
  mse <- mean((true_val-mean(value))^2)
  return(mse)
}

coverage <- function(value, CI.lower, CI.upper){
  ifelse(CI.lower < value & CI.upper > value,1,0)
}

# Function for analysis

MI.analysis <- function(theta,within,m) {
  theta_est <- mean(theta)
  # Within_variance
  within_variance <- mean(within)
  # Between_variance
  between_variance <- var(theta)
  # Total variance
  total_variance <- within_variance + (1+1/m)*between_variance 
  # Degrees of Freedom
  dfreedom <- (m-1)*(1+within_variance/((1+1/m)*between_variance))^2
  # Confidence interval
  CI.lower <- theta_est-qt(0.975, dfreedom)*sqrt(total_variance)
  CI.upper <- theta_est+qt(0.975, dfreedom)*sqrt(total_variance)
  return(cbind(theta_est,CI.lower,CI.upper))
}




##### Missing Completely At Random #####

BD_MCAR <- MImean_MCAR <- MInorm_MCAR <- MIpmm_MCAR <- MImidas_MCAR <- array(dim=c(R, 2, 2))



for (r in 1:R) {
  
  ## Sample
  
  sample_dat <- sample(nrow(dat_2), n)
  dat <- dat_2[sample_dat,]
  
  ## Generating Missing Values
  
  #MCAR
  
  p.mis <- 0.30 
  dat.mcar <- dat
  mis.mcar_1 <- sample(1:n,p.mis*n,replace=FALSE)
  mis.mcar_2 <- sample(1:n,p.mis*n,replace=FALSE)
  dat.mcar[mis.mcar_1,1] <- NA
  dat.mcar[mis.mcar_2,2] <- NA
  
  # marginplot(dat.mcar[, c("miete", "flaeche")], col = mdc(1:2), cex = 1.2, cex.lab = 1.2, cex.numbers = 1.3, pch = 19)

  ### Before deletion
  
  # mean:
  mean_est_miete <- mean(dat$miete)
  mean_est_flaeche <- mean(dat$flaeche)
  
  ci <- t.test(dat$miete)$conf.int
  ci_low_1 <- ci[1]
  ci_upp_1 <- ci[2]
  
  ci <- t.test(dat$flaeche)$conf.int
  ci_low_2 <- ci[1]
  ci_upp_2 <- ci[2]
  

  BD_MCAR[r, , 1] <- c(bias(true_val_miete,mean_est_miete),coverage(true_val_miete, ci_low_1, ci_upp_1))
  BD_MCAR[r, , 2] <- c(bias(true_val_flaeche,mean_est_flaeche),coverage(true_val_flaeche, ci_low_2, ci_upp_2))
  
  # initialization:
  ini <- mice(dat.mcar,m=1,maxit=0)
  pred <- quickpred(dat.mcar)
  
  
  ### Mean imputation ###
  meth <- ini$method
  meth[c("miete","flaeche")] <- "mean"
  
  imp_mean <- mice(dat.mcar,m=M,maxit=runs, method = meth, predictorMatrix = pred, print=FALSE)
  dat.mcar_mean <- complete(imp_mean,  action="long", include=FALSE)
  ### Mean nicht sinvoll bei "flaeche", da 65.94857 = unplausibler Wert
  
  ## Analysis ##
  
  
  # Mean for each imputed data set:
  theta <- aggregate(miete ~.imp, data=dat.mcar_mean, mean)$miete
  # Variance of mean for each imputed data set:
  within_variance <- aggregate(miete ~.imp, data = dat.mcar_mean, var)$miete/n
  # MI-Analysis
  MI_mean_1 <- MI.analysis(theta,within_variance,M)
  
  # Mean for each imputed data set:
  theta <- aggregate(flaeche ~.imp, data=dat.mcar_mean, mean)$flaeche
  # Variance of mean for each imputed data set:
  within_variance <- aggregate(flaeche ~.imp, data = dat.mcar_mean, var)$flaeche/n
  # MI-Analysis
  MI_mean_2 <- MI.analysis(theta,within_variance,M)
  
  # Save results in matrix
  MImean_MCAR[r, , 1] <-  c(bias(true_val_miete,MI_mean_1[1]), coverage(true_val_miete, MI_mean_1[2], MI_mean_1[3]))
  MImean_MCAR[r, , 2] <-  c(bias(true_val_flaeche,MI_mean_2[1]), coverage(true_val_flaeche, MI_mean_2[2], MI_mean_2[3]))
  
  
  
  ### Norm imputation ###
  
  meth <- ini$method
  meth[c("miete","flaeche")] <- "norm"
  pred <- quickpred(dat.mcar)
  
  
  imp_norm <- mice(dat.mcar, m=M, maxit=runs, method=meth, predictorMatrix = pred, print=F)
  dat.mcar_norm <- complete(imp_norm, action="long", include=FALSE) # Stack completed data sets
  
  ## Analysis ##
  
  # Mean for each imputed data set:
  theta <- aggregate(miete ~.imp, data = dat.mcar_norm, mean)$miete
  # Variance of mean for each imputed data set:
  within_variance <- aggregate(miete ~.imp, data = dat.mcar_norm, var)$miete/n
  # MI-Analysis
  MI_mean_1 <- MI.analysis(theta,within_variance,M)
  
  # Mean for each imputed data set:
  theta <- aggregate(flaeche ~.imp, data = dat.mcar_norm, mean)$flaeche
  # Variance of mean for each imputed data set:
  within_variance <- aggregate(flaeche ~.imp, data = dat.mcar_norm, var)$flaeche/n
  # MI-Analysis
  MI_mean_2 <- MI.analysis(theta,within_variance,M)
  
  # Save results in matrix
  MInorm_MCAR[r, , 1] <-  c(bias(true_val_miete,MI_mean_1[1]), coverage(true_val_miete, MI_mean_1[2], MI_mean_1[3]))
  MInorm_MCAR[r, , 2] <-  c(bias(true_val_flaeche,MI_mean_2[1]), coverage(true_val_flaeche, MI_mean_2[2], MI_mean_2[3]))
  
  
  
  ### PMM ###
  
  meth <- ini$method
  meth[c("miete","flaeche")] <- "pmm"
  
  
  imp_pmm <- mice(dat.mcar, m=M, maxit=runs, predictorMatrix = pred, print=F)
  dat.mcar_pmm <- complete(imp_pmm, action="long", include=FALSE) # Stack completed data sets
  
  # stripplot(imp_pmm, pch = c(21, 20), cex = c(1, 1.5))
  # bwplot(imp_pmm)
  # densityplot(imp_pmm, layout = c(2, 1))

  
  ## Analysis ##
  
  # Mean for each imputed data set:
  theta <- aggregate(miete ~.imp, data = dat.mcar_pmm, mean)$miete
  # Variance of mean for each imputed data set:
  within_variance <- aggregate(miete ~.imp, data = dat.mcar_pmm, var)$miete/n
  # MI-Analysis
  MI_mean_1 <- MI.analysis(theta,within_variance,M)
  
  # Mean for each imputed data set:
  theta <- aggregate(flaeche ~.imp, data = dat.mcar_pmm, mean)$flaeche
  # Variance of mean for each imputed data set:
  within_variance <- aggregate(flaeche ~.imp, data = dat.mcar_pmm, var)$flaeche/n
  # MI-Analysis
  MI_mean_2 <- MI.analysis(theta,within_variance,M)
  
  # Save results in matrix
  MIpmm_MCAR[r, , 1] <- c(bias(true_val_miete,MI_mean_1[1]), coverage(true_val_miete, MI_mean_1[2], MI_mean_1[3]))
  MIpmm_MCAR[r, , 2] <- c(bias(true_val_flaeche,MI_mean_2[1]), coverage(true_val_flaeche, MI_mean_2[2], MI_mean_2[3]))
  
  ### Midastouch
  
  meth <- ini$method
  meth[c("miete","flaeche")] <- "midastouch"  
  
  imp_midas <- mice(dat.mcar, m=M, maxit=runs, predictorMatrix = pred, print=F)
  dat.mcar_midas <- complete(imp_midas, action="long", include=FALSE) 
  
  # plot(imp_pmm)
  
  ## Analysis ##
  
  # Mean for each imputed data set:
  theta <- aggregate(miete ~.imp, data = dat.mcar_midas, mean)$miete
  # Variance 
  within_variance <- aggregate(miete ~.imp, data = dat.mcar_midas, var)$miete/n
  # MI-Analysis
  MI_mean_1 <- MI.analysis(theta,within_variance,M)
  
  # Mean for each imputed data set:
  theta <- aggregate(flaeche ~.imp, data = dat.mcar_midas, mean)$flaeche
  # Variance 
  within_variance <- aggregate(flaeche ~.imp, data = dat.mcar_midas, var)$flaeche/n
  # MI-Analysis
  MI_mean_2 <- MI.analysis(theta,within_variance,M)
  
  # Save results in array
  MImidas_MCAR[r, , 1] <-  c(bias(true_val_miete,MI_mean_1[1]), coverage(true_val_miete, MI_mean_1[2], MI_mean_1[3]))
  MImidas_MCAR[r, , 2] <-  c(bias(true_val_flaeche,MI_mean_2[1]), coverage(true_val_flaeche, MI_mean_2[2], MI_mean_2[3]))
  
  }





##### Missing at random #####

BD_MAR <- MInorm_MAR <- MIpmm_MAR <- MImean_MAR <- MImidas_MAR <-  array(dim=c(R, 2, 2))



for (r in 1:R) {
  
  # Subsample of 500 observation
  sample_dat <- sample(nrow(dat_2), n)
  dat <- dat_2[sample_dat,]
  
  # MAR
  
  dat.mar <- dat
  
  # missingness in "miete" depend on the observed variable "lage"
  mis.mar_1 <- order(dat$lage + rnorm(n, 0, 0.5*sd(dat$lage)), 
                   decreasing=TRUE)[1:(round(0.3*n))] 
  
  # missingness in "flaeche" depend on the observed variable "zh"
  mis.mar_2 <- order(dat$zh + rnorm(n, 0, 0.5*sd(dat$zh)), 
                   decreasing=TRUE)[1:(round(0.3*n))] 
  
  # missingness depends on first variable + random component
  is.na(dat.mar[mis.mar_1, 1]) <- TRUE
  is.na(dat.mar[mis.mar_2, 2]) <- TRUE
  
  
  # sum(is.na(dat.mar))/nrow(dat.mar)
 
  
  #### Before deletion
  
  # mean analysis:
  mean_est_miete <- mean(dat$miete)
  mean_est_flaeche <- mean(dat$flaeche)
  
  ci <- t.test(dat$miete)$conf.int
  ci_low_1 <- ci[1]
  ci_upp_1 <- ci[2]
  
  ci <- t.test(dat$flaeche)$conf.int
  ci_low_2 <- ci[1]
  ci_upp_2 <- ci[2]
  
  BD_MAR[r, , 1] <- c(bias(true_val_miete,mean_est_miete),coverage(true_val_miete, ci_low_1, ci_upp_1))
  BD_MAR[r, , 2] <- c(bias(true_val_flaeche,mean_est_flaeche),coverage(true_val_flaeche, ci_low_2, ci_upp_2))
  
  
  
  #### Multiple imputation and Analysis
  
  ini <- mice(dat.mar,m=1,maxit=0)
  pred <- quickpred(dat.mar)
  
  
  
  ### Mean Imputation 
  
  meth <- ini$method
  meth[c("miete","flaeche")] <- "mean"
  
  imp_mean <- mice(dat.mar,m=M,maxit=runs, method = meth, predictorMatrix = pred, print=FALSE) 
  dat.mar_mean <- complete(imp_mean,  action="long", include=FALSE)
  
  ## Analysis
  
  # Mean for each imputed data set:
  theta <- aggregate(miete ~.imp, data=dat.mar_mean, mean)$miete
  # Variance of mean for each imputed data set:
  within_variance <- aggregate(miete ~.imp, data = dat.mar_mean, var)$miete/n
  # MI-Analysis
  MI_mean_1 <- MI.analysis(theta,within_variance,M)
  
  # Mean for each imputed data set:
  theta <- aggregate(flaeche ~.imp, data=dat.mar_mean, mean)$flaeche
  # Variance of mean for each imputed data set:
  within_variance <- aggregate(flaeche ~.imp, data = dat.mar_mean, var)$flaeche/n
  # MI-Analysis
  MI_mean_2 <- MI.analysis(theta,within_variance,M)
  
  # Save results in matrix
  MImean_MAR[r, , 1] <-  c(bias(true_val_miete,MI_mean_1[1]), coverage(true_val_miete, MI_mean_1[2], MI_mean_1[3]))
  MImean_MAR[r, , 2] <-  c(bias(true_val_flaeche,MI_mean_2[1]), coverage(true_val_flaeche, MI_mean_2[2], MI_mean_2[3]))
  
  
  
  ### Norm (BLR) Imputation
  
  meth <- ini$method
  meth[c("miete","flaeche")] <- "norm"
  
  imp_norm <- mice(dat.mar, m=M, maxit=runs, method=meth, predictorMatrix = pred, print=F)
  dat.mar_norm <- complete(imp_norm, action="long", include=FALSE) # Stack completed data sets
  
  ## Analysis ##
  
  # Mean for each imputed data set:
  theta <- aggregate(miete ~.imp, data = dat.mar_norm, mean)$miete
  # Variance of mean for each imputed data set:
  within_variance <- aggregate(miete ~.imp, data = dat.mar_norm, var)$miete/n
  # MI-Analysis
  MI_mean_1 <- MI.analysis(theta,within_variance,M)
  
  # Mean for each imputed data set:
  theta <- aggregate(flaeche ~.imp, data = dat.mar_norm, mean)$flaeche
  # Variance of mean for each imputed data set:
  within_variance <- aggregate(flaeche ~.imp, data = dat.mar_norm, var)$flaeche/n
  # MI-Analysis
  MI_mean_2 <- MI.analysis(theta,within_variance,M)
  
  # Save results in matrix
  MInorm_MAR[r, , 1] <-  c(bias(true_val_miete,MI_mean_1[1]), coverage(true_val_miete, MI_mean_1[2], MI_mean_1[3]))
  MInorm_MAR[r, , 2] <-  c(bias(true_val_flaeche,MI_mean_2[1]), coverage(true_val_flaeche, MI_mean_2[2], MI_mean_2[3]))
  
  
  
  ### PMM Imputation
  
  meth <- ini$method
  meth[c("miete","flaeche")] <- "pmm"
  
  imp_pmm <- mice(dat.mar, m=M, maxit=runs, predictorMatrix = pred, print=F)
  plot(imp_pmm)
  dat.mar_pmm <- complete(imp_pmm, action="long", include=FALSE) 
  

  ## Analysis ##
  
  # Mean for each imputed data set:
  theta <- aggregate(miete ~.imp, data = dat.mar_pmm, mean)$miete
  # Variance 
  within_variance <- aggregate(miete ~.imp, data = dat.mar_pmm, var)$miete/n
  # MI-Analysis
  MI_mean_1 <- MI.analysis(theta,within_variance,M)
  
  # Mean for each imputed data set:
  theta <- aggregate(flaeche ~.imp, data = dat.mar_pmm, mean)$flaeche
  # Variance 
  within_variance <- aggregate(flaeche ~.imp, data = dat.mar_pmm, var)$flaeche/n
  # MI-Analysis
  MI_mean_2 <- MI.analysis(theta,within_variance,M)
  
  # Save results in array
  MIpmm_MAR[r, , 1] <-  c(bias(true_val_miete,MI_mean_1[1]), coverage(true_val_miete, MI_mean_1[2], MI_mean_1[3]))
  MIpmm_MAR[r, , 2] <-  c(bias(true_val_flaeche,MI_mean_2[1]), coverage(true_val_flaeche, MI_mean_2[2], MI_mean_2[3]))
  

  ### Midastouch Imputation
  
  meth <- ini$method
  meth[c("miete","flaeche")] <- "midastouch"  
  
  imp_midas <- mice(dat.mar, m=M, maxit=runs, predictorMatrix = pred, print=F)
  dat.mar_midas <- complete(imp_midas, action="long", include=FALSE) 
  
  # plot(imp_pmm)
  
  ## Analysis ##
  
  # Mean for each imputed data set:
  theta <- aggregate(miete ~.imp, data = dat.mar_midas, mean)$miete
  # Variance 
  within_variance <- aggregate(miete ~.imp, data = dat.mar_midas, var)$miete/n
  # MI-Analysis
  MI_mean_1 <- MI.analysis(theta,within_variance,M)
  
  # Mean for each imputed data set:
  theta <- aggregate(flaeche ~.imp, data = dat.mar_midas, mean)$flaeche
  # Variance 
  within_variance <- aggregate(flaeche ~.imp, data = dat.mar_midas, var)$flaeche/n
  # MI-Analysis
  MI_mean_2 <- MI.analysis(theta,within_variance,M)
  
  # Save results in array
  MImidas_MAR[r, , 1] <-  c(bias(true_val_miete,MI_mean_1[1]), coverage(true_val_miete, MI_mean_1[2], MI_mean_1[3]))
  MImidas_MAR[r, , 2] <-  c(bias(true_val_flaeche,MI_mean_2[1]), coverage(true_val_flaeche, MI_mean_2[2], MI_mean_2[3]))
  
  }




### Diagnostics ###

Bias <- Coverage <- array(dim=c(2, 5, 2))

## MCAR
# Bias: miete
Bias[1,1,1] <- mean(BD_MCAR[,1,1]) 
Bias[1,2,1] <- mean(MImean_MCAR[,1,1]) 
Bias[1,3,1] <- mean(MInorm_MCAR[,1,1])
Bias[1,4,1] <- mean(MIpmm_MCAR[,1,1])
Bias[1,5,1] <- mean(MImidas_MCAR[,1,1])

# Bias: flaeche
Bias[1,1,2] <- mean(BD_MCAR[,1,2])
Bias[1,2,2] <- mean(MImean_MCAR[,1,2])
Bias[1,3,2] <- mean(MInorm_MCAR[,1,2]) 
Bias[1,4,2] <- mean(MIpmm_MCAR[,1,2]) 
Bias[1,5,2] <- mean(MImidas_MCAR[,1,2])

## MAR
# Bias: miete
Bias[2,1,1] <- mean(BD_MAR[,1,1]) 
Bias[2,2,1] <- mean(MImean_MAR[,1,1]) 
Bias[2,3,1] <- mean(MInorm_MAR[,1,1])
Bias[2,4,1] <- mean(MIpmm_MAR[,1,1])
Bias[2,5,1] <- mean(MImidas_MAR[,1,1])

# Bias: flaeche
Bias[2,1,2] <- mean(BD_MAR[,1,2])
Bias[2,2,2] <- mean(MImean_MAR[,1,2])
Bias[2,3,2] <- mean(MInorm_MAR[,1,2]) 
Bias[2,4,2] <- mean(MIpmm_MAR[,1,2]) 
Bias[2,5,2] <- mean(MImidas_MAR[,1,2])


## MCAR
# Coverage: miete
Coverage[1,1,1] <- mean(BD_MCAR[,2,1]) 
Coverage[1,2,1] <- mean(MImean_MCAR[,2,1]) 
Coverage[1,3,1] <- mean(MInorm_MCAR[,2,1]) 
Coverage[1,4,1] <- mean(MIpmm_MCAR[,2,1])
Coverage[1,5,1] <- mean(MImidas_MCAR[,2,1])

# Coverage: flache
Coverage[1,1,2] <- mean(BD_MCAR[,2,2]) 
Coverage[1,2,2] <- mean(MImean_MCAR[,2,2]) 
Coverage[1,3,2] <- mean(MInorm_MCAR[,2,2]) 
Coverage[1,4,2] <- mean(MIpmm_MCAR[,2,2]) 
Coverage[1,5,2] <- mean(MImidas_MCAR[,2,2])


## MAR
# Coverage: miete
Coverage[2,1,1] <- mean(BD_MAR[,2,1]) 
Coverage[2,2,1] <- mean(MImean_MAR[,2,1]) 
Coverage[2,3,1] <- mean(MInorm_MAR[,2,1]) 
Coverage[2,4,1] <- mean(MIpmm_MAR[,2,1])
Coverage[2,5,1] <- mean(MImidas_MAR[,2,1])


# Coverage: flaeche
Coverage[2,1,2] <- mean(BD_MAR[,2,2]) 
Coverage[2,2,2] <- mean(MImean_MAR[,2,2]) 
Coverage[2,3,2] <- mean(MInorm_MAR[,2,2]) 
Coverage[2,4,2] <- mean(MIpmm_MAR[,2,2]) 
Coverage[2,5,2] <- mean(MImidas_MAR[,2,2]) 


colnames(Bias) <- colnames(Coverage) <- c("bd","mean","norm","pmm","mida")
rownames(Bias) <- rownames(Coverage) <- c("MCAR_mean","MAR_mean")
dimnames(Bias)[[3]] <- dimnames(Coverage)[[3]]  <- c("Miete", "Flaeche")

Bias
Coverage


#### Graphical diagnostics ####


###########################
## PMM

windows()
par(mfrow=c(1,4))
hist(dat$miete, col="slategray4", freq = T, main = "BD", xlab = "Miete")
hist(dat.mar_mean$miete, col="slategray3", freq = T, main = "Mean", xlab = "Miete")
hist(dat.mar_norm$miete, col="slategray2", freq = T, main = "Norm", xlab = "Miete")
hist(dat.mar_pmm$miete, col="slategray1", freq = T, main = "PMM", xlab = "Miete")
dev.off()
###########################

### Boxplot
windows()
par(mfrow=c(2,2))
## MCAR
# miete
boxplot(MImean_MCAR[,1,1], MInorm_MCAR[,1,1], MIpmm_MCAR[,1,1], MImidas_MCAR[,1,1], 
        col=c("slategray4","slategray3","slategray2","slategray1"),
        main="Comparison: rent", xlab="methods", ylab="Relative Bias", names = c("mean","norm", "pmm", "mida"))
abline(h=0, lty = 5, col="red")

# flaeche
boxplot(MImean_MCAR[,1,1], MInorm_MCAR[,1,2], MIpmm_MCAR[,1,2], MImean_MCAR[,1,2], 
        col=c("slategray4","slategray3","slategray2","slategray1"),
        main="Comparison: rental space", xlab="methods", ylab="Relative Bias", names = c("mean","norm", "pmm", "mida"))
abline(h=0, lty = 5, col="red")


## MAR
# miete
boxplot(MImean_MAR[,1,1], MInorm_MAR[,1,1], MIpmm_MAR[,1,1], MImidas_MAR[,1,1], 
        col=c("slategray4","slategray3","slategray2","slategray1"),
        main="Comparison: rent", xlab="methods", ylab="Relative Bias", names = c("mean","norm", "pmm", "mida"))
abline(h=0, lty = 5, col="red")

# flaeche
boxplot(MImean_MAR[,1,1], MInorm_MAR[,1,2], MIpmm_MAR[,1,2], MImean_MAR[,1,2], 
        col=c("slategray4","slategray3","slategray2","slategray1"),
        main="Comparison: rental space", xlab="methods", ylab="Relative Bias", names = c("mean","norm", "pmm", "mida"))
abline(h=0, lty = 5, col="red")
dev.off()


### Barplot
?barplot
windows()
par(mfrow=c(2,2))
## MCAR
# miete
barplot(c(mean(MImean_MCAR[,2,1]),mean(MInorm_MCAR[,2,1]),mean(MIpmm_MCAR[,2,1]),mean(MImidas_MCAR[,2,1])),
        main = "rent under MCAR",col=c("slategray4","slategray3","slategray2","slategray1"),
        beside = TRUE,ylim = c(0,1),xlab ="methods",ylab = "coverage",names = c("mean","norm", "pmm", "mida"))
abline(h=0.95,col="green",lwd=2)

# flaeche
barplot(c(mean(MImean_MCAR[,2,2]),mean(MInorm_MCAR[,2,2]),mean(MIpmm_MCAR[,2,2]),mean(MImidas_MCAR[,2,2])),
        main = "rental space under MCAR",col=c("slategray4","slategray3","slategray2","slategray1"),
        beside = TRUE,ylim = c(0,1),xlab ="methods",ylab = "coverage",names = c("mean","norm", "pmm", "mida"))
abline(h=0.95,col="green",lwd=2)


## MAR
# miete
barplot(c(mean(MImean_MAR[,2,1]),mean(MInorm_MAR[,2,1]),mean(MIpmm_MAR[,2,1]),mean(MImidas_MAR[,2,1])),
        main = "rent under MAR",col=c("slategray4","slategray3","slategray2","slategray1"),
        beside = TRUE,ylim = c(0,1),xlab ="methods",ylab = "coverage",names = c("mean","norm", "pmm", "mida"))
abline(h=0.95,col="green",lwd=2)

# flaeche
barplot(c(mean(MImean_MAR[,2,2]),mean(MInorm_MAR[,2,2]),mean(MIpmm_MAR[,2,2]),mean(MImidas_MAR[,2,2])),
        main = "rental space under MAR",col=c("slategray4","slategray3","slategray2","slategray1"),
        beside = TRUE,ylim = c(0,1),xlab ="methods",ylab = "coverage",names = c("mean","norm", "pmm", "mida"))
abline(h=0.95,col="green",lwd=2)
dev.off()

#########################################
attach(dat.mcar)
a <- 0.1
b <- 0.2
c <- 0.3
d <- 4

x1 <- dat.mcar$flaeche
x2 <- dat.mcar$bjahr
x3 <- dat.mcar$zh
x4 <- dat.mcar$mieteqm
eps <- rchisq(n,107/96)
y <-a*x1^2 + b*x2 + c*as.numeric(x3) + eps 
y
plot(x1,y)
lin.mod <- lm(miete~ x1^2 + x2 + x3,data = dat.mcar)
abline(lin.mod)
##
g <- a*mietspiegel$flaeche^2+ b*mietspiegel$bjahr + c*as.numeric(mietspiegel$zh)
g
plot(mietspiegel$flaeche,g)
lin.mod2 <- lm(miete~ mietspiegel$flaeche^2 + mietspiegel$bjahr + mietspiegel$zh,data = mietspiegel)
abline(lin.mod2)

