# Ustawienie prywatnego katalogu na pakiety R w środowisku Streamlit Cloud
user_lib <- Sys.getenv("R_LIBS_USER")
if (!dir.exists(user_lib)) dir.create(user_lib, recursive = TRUE)
.libPaths(c(user_lib, .libPaths()))

# Lista wymaganych pakietów
needed_packages <- c("bmrm", "dplyr", "future", "future.apply")
missing_packages <- needed_packages[!(needed_packages %in% installed.packages()[,"Package"])]

# Automatyczna instalacja bez pytania o interaktywne potwierdzenie
if(length(missing_packages) > 0) {
  install.packages(missing_packages, 
                   lib = user_lib, 
                   repos = "https://cloud.r-project.org", 
                   dependencies = TRUE)
}

library(bmrm)
library(dplyr)
library(future)
library(future.apply)

# Wczytanie argumentów wiersza poleceń
args <- commandArgs(trailingOnly = TRUE)

n           <- as.integer(args[1])
m           <- as.integer(args[2])
p           <- as.integer(args[3])
rho         <- as.numeric(args[4])
reps        <- as.integer(args[5])
model_type  <- args[6]
beta_type   <- args[7]
output_file <- args[8]

plan(multisession)

if(beta_type == "betaA") {
  beta <- c(3, 1.5, 0, 0, 2, rep(0, max(0, p - 5)))[1:p]
} else {
  beta <- c(rep(1, 7), rep(-1, 3), rep(0, max(0, p - 10)))[1:p]
}

Covar <- outer(1:p, 1:p, function(i, j) rho^abs(i - j))
Ch <- t(chol(Covar))

results_list <- future_lapply(1:reps, function(r) {
  XX_q <- matrix(rnorm(m * p), m) %*% t(Ch)
  if(model_type == "M1") YY_q <- XX_q %*% beta + rnorm(m)
  else if(model_type == "M2") YY_q <- XX_q %*% beta + rcauchy(m)
  else YY_q <- exp(XX_q %*% beta) + rnorm(m)
  
  cuts <- quantile(YY_q, probs = c(0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875))
  
  Xtrain <- matrix(rnorm(n * p), n) %*% t(Ch)
  Xtest  <- matrix(rnorm(m * p), m) %*% t(Ch)
  
  if(model_type == "M1") {
    Ytrain1 <- Xtrain %*% beta + rnorm(n)
    Ytest1  <- Xtest %*% beta + rnorm(m)
  } else if(model_type == "M2") {
    Ytrain1 <- Xtrain %*% beta + rcauchy(n)
    Ytest1  <- Xtest %*% beta + rcauchy(m)
  } else {
    Ytrain1 <- exp(Xtrain %*% beta) + rnorm(n)
    Ytest1  <- exp(Xtest %*% beta) + rnorm(m)
  }
  
  Ytrain <- cut(Ytrain1, c(-Inf, cuts, Inf), labels = FALSE)
  Ytest  <- cut(Ytest1, c(-Inf, cuts, Inf), labels = FALSE)
  
  mu <- colMeans(Xtrain)
  sdv <- apply(Xtrain, 2, sd)
  sdv[sdv == 0] <- 1
  Xtrain <- scale(Xtrain, mu, sdv)
  Xtest  <- scale(Xtest, mu, sdv)
  
  ij <- expand.grid(i = 1:n, j = 1:n)
  ijtest <- expand.grid(i = 1:m, j = 1:m)
  btrain <- which(Ytrain[ij$i] - Ytrain[ij$j] > 0)
  btest  <- which(Ytest[ijtest$i] - Ytest[ijtest$j] > 0)
  
  lossfun <- ordinalRegressionLoss(Xtrain, Ytrain, impl = "loglin")
  
  # Ridge
  beta_R <- nrbm(lossfun, LAMBDA = 0.001, w0 = rep(0, p), EPSILON_TOL = 0.01, MAX_ITER = 100L)
  
  # LASSO
  a <- 2^(-10:10)
  lambda_grid <- a * sqrt(log(p)/n)
  nlambda <- length(lambda_grid)
  U <- matrix(0, nlambda, p)
  BIC <- numeric(nlambda)
  
  for(i in 1:nlambda){
    wstart <- if(i == 1) rep(0, p) else U[i-1,]
    U[i,] <- nrbmL1(lossfun, LAMBDA = lambda_grid[i], w0 = wstart, EPSILON_TOL = 0.01, MAX_ITER = 100L)
    ftrain <- as.vector(Xtrain %*% U[i,])
    diff <- ftrain[ij$i] - ftrain[ij$j]
    loss <- sum(pmax(0, 1 - diff[btrain])) / (n * (n - 1))
    df <- sum(abs(U[i,]) != 0)
    BIC[i] <- loss + log(n)/(2*n)*df
  }
  bestL <- which.min(BIC)
  beta_L <- U[bestL,]
  
  # Adaptive LASSO
  active <- which(beta_L != 0)
  beta_AL <- rep(0, p)
  if(length(active) > 0){
    weights <- 1 / abs(beta_L[active])
    XtrainAL <- sweep(Xtrain[, active, drop = FALSE], 2, weights, "/")
    lossfunAL <- ordinalRegressionLoss(XtrainAL, Ytrain, impl = "loglin")
    UAL <- matrix(0, nlambda, ncol(XtrainAL))
    BICAL <- numeric(nlambda)
    
    for(i in 1:nlambda){
      wstart <- if(i == 1) rep(0, ncol(XtrainAL)) else UAL[i-1,]
      UAL[i,] <- nrbmL1(lossfunAL, LAMBDA = lambda_grid[i], w0 = wstart, EPSILON_TOL = 0.01, MAX_ITER = 100L)
      ftrain <- as.vector(XtrainAL %*% UAL[i,])
      diff <- ftrain[ij$i] - ftrain[ij$j]
      loss <- sum(pmax(0, 1 - diff[btrain])) / (n * (n - 1))
      df <- sum(abs(UAL[i,]) != 0)
      BICAL[i] <- loss + log(n)/(2*n)*df
    }
    bestAL <- which.min(BICAL)
    beta_AL[active] <- UAL[bestAL,] / weights
  }
  
  ranking_acc <- function(X, b) {
    f <- as.vector(X %*% b)
    diff <- f[ijtest$i] - f[ijtest$j]
    sum(diff[btest] > 0) / length(btest)
  }
  
  err_R <- 1 - ranking_acc(Xtest, beta_R)
  err_L <- 1 - ranking_acc(Xtest, beta_L)
  err_AL <- 1 - ranking_acc(Xtest, beta_AL)
  
  b_true_n <- beta / sqrt(sum(beta^2))
  norm_R <- sqrt(sum((beta_R / sdv)^2))
  norm_L <- sqrt(sum((beta_L / sdv)^2))
  norm_AL <- sqrt(sum((beta_AL / sdv)^2))
  
  b_Rn <- if(norm_R > 0) (beta_R / sdv) / norm_R else rep(0, p)
  b_Ln <- if(norm_L > 0) (beta_L / sdv) / norm_L else rep(0, p)
  b_ALn <- if(norm_AL > 0) (beta_AL / sdv) / norm_AL else rep(0, p)
  
  est_R <- sqrt(sum((b_Rn - b_true_n)^2))
  est_L <- sqrt(sum((b_Ln - b_true_n)^2))
  est_AL <- sqrt(sum((b_ALn - b_true_n)^2))
  
  list(
    pred = data.frame(model = c("Ridge", "LASSO", "Adaptive LASSO"), error = c(err_R, err_L, err_AL)),
    est  = data.frame(model = c("Ridge", "LASSO", "Adaptive LASSO"), error = c(est_R, est_L, est_AL))
  )
}, future.seed = TRUE)

boxplot_prediction <- bind_rows(lapply(results_list, function(x) x$pred))
boxplot_estimation <- bind_rows(lapply(results_list, function(x) x$est))

save(boxplot_prediction, boxplot_estimation, file = output_file)
