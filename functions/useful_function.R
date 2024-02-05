#' Load required libraries
#'
#' This function loads a set of libraries commonly used in statistical
#' analysis and data manipulation.
#'
#' @return None
#' @import data.table
#' @import roahd
#' @import twinning
#' @import progress
#' @import plotrix
#' @import car
#' @import glmnet
#' @import caret
#' @import MASS
#' @import fda
#' @import KernSmooth
#' @import rgl
#' @import fields
#'
#' @examples
#' load_libraries()
#'
load_libraries = function(){
  library(data.table)
  library(roahd)
  library(twinning)
  library(progress)
  library("plotrix") 
  library(data.table)
  library(car)
  library(glmnet)
  library(caret)
  library(MASS)
  library(fda)
  library(KernSmooth)
  library(rgl)
  library(fields)
  library(GGally)
  library(ggplot2)
  library(maps)
  library(cowplot)
  library(sp)           
  library(lattice)      
  library(gstat)
  library(fdagstat)
  library(expm)
  library(readr)
  library(fdatest)
}

#' Prepare the dataset
#'
#' This function prepares the dataset by performing various operations, 
#' including converting certain columns to factors and removing rows with 
#' missing values.
#'
#' @param dataset The dataset to be prepared.
#' @return The prepared dataset.
#' @export
#' 
#' @examples
#' dataset <- data_prepare(dataset)
#'
data_preparation = function(dataset){
  
  # Convert specified columns to factors
  dataset$event_id = as.factor(dataset$event_id)
  dataset$ev_nation_code = as.factor(dataset$ev_nation_code)
  dataset$sof = as.factor(dataset$sof)
  dataset$station_id = as.factor(dataset$station_id)
  
  # Remove rows where the magnitude column is missing
  dataset = dataset[!is.na(dataset$magnitude), ]
  
  # Remove rows where the ev_nation_code column is missing
  dataset = dataset[!is.na(dataset$ev_nation_code), ]
  
  # Remove rows where the SA_0 column is missing
  dataset = dataset[!is.na(dataset$SA_0), ]
  
  # Return the modified dataset
  return(dataset)
}


#' Hypothesis Testing with Permutation Method
#'
#' This function performs hypothesis testing using the permutation method to compare means
#' between two subpopulations based on a magnitude threshold.
#'
#' @param SA_subpop1 Data for the first subpopulation.
#' @param SA_subpop2 Data for the second subpopulation.
#' @param mh0 Threshold value for magnitude.
#' @param periods Time periods.
#' @param B Number of permutations.
#' @param seed Seed for randomization.
#'
#' @return A list containing the observed test statistic (`t_obs`)
#'         and the permutational distribution (`t_dist`).
#'
#' @examples
#' result <- anova_magnitude_testing(SA_subpop1, SA_subpop2, mh0, periods, B, seed)
#'
anova_magnitude_testing = function(SA_subpop1, SA_subpop2, mh0, periods, B, seed){
  # Extract mean values
  t1.mean = colMeans(SA_subpop1)  
  t2.mean = colMeans(SA_subpop2)
  
  # Reference line
  Mh <- rep(mh0, times = 36)
  
  # Plotting the difference between means for each group
  matplot(log10(periods[2:37]), t(rbind(t1.mean - Mh, t2.mean - Mh)), type='l', col=c('#FFC700', '#9D260C'), lty=1, xlab='Log(T)', ylab='Magnitude')
  legend('topright', c('Magnitude < 5.7', 'Magnitude > 5.7'), col=c('#FFC700', '#9D260C'), lty=1, lwd=2)
  
  # Compute the test statistics
  n1 = dim(SA_subpop1)[1]
  n2 = dim(SA_subpop2)[1]
  n  = n1 + n2
  
  t1 <- as.numeric(t1.mean) - Mh
  t2 <- as.numeric(t2.mean) - Mh
  
  T_test = as.numeric(t(t1 - t2) %*% (t1 - t2))
  
  # Estimating the permutational distribution T2 under H0
  T2 = numeric(B)
  set.seed(seed)
  for (perm in 1:B){
    # Permute the rows of the data matrix
    t_pooled = rbind(SA_subpop1, SA_subpop2)
    permutation = sample(n)
    t_perm = t_pooled[permutation,]
    t1_perm = t_perm[1:n1,]
    t2_perm = t_perm[(n1 + 1):n,]
    
    # Evaluation of the test statistic on permuted data
    t1.mean_perm = colMeans(t1_perm)
    t2.mean_perm = colMeans(t2_perm)
    T2[perm]  = t(t1.mean_perm - t2.mean_perm) %*% (t1.mean_perm - t2.mean_perm)
  }
  
  # Compute the test p-value
  p_val = sum(T2 >= T_test)/B
  
  return(list(t_obs = T_test, t_dist = T2, p_val = p_val))
}


#' Hypothesis Testing with Permutation Method for MANOVA
#'
#' This function performs hypothesis testing using the permutation 
#' method for Multivariate Analysis of Variance (MANOVA).
#'
#' @param dataset The dataset for MANOVA.
#' @param fit_manova The MANOVA fit obtained using the `manova` function.
#' @param B Number of permutations.
#' @param seed Seed for randomization.
#'
#' @return A list containing the observed test statistic (`t_obs`), 
#'         the permutational distribution (`t_dist`)
#'         and the p-value (`t_p_val`).
#'
#' @examples
#' result <- manova_sof_testing(dataset, fit_manova, B, seed)
#'
manova_sof_testing = function(dataset, fit_manova, B, seed){
  
  # Extract the test statistic
  T0 <- -summary.manova(fit_manova, test = 'Wilks')$stats[1, 2]
  
  # Extract indices for different levels of the factor 'sof'
  i1 <- which(dataset$sof == 'SS')
  i2 <- which(dataset$sof == 'TF')
  i3 <- which(dataset$sof == 'NF')
  n1 <- length(i1)
  n2 <- length(i2)
  n3 <- length(i3)
  n  <- n1 + n2 + n3
  
  # Log-transformed SA data
  SA_log = log10(dataset[, 12:47])
  
  # Initialize an empty vector for the test statistics
  T_stat <- numeric(B) 
  
  set.seed(seed)
  for(perm in 1:B){
    # Permutation:
    permutation <- sample(1:n)
    SA_perm <- as.matrix(SA_log[permutation, ])
    fit_perm <- manova(SA_perm ~ dataset$sof)
    
    # Test statistic:
    T_stat[perm] <- -summary.manova(fit_perm, test = 'Wilks')$stats[1, 2]
  }
  
  # Compute the p-value
  p_val <- sum(T_stat >= T0) / B
  
  # Extract sub-populations for plotting
  SA_subpop1 <- dataset[which(dataset$sof == 'SS'), 12:47]
  SA_subpop2 <- dataset[which(dataset$sof == 'TF'), 12:47]   
  SA_subpop3 <- dataset[which(dataset$sof == 'NF'), 12:47]
  
  # Calculate mean values for each sub-population
  t1.mean = colMeans(SA_subpop1)
  t2.mean = colMeans(SA_subpop2)
  t3.mean = colMeans(SA_subpop3)
  
  # Plot the means
  matplot(logperiods, t(rbind(t1.mean, t2.mean, t3.mean)), type='l', col=c("#FFC700", "#FF9900", "#9D260C"), lty=1, xlab='Log(T)', ylab='log(SA)')
  legend('topleft', c('SS', 'TF', 'NF'), col=c("#FFC700", "#FF9900", "#9D260C"), lty=1, lwd=2)
  
  # Return the results as a list
  return(list(t_obs = T0, t_dist = T_stat, p_val = p_val))
}


#' Build Design Matrix
#'
#' This function constructs a design matrix based on the given dataset and specified parameters.
#'
#' @param dataset The input dataset.
#' @param Mh Threshold value for magnitude.
#' @param Mref Reference value for magnitude.
#' @param with_sof Logical indicating whether to include the 'sof' variable.
#'
#' @return A list containing the design matrix (`design_matrix`) and 
#'         all its components.
#'
#' @examples
#' design_matrix <- build_design_matrix(dataset, Mh, Mref, with_sof)
#'
build_design_matrix = function(dataset, Mh, Mref, with_sof){
  
  # Remove rows with distance == 0 and save the new dimension
  dataset <- dataset[-which(dataset$distance == 0), ]
  n = dim(dataset)[1]
  
  # Source function
  source <- dataset$magnitude - Mh
  
  # Dummy variables
  dummy1 <- dummy2 <- rep(0, n)
  for (i in c(1:n)){
    if (dataset$magnitude[i] <= Mh){
      dummy1[i] <- 1
    } else { dummy2[i] <- 1 }
  }
  
  # SoF
  if(with_sof){
    sof = as.numeric(dataset$sof)
  }else{
    sof = NULL
  }
  
  # Path functions
  # c1
  path1 <- (dataset$magnitude - Mref) * log(dataset$distance, base=10)
  # c2
  path2 <- log(dataset$distance, base = 10) 
  
  # Site function
  site <- rep(0, n)
  site[which(dataset$vs30 < 1500)] <- log(dataset[which(dataset$vs30 < 1500), 4]/800)
  site[which(dataset$vs30 >= 1500)] <- log(1500/800)
  
  # Distance
  distance <- dataset$distance
  
  # Construct the design matrix
  if (with_sof){
    Z <- data.frame(
      source, 
      dummy1, dummy2, 
      sof, 
      distance, 
      path1, path2, 
      site, 
      dataset[, 5:40]
    )
  } else {
    Z <- data.frame(
      source, 
      dummy1, dummy2,
      distance, 
      path1, path2, 
      site, 
      dataset[, 5:40]
    )
  }
  
  res = list(
    design_matrix = Z,
    source = source,
    dummy1 = dummy1,
    dummy2 = dummy2,
    sof = sof,
    distance = distance, 
    path1 = path1, 
    path2 = path2, 
    site = site, 
    log10_SA = log10(dataset[, 5:40]),
    n = dim(dataset)[1]
  )
  
  return(res)
}


#' Bootstrap Confidence Intervals for Regression Coefficients
#'
#' This function performs bootstrapping to compute confidence intervals for regression coefficients.
#'
#' @param design_matrix The design matrix for the regression.
#' @param names Names of the regression coefficients.
#' @param alpha Significance level.
#' @param B Number of bootstrap samples.
#' @param seed Seed for randomization.
#' @param with_sof Logical indicating whether 'sof' variable is included.
#'
#' @return A list of bootstrap confidence intervals for each regression coefficient.
#'
#' @examples
#' res = bootstrap_CI(design_matrix, names, alpha, B, seed, with_sof)
#'
bootstrap_CI = function(design_matrix, names, alpha, B, seed, with_sof){
  
  # Time periods
  periods = c(
    0, 0.01, 0.025, 0.04, 0.05, 0.07, 0.1, 0.15, 0.2, 0.25, 0.3,
    0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.75, 0.8, 0.9, 1, 1.2, 1.4,
    1.6, 1.8, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10
  )
  
  # Number of predictors
  n_coeffs = length(names)
  p = dim(design_matrix)[2]
  
  # Set up plot layout
  par(mfrow = c(3, 3), mar = c(3, 3, 2, 1))
  
  results = list()
  
  for (j in c(1:n_coeffs)){
    
    CI.RP.L <- c()
    
    for(i in c(n_coeffs:p)){
      # Build the linear model
      if(with_sof){
        model <- lm(
          design_matrix[, i] ~ 
            source:dummy1 + source:dummy2 + as.factor(sof) +
            distance + path1 + path2 + site, 
          data = design_matrix[, c(1, 2, 3, 4, 5, 6, 7, 8, i)]
        )
      } else {
        model <- lm(
          design_matrix[, i] ~ 
            source:dummy1 + source:dummy2 +
            distance + path1 + path2 + site, 
          data = design_matrix[, c(1, 2, 3, 4, 5, 6, 7, i)]
        )
      }
      
      fitted.obs <- predict(model)
      res.obs <- design_matrix[, i] - fitted.obs
      
      L.obs <- summary(model)$coefficients[j, 1]
      T.boot.L <- numeric(B)
      set.seed(seed)
      
      # Bootstrap loop
      for (b in 1:B) {
        response.b <- fitted.obs + sample(res.obs, replace = TRUE)
        if(with_sof){
          dataset <- cbind(design_matrix[, c(1, 2, 3, 4, 5, 6, 7, 8)], response.b)
          fm.b <- lm(response.b ~ source:dummy1 + source:dummy2 + as.factor(sof) + distance + 
                       path1 + path2 + site, data = dataset)
        } else {
          dataset <- cbind(design_matrix[, c(1, 2, 3, 4, 5, 6, 7)], response.b)
          fm.b <- lm(response.b ~ source:dummy1 + source:dummy2 + distance + 
                       path1 + path2 + site, data = dataset)
        }
        T.boot.L[b] <- summary(fm.b)$coefficients[j, 1]
      }
      
      # Calculate quantiles for confidence interval
      right.quantile.L <- quantile(T.boot.L, 1 - alpha/2)
      left.quantile.L <- quantile(T.boot.L, alpha/2)
      
      # Create a matrix for confidence intervals
      CI.RP.L <- rbind(
        CI.RP.L, 
        c(L.obs - (right.quantile.L - L.obs), L.obs, L.obs - (left.quantile.L - L.obs))
      )
    }
    
    # Plot the confidence intervals
    plotCI(x = log10(periods[2:37]), y = CI.RP.L[1:36, 2], li = CI.RP.L[1:36, 1], ui = CI.RP.L[1:36, 3],
           col = "black", pch=19, main = paste(names[j]), xlab = 'log(T)')
    abline(h = 0)
    
    results[[names[j]]] <- list(lwr = CI.RP.L[, 1], lvl = CI.RP.L[, 2], upr = CI.RP.L[, 3])
  }
  
  return(results)
}



#' Get Optimal Lambda for Regularization
#'
#' This function calculates the optimal lambda for regularization using the Generalized Cross Validation (GCV) criterion.
#'
#' @param grid0 Initial value of the exponent for lambda.
#' @param gridN Final value of the exponent for lambda.
#' @param by Step size for the exponent.
#' @param basis Basis functions.
#' @param Xobs0 Observed data.
#' @param periods Time periods.
#'
#' @return The optimal lambda value.
#'
#' @examples
#' optimal_lambda <- get_optimal_lambda(-5, 2, 0.1, my_basis, my_data)
#'
get_optimal_lambda = function(grid0, gridN, by, basis, Xobs0, periods){
  
  # Generate a sequence of lambda values
  lambda <- 10^seq(grid0, gridN, by = by)
  n = length(lambda)
  gcv <- numeric(n)
  
  # Loop to calculate GCV for each lambda
  for (i in 1:n){
    
    # Set up functional parameters for smoothing
    functionalPar <- fdPar(fdobj = basis, Lfdobj = 2, lambda = lambda[i])  
    
    # Calculate GCV
    gcv[i] <- smooth.basis(periods, Xobs0, functionalPar)$gcv
  }
  
  # Find the optimal lambda (minimum GCV)
  opt.lambda <- lambda[which.min(gcv)]
  
  # Plot GCV values over log(lambda)
  plot(log10(lambda), gcv, pch=19, cex=0.75, col='black', xlab='log10(lambda)', ylab='GCV')
  points(opt.lambda, min(gcv), col='darkred', cex=2)
  
  return(opt.lambda)
}


#' Plot Functional Coefficients
#'
#' This function plots the functional coefficients obtained from a model for visualization.
#'
#' @param mod The model object.
#' @param t0 Starting point for the time grid.
#' @param tN Ending point for the time grid.
#' @param length.out Number of points in the time grid.
#' @param names Names of the coefficients.
#'
#' @examples
#' plot_functional_coeff(my_model, 0, 10, 100, c('coeff1', 'coeff2'))
#'
plot_functional_coeff = function(mod, t0, tN, length.out, names){
  
  # Set up the plot layout
  par(mfrow = c(2, 4))
  
  # Number of coefficients
  n_coeff = length(names)
  
  # Generate the time grid
  t_grid <- seq(t0, tN, length.out = length.out)
  
  # Loop through each coefficient and plot
  for(i in 1:n_coeff){
    # Extract the functional coefficient
    beta.hat.fd <- mod$betaestlist[[i]]
    beta.vals <- eval.fd(t_grid, beta.hat.fd$fd)
    
    # Plot the estimated coefficients beta_hat
    plot(t_grid, beta.vals, type = 'l', lwd = 3, main = paste(names[i]))
  }
}


#' Residual Decorrelation Analysis
#'
#' This function performs residual decorrelation analysis to estimate and model the covariance structure of the residuals.
#'
#' @param train_set The training dataset.
#' @param mod The model object.
#' @param t0 Starting point for the time grid.
#' @param tN Ending point for the time grid.
#' @param length.out Number of points in the time grid.
#' @param vgm Variogram model for fitting.
#' @param file_to_save_res File path to save the functional residuals.
#' @param file_to_save_sigma File path to save the estimated covariance matrix.
#'
#' @examples
#' residual_decorrelation(train_data, my_model, 0, 10, 100, "residuals.csv", "sigma.csv", my_variogram)
#'
residual_decorrelation = function(train_set, mod, t0, tN, length.out, vgm, file_to_save_res, file_to_save_sigma){
  
  # Define the evaluation grid
  t_grid <- seq(t0, tN, length.out = length.out)
  
  # Compute the functional residuals
  res_fd <- mod$yfdobj - mod$yhatfdobj
  res_fd_vals <- data.frame(eval.fd(t_grid, res_fd))
  
  # Save the functional residuals to a file
  write.table(res_fd_vals, file = file_to_save_res, sep = ",", row.names = FALSE, col.names = FALSE)
  
  # Remove rows with zero distance
  train_set = train_set[-which(train_set$distance == 0), ]
  
  # Create a data frame for station coordinates
  train_set_station <- data.frame(
    station_lat = jitter(train_set$station_lat),
    station_lon = jitter(train_set$station_lon)
  )
  
  # Create a functional data object for residuals
  g <- fstat(NULL, vName = "res_fd", Coordinates = train_set_station, Functions = res_fd_vals, scalar = FALSE)
  
  # Estimate the drift
  g <- estimateDrift("~.", g, Intercept = TRUE)
  
  # Estimate the functional variogram
  g <- fvariogram("~ station_lat + station_lon", g, Nlags = 100, LagMax = 1, ArgStep = 1, comments = FALSE)
  
  # Fit the variogram model
  g <- fitVariograms(g, vgm, fitRanges = FALSE, forceNugget = FALSE)
  
  # Plot the variogram
  plotVariogram(g)
  
  # Add the estimated covariance matrix
  g <- addCovariance(g)
  
  # Extract the estimated covariance matrix
  sigma_res <- g$covariance$omni$res_fd
  
  # Save the estimated covariance matrix to a file
  write.table(sigma_res, file = file_to_save_sigma, sep = ",", row.names = FALSE, col.names = FALSE)
  
  return(list(sigma = sigma_res, res = res_fd_vals, t_grid=t_grid))
}


#' Functional Bootstrap Significance Test
#'
#' This function performs a significance test using functional bootstrap resampling.
#'
#' @param logperiods Logarithm of time periods.
#' @param basis Basis matrix.
#' @param mod The model object.
#' @param opt.lambda Optimal lambda.
#' @param res Residuals.
#' @param sqrt_root_sigma Square root of the estimated covariance matrix.
#' @param design_matrix Design matrix.
#' @param n_coeffs Number of coefficients.
#' @param sel_coeff Selected coefficient for the test.
#' @param value Values containing the grid, sigma, and residuals.
#' @param B Number of bootstrap samples.
#' @param seed Seed for randomization.
#'
#' @return A list containing the bootstrap test statistics (`t_stat`),
#'         observed test statistic (`t_obs`) and p-value (`p_val`).
#'
#' @examples
#' result <- functional_bootstrap_significance_test(logperiods, basis, my_model, 
#'    my_lambda, residuals, sqrt_sigma, design_matrix, n_coeffs, sel_coeff, my_value,
#'    B, seed)
#'
functional_bootstrap_significance_test = function(
    logperiods, basis, mod, opt.lambda, res, sqrt_root_sigma, design_matrix, n_coeffs, sel_coeff, value, B, seed
  ){
  
  # Create the basis matrix
  basismat <- eval.basis(logperiods, basis)
  
  t_grid <- value$t_grid
  sigma <- as.matrix(value$sigma)
  res <- value$res
  design_matrix <- design_matrix[, 1:n_coeffs]
  
  # Evaluate the fitted values
  y_true <- eval.fd(t_grid, mod$yfdobj)
  
  # Calculate the inverse of the covariance matrix
  inverse_sigma <- solve(sigma)
  temp <- t(design_matrix) %*% inverse_sigma %*% design_matrix
  
  # Initialize the coefficient vector
  coeff = rep(0, n_coeffs)
  coeff[sel_coeff] = 1
  
  # Calculate the denominator for the test statistic
  den <- (t(coeff) %*% solve(temp) %*% coeff)
  
  # Define the integrand for the numerator of the test statistic
  integrand <- function(w) (eval.fd(w, mod$betaestlist[[sel_coeff]]$fd[1, 1]))^2
  
  # Calculate the numerator of the test statistic
  num <- integrate(integrand, lower = min(t_grid), upper = max(t_grid))$value
  
  # Calculate the observed test statistic
  T0 <- num / den
  
  # Initialize vector for storing bootstrap test statistics
  T.perm.res <- numeric(B)
  n <- dim(res)[1]
  
  sum = 0
  set.seed(seed)
  
  # Perform bootstrap resampling
  for (perm in 1:B) {
    
    # Generate a random permutation
    permutation <- sample(1:n)
    
    # Apply the permutation to the residuals
    res.b <- as.matrix(res[permutation, ])
    res_correlati.b = as.matrix(sqrt_root_sigma) %*% t(res.b)
    
    # Generate a response function with the permuted residuals
    response.b <- as.matrix(y_true + t(res_correlati.b))
    
    # Create a functional data object
    data.fd.b <- Data2fd(y = response.b, argvals = t_grid, basisobj = basis, lambda = opt.lambda)
    
    # Fit the model with permuted data
    mod.b <- fRegress(y = data.fd.b, xfdlist = xlist, betalist = blist)
    
    # Calculate the numerator of the test statistic for bootstrap sample
    integrand <- function(w) (eval.fd(w, mod.b$betaestlist[[sel_coeff]]$fd[1, 1]))^2
    num <- integrate(integrand, lower = min(t_grid), upper = max(t_grid))$value
    
    # Calculate the bootstrap test statistic
    T.perm.res[perm] <- num / den
    
    # Count the number of bootstrap samples with test statistic greater than observed
    if (T.perm.res[perm] > T0) {
      sum = sum + 1
    }
  }
  
  # Calculate the p-value
  p_val = sum / B
  
  # Return the results
  return(list(t_stat = T.perm.res, t_obs = T0, p_val = p_val))
}

#' Functional Bootstrap Confidence Intervals
#'
#' This function performs functional bootstrap resampling to compute confidence intervals for functional coefficients.
#'
#' @param logperiods Logarithm of time periods.
#' @param basis Basis matrix.
#' @param mod The model object.
#' @param opt.lambda Optimal lambda.
#' @param blist List of functional coefficients.
#' @param t0 Starting time.
#' @param tN Ending time.
#' @param length.out Number of grid points.
#' @param B Number of bootstrap samples.
#' @param seed Seed for randomization.
#'
#' @return Matrix of bootstrap functional coefficients.
#'
#' @examples
#' result <- functional_bootstrap_CI(logperiods, basis, my_model, my_lambda, my_blist, 0, 10, 100, 1000, 123)
#'
functional_bootstrap_CI = function(t_grid, basis, mod, opt.lambda, xlist, blist, coeff, sqrt_root_sigma, names, B, seed){
  
  # Define the time grid
  length.out = length(t_grid)
  
  colors <- colorRampPalette(c("#FFC700", "#9D260C"))(B)
  
  # Extract the real y
  y.fd <- eval.fd(t_grid, mod.no_sof$yfdobj)
  
  # Compute the functional residuals
  res.fd <- eval.fd(t_grid, mod$yfdobj - mod$yhatfdobj)
  
  # Set seed for reproducibility
  set.seed(seed)
    
  # Initialize a vector for storing bootstrap CI
  T.boot.L <- matrix(, nrow = length.out, ncol=B)
  
  # Perform Bootstrap inference
  for (b in 1:B) {
    
    # Generate a random permutation
    permutation <- sample(1:length.out)
    
    # Apply the permutation to the residuals
    res.b <- as.matrix(res.fd[permutation, ])
    res_correlati.b = as.matrix(sqrt_root_sigma) %*% t(res.b)
    
    # Generate a response function with the permuted residuals
    response.b <- as.matrix(y.fd + t(res_correlati.b))
    
    # Create a functional data object
    data.fd.b <- Data2fd(y = response.b, argvals = t_grid, basisobj = basis, lambda = opt.lambda)
    
    # Fit the model with permuted data
    mod.b <- fRegress(y = data.fd.b, xfdlist = xlist, betalist = blist)
    
    # Store the bootstrap functional coefficients
    T.boot.L[, b] <- eval.fd(t_grid, mod.b$betaestlist[[coeff]]$fd)
    
  }
    
  # Calculate quantiles for confidence interval
  boxplot.fd(T.boot.L, col = colors, barcol = "#9D260C", outliercol = 'black', main = names[coeff], xlab='log(T)', ylab='coefficient')
}

