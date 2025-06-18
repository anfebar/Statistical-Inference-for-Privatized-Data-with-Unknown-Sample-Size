library(MASS)        ### For multivariate normal simulation
library(parallel)    ### For parallel computing
library(matrixcalc)  ### For checking if a matrix is positive definite

#---------------------------------------------
# Function to generate simulated data
#---------------------------------------------
rinit_data <- function(N, params) {
  
  ### Generate X from multivariate normal distribution
  x <- mvrnorm(n = N, mu = params$mu, Sigma = params$phi_inv)
  X <- cbind(1, x)
  
  ### Generate Y from linear model with Gaussian noise
  Y <- as.numeric(X %*% params$beta + rnorm(N, 0, sd = sqrt(1 / params$tau)))
  
  ### Return data list including sample size
  return(list(N = N, X = X, Y = Y))
}

#---------------------------------------------
# Load DP bootstrap functions
#---------------------------------------------
# This file includes the implementation of the parametric 
# bootstrap by Ferrando et al. (2022), adapted by Barrientos et al. (2024)
source("DPBootsLMFunctions.R")

#---------------------------------------------
# Define design matrix and response
#---------------------------------------------
### Variable names (last one is the response)
Names <- c("X1", "X2", "Y")

### Type of each variable
Type <- c(X1 = "numeric", 
          X2 = "numeric", 
          Y  = "numeric")

### Bounds for each variable
Bounds <- list(X1 = c(-5, 5), 
               X2 = c(-5, 5), 
               Y  = c(-5, 5))

#---------------------------------------------
# Define values for privacy parameters
#---------------------------------------------
### Rows of Eps contain (ep_N, eps)
Eps <- rbind(
  cbind(c(1e7, 0.0001), 1e7),
  cbind(c(1e7, 1, 0.1, 0.01, 0.001, 0.0001), 1),
  cbind(c(1e7, 1, 0.1, 0.01, 0.001, 0.0001), 0.1)
)
colnames(Eps) <- c("ep_N", "eps")

#---------------------------------------------
# Set true parameter values and sample size
#---------------------------------------------
N <- 1000

params_true <- list(
  beta    = c(0, -1, 1),
  tau     = 1,
  mu      = c(1, -1),
  phi     = diag(2),
  phi_inv = diag(2)
)


#-------------------------------------------------------------
# f_Eps evaluates the DP bootstrap for a given (ep_N, eps) pair.
# It simulates data, runs the DP bootstrap algorithm, and returns
# summary statistics including mean and variance of beta estimates,
# privatized sample size, and estimated error variance.
#-------------------------------------------------------------
f_Eps <- function(j, iteration) {
  
  ep_N <- Eps[j, 1]
  eps  <- Eps[j, 2]
  
  ### Fix seed for generating the true (sensitive) dataset
  set.seed(0)
  data <- rinit_data(N, params_true)
  
  ### Use iteration as seed for privatized statistics
  set.seed(iteration)
  A_data <- data.frame(data$X[,-1], Y = data$Y)
  
  ### Run the DP bootstrap algorithm
  out <- DpBoots.lm(A_data, Names, Type, Bounds, Ref = NULL, 
                    epsilon = list(eps = eps, ep_N = ep_N), 
                    alpha = 0.05)
  
  if (j == nrow(Eps)) {
    cat("iteration =", iteration, "\n")
  }
  
  ### Return summary statistics
  c(ep_n = ep_N, eps = eps, 
    beta_m = apply(out$bootsBeta, 2, mean),
    beta_var = apply(out$bootsBeta, 2, var),
    n_dp = out$NAtA[1, 1],
    Rn_dp = out$Rndp,
    sigma2 = out$hat.sigma2,
    iteration = iteration)
}

#-------------------------------------------------------------
# Run DP bootstrap simulations across multiple iterations and 
# multiple (ep_N, eps) configurations using parallel processing.
#-------------------------------------------------------------

out <- do.call(rbind, 
               mclapply(1:100, function(iteration) {
                 
                 ### For each iteration, evaluate all (ep_N, eps) pairs in Eps
                 do.call(rbind, 
                         mclapply(1:nrow(Eps), f_Eps, 
                                  iteration = iteration, 
                                  mc.cores = 3,               # Use 3 cores for inner loop
                                  mc.preschedule = FALSE))   # Avoid pre-scheduling to balance load better
                 
               },
               mc.cores = 5,                 # Use 5 cores for outer loop (iterations)
               mc.preschedule = FALSE))      # Again, disable prescheduling for load balancing

#-------------------------------------------------------------
# Save the results to a CSV file
#-------------------------------------------------------------
write.csv(out, file = "bootstrap_outputs.csv")
