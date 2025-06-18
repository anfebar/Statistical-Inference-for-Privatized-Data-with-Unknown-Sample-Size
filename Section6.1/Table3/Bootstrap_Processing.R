#############
# Load required package and read results
#############
library(xtable)

# Read bootstrap output
output <- read.csv("bootstrap_outputs.csv")

# Drop the first column (usually row numbers from write.csv)
output <- output[,-1]

# Rename columns
names(output) <- c(
  "epn", "eps", 
  "beta0_m", "beta1_m", "beta2_m", 
  "beta0_v", "beta1_v", "beta2_v", 
  "n_dp", "Rn_dp", "inv_tau", "iteration"
)

# Add column: max(n_dp, p + 1)
p <- 3
output$n_dp_max <- apply(cbind(output$n_dp, p + 1), 1, max)

#############
# Define trimmed mean function - trim = 0.01 removes the maximum value 
# since the mean is computed over 100 values
#############
trim_mean <- function(x, trim = 0.01) {
  n <- length(x)
  cutoff <- floor(n * (1 - trim))
  mean(sort(x)[1:cutoff])
}

#############
# Aggregate trimmed means by (epn, eps) for each eps level
#############
tables <- list()

tables[["1_estimates"]] <- 
  t(aggregate(cbind(beta0_m, beta0_v, beta1_m, beta1_v,
                    beta2_m, beta2_v, inv_tau, n_dp, n_dp_max, Rn_dp) 
              ~ epn + eps, 
              FUN = trim_mean, 
              data = output[output$eps == 1, ]))[-c(1,2),]

tables[["0.1_estimates"]] <- 
  t(aggregate(cbind(beta0_m, beta0_v, beta1_m, beta1_v,
                    beta2_m, beta2_v, inv_tau, n_dp, n_dp_max, Rn_dp) 
              ~ epn + eps, 
              FUN = trim_mean, 
              data = output[output$eps == 0.1, ]))[-c(1,2),]

tables[["10^7_estimates"]] <- 
  t(aggregate(cbind(beta0_m, beta0_v, beta1_m, beta1_v,
                    beta2_m, beta2_v, inv_tau, n_dp, n_dp_max, Rn_dp) 
              ~ epn + eps, 
              FUN = trim_mean, 
              data = output[output$eps == 1e7, ]))[-c(1,2),]

#############
# Display tables using xtable
#############
xtable(tables[["0.1_estimates"]], digits = 3, align = "lrrrrrr")
xtable(tables[["1_estimates"]], digits = 3, align = "lrrrrrr")
xtable(tables[["10^7_estimates"]], digits = 3, align = "lrr")
