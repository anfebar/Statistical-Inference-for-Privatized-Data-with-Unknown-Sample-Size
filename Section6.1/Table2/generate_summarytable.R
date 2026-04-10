library(parallel)
library(MASS)
library(VGAM)

source("LRMCEM_unbounded_helper.R")

p <- 2
stepsize <- 0.01
eps <- 1
N <- 1000
num_cycles <- 10000
iter_vec <- 1:100
ep_N_vec <- c(0.001, 0.01, 0.1, 1, 10, Inf)
mc_cores <- 10

hyperparams <- list(
  p = 2,
  m = rep(0, p + 1),
  V_inv = diag(p + 1),
  a = 2,
  b = 2,
  theta = rep(0, p),
  Sigma_inv = diag(p),
  d = 2,
  W_inv = diag(p)
)

run_one_case <- function(job) {
  iteration <- job$iteration
  epn <- job$epn

  cat("starting iteration =", iteration, "ep_N =", epn, "\n")

  privacy_params <- list(
    eps = eps,
    ep_N = epn,
    deltaa = p^2 / 2 + (2.5) * p + 2,
    clbounds = c(-5, -5, 5, 5)
  )

  out <- run_mcem(
    iteration = iteration,
    num_cycles = num_cycles,
    N = N,
    ep_N = epn,
    eps = eps,
    N_guess = NA,
    hyperparams = hyperparams,
    privacy_params = privacy_params,
    swaps = 1,
    jumps = 1,
    fullUpdate = TRUE
  )

  data.frame(
    iter = iteration,
    eps = eps,
    epn = epn,
    b1 = out$state$params$beta[1],
    b2 = out$state$params$beta[2],
    b3 = out$state$params$beta[3],
    tau = out$state$params$tau,
    m1 = out$state$params$mu[1],
    m2 = out$state$params$mu[2]
  )
}

jobs <- expand.grid(
  iteration = iter_vec,
  epn = ep_N_vec,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

results <- do.call(
  rbind,
  mclapply(
    seq_len(nrow(jobs)),
    function(i) run_one_case(jobs[i, ]),
    mc.cores = mc_cores,
    mc.preschedule = FALSE
  )
)

write.table(
  results,
  file = "summarytable.txt",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)
