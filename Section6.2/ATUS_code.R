library(coda)
library(LaplacesDemon)
library(parallel)
library(Rcpp)

source("aux.R") # File with auxiliar functions
sourceCpp("slice_final.cpp") # File with MCMC implementations

for(iteration in 1:10)
{
  # MCMC length, number of chains, and retained iterations
  nsave <- 2000
  nch <- 100
  draw_seq <- 1000+(1:10)*100
  

  for(epsilon_ss in c(1,10))  # Privacy budget for s
  {
    # Read data
    xF <- read.csv("female.csv")
    xM <- read.csv("male.csv")
    
    x_true <- rbind(xF,xM)[sample(1:(nrow(xF)+nrow(xM)),6656),]
    a <- min(x_true) # Find truncation level
    k <- dim(x_true)[2]
    n <- dim(x_true)[1]
    
    
    #####################################################################
    ############################## BOUNDED DP ###########################
    #####################################################################
    set.seed(iteration)
    # Compute noisy summaries
    b <- -k*log(a)/(n*epsilon_ss)
    noise_SS <- rlaplace(k, 0, b)
    nSS <- colMeans(log(a*(x_true<=a)+x_true*(x_true>a))) + noise_SS
    
    print("Running Bounded DPMCMCp1")
    # Prior 1
    # Prior parameters
    shape_prior <- 1
    rate_prior <- 0.1  # mean = shape / rate, var = shape / (rate)^2
    scale_prior <- 1/rate_prior
    # DPMCMCp1
    Alpha0 <- sapply(1:nch, initial_select, n = n, nSS = nSS, b = b) # Find starting points
    fch <- function(ch)
    {
      set.seed(ch)
      alpha0 <- Alpha0[,ch]
      
      # initial values for mcmc
      Yaug <- rdirichlet(n,alpha0)
      S0 <- colMeans(log(ifelse(Yaug < a, a, Yaug)))
      K <- 1
      Yaug <- Yaug[rep(1:(n/K),1),]
      # Run MCMC
      fitmcmc <- runGibbsSampler(nsave, Yaug, alpha0,
                                 S0, nSS,
                                 shape_est = rep(shape_prior,k), 
                                 scale_est = rep(scale_prior,k),
                                 n, K, a, b, ch)
      fitmcmc[[1]][draw_seq,]
    } 
    alpha_BoundedDPMCMCp1 <- mclapply(1:nch, fch, mc.cores = 15, mc.preschedule = F)
    summary(mcmc.list(lapply(alpha_BoundedDPMCMCp1,mcmc)))[[1]]
    
    #####################################################################
    ############################## UNBOUNDED DP #########################
    #####################################################################
    for(iter in 1:10)
    {
      OUT <- list()
      J <- 0      
      for(epsilon_n in  c(.01, .1, 1, 10)) # Privacy budget for n
      {
        # Compute noisy summaries
        b <- -k*log(a)/epsilon_ss
        b_n <- 1/epsilon_n
        nSS <- colSums(log(a*(x_true<=a)+x_true*(x_true>a))) + n*noise_SS
        set.seed(iter)
        n_dp <- round(rlaplace(1, n, b_n))

        print("Running UnBounded DPMCMCp1")
        # Prior 1
        # Prior parameters
        shape_prior <- 1
        rate_prior <- 0.1  # mean = shape / rate, var = shape / (rate)^2
        scale_prior <- 1/rate_prior
        # DPMCMCp1
        fch <- function(ch)
        {
          set.seed(ch)
          nn_dp <- round(qlaplace(ch/(nch+1), n_dp, b_n))
          n_dp <- nn_dp
          set.seed(ch)
          alpha0 <- Alpha0[,ch] 
          
          # Initial values for mcmc
          Yaug <- rdirichlet(n_dp,alpha0)
          S0 <- colSums(log(ifelse(Yaug < a, a, Yaug)))
          K <- 1
          Yaug <- rbind(Yaug, matrix(1, nrow = 2*n_dp, ncol = ncol(Yaug)))
          # Run MCMC
          fitmcmc <- runUnboundedGibbsSampler(nsave, Yaug, alpha0,
                                              S0, nSS, n_dp,
                                              shape_est = rep(shape_prior,k), 
                                              scale_est = rep(scale_prior,k),
                                              n = n_dp, K, a, b, b_n, ch)
          list(fitmcmc[[1]][draw_seq,], fitmcmc[[3]][draw_seq])
        }  
        out <- mclapply(1:nch, fch, mc.cores = 15, mc.preschedule = F)

        J <- J + 1
        OUT[[J]] <- out
      }
      
      save(alpha_BoundedDPMCMCp1, OUT, 
           file = paste("output/epsilonSS",epsilon_ss,
                        "_iteration",iteration,
                        "_iter",iter,".RData",sep = ""))
      
    }
    
  }
}
