################################################################################
# Functions for DP parametric bootstrap
################################################################################
# Function for calculating the derivative of the digamma function
sirt_digamma1 <- function(x, h=1e-3)
{
  ( digamma(x+h) - digamma(x-h) ) / (2*h)
}

# Function for finding the maximum likelihood estimate of alpha using the mle_dirichlet algorithm
mle_dirichlet = function(f_D, alpha, maxit, convcrit=1e-05, oldfac=.3){
  # convcrit: Convergence criterion
  # maxit: Maximum number of iterations
  # oldfac: Convergence acceleration factor. It must be a parameter between 0 and 1.
  conv = 1
  iter = 1
  while( ( conv > convcrit ) & (iter < maxit) ){
    alpha0 = alpha
    g = digamma( sum(alpha ) ) - digamma(alpha) + f_D
    z = sirt_digamma1( sum(alpha ))
    if (z==0){
      alpha='try-error'
      break
    }
    H = diag( -sirt_digamma1( alpha ) ) + z
    alpha = try(alpha0 - Matrix::solve(H, g ))
    if ('try-error' %in% class(alpha)){
      alpha='try-error'
      break
    }else{
      alpha[ alpha < 0 ] = 1e-10
      alpha = alpha0 + oldfac*( alpha - alpha0 )
      conv = max( abs( alpha0 - alpha ) )
      iter = iter+1
    }
  }
  return(alpha)
}


#-- DPMLEBoots_iter -- This function generates a single draw using DPMLEBoots algorithm
# nS: DP sufficient statistic
# alpha0: starting alpha
# b: variance Laplace mechanism
DPMLEBoots_iter = function(nS, alpha_initial, epsilon, n)
{
  # Sampling using a truncated Laplace distribution that guarantee
  # nSS = nS + error is nonpositive
  b = -k*log(a)/(n*epsilon)
  error = sapply(1:k, function(j) rtrunc(1, "laplace", a=-Inf, b=-nS[j], location=0, scale=b))
  nSS0 = nS + error
  # nSS = nS + error is in the range of S_0 (original sufficient statistic)
  crit = sum(exp(nSS0))>=0.999
  while(crit > 0)
  {
    error = sapply(1:k, function(j) rtrunc(1, "laplace", a=-Inf, b=-nS[j], location=0, scale=b))
    error = rlaplace(k,0,b)
    nSS0 = nS + error
    crit = sum(exp(nSS0))>=0.999
  }
  alpha = mle_dirichlet(f_D=nSS0, alpha=alpha_initial, maxit = 100)
  
  if(length(alpha)==length(nSS0))
  {
    Ytmp = rdirichlet(n, alpha = alpha)
    return(dirichlet.mle(Ytmp)$alpha)
  }
}
##-- DPMLEBoots -- This function generates a multiple draws using DPMLEBoots_iter
## algorithm
# N: Number of draws
# nS: DP sufficient statistic
# alpha0: starting alpha
# b: variance Laplace mechanism
DPMLEBoots = function(N, nS, alpha_initial, epsilon,n)
{
  #start_time <- Sys.time()
  # Generate N draws
  ALPHA = lapply(1:N,function(i) try(DPMLEBoots_iter(nS, alpha_initial, epsilon,n)))
  Ind = which(!(simplify2array(lapply(ALPHA, is.vector))))
  while(length(Ind)>0)
  {
    ALPHAtmp = lapply(1:length(Ind),function(i) try(DPMLEBoots(i, nS, alpha_initial, epsilon,n)))
    for(i in 1:length(Ind)) ALPHA[[Ind[i]]] = ALPHAtmp[[i]]
    Ind = which(!(simplify2array(lapply(ALPHA, is.vector))))
  }
  ALPHA = t(simplify2array(ALPHA))
  #end_time <- Sys.time()
  #list(ALPHA = ALPHA,  time = end_time - start_time)
  return(ALPHA)
}

################################################################################
# Function that helps to select the initial alpha
################################################################################
initial_select <- function(i, n, nSS, b) {
  repeat{
    b=-k*log(a)/(n*epsilon_ss)
    auxS0 <- nSS+rlaplace(k,0,b)
    crit = sum(exp(auxS0))>=0.999
    while(crit > 0)
    {
      auxS0 = nSS+rlaplace(k,0,b)
      crit = sum(exp(auxS0))>=0.999
    }
    alpha_initial = mle_dirichlet(f_D=auxS0,alpha=runif(k,0.1,20), maxit = 100)
    if (is.numeric(alpha_initial)&(all(alpha_initial<1e5))){
      break
    }
  }
  return(alpha_initial=alpha_initial)
}