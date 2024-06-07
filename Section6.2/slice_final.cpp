#include <Rcpp.h>
using namespace Rcpp;

// Log-likelihood function for the Dirichlet distribution
// [[Rcpp::export]]
double loglikelihood(NumericMatrix Yaug, NumericVector alpha) {
  int n2 = Yaug.nrow();
  int k = Yaug.ncol();
  double loglikelihood = 0.0;
  
  double alpha_sum = sum(alpha);
  
  for (int i = 0; i < n2; i++) {
    loglikelihood += lgamma(alpha_sum);
    for (int j = 0; j < k; j++) {
      loglikelihood += (alpha[j] - 1) * log(Yaug(i, j)) - lgamma(alpha[j]);
    }
  }
  return loglikelihood;
}

// Log-prior function for the gamma distribution
// [[Rcpp::export]]
double logprior(NumericVector alpha, NumericVector shape_est, NumericVector scale_est) {
  int k = alpha.size();
  double logprior = 0.0;
  
  for (int i = 0; i < k; i++) {
    logprior += R::dgamma(alpha[i], shape_est[i], scale_est[i], true);
  }
  
  return logprior;
}

// Log-posterior function
// [[Rcpp::export]]
double logpost(NumericMatrix Yaug, NumericVector alpha, NumericVector shape_est, NumericVector scale_est, int K) {
  double logpost = K * loglikelihood(Yaug, alpha) + logprior(alpha, shape_est, scale_est);
  
  return logpost;
}

// Slice sampler function
// [[Rcpp::export]]
NumericVector slice_sampler(NumericVector x0, NumericMatrix Yaug, NumericVector shape_est, NumericVector scale_est, int K) {
  int k = x0.size();
  NumericVector x1 = NumericVector(k, 0.0);
  NumericVector w = NumericVector(k, 1.0);
  
  double fx0 = log(R::runif(0, 1)) + logpost(Yaug, x0, shape_est, scale_est, K);
  
  NumericVector u = Rcpp::runif(k);
  NumericVector LL = x0 - w * u;
  NumericVector RR = LL + w;
  LL = ifelse(LL < 0.0, 0.0, LL);
  RR = ifelse(RR < 0.0, 0.0, RR);
  
  for (int i1 = 0; i1 < k; i1++) {
    x1[i1] = R::runif(LL[i1], RR[i1]);
  }
  int z = 0;
  
  while (fx0 >= logpost(Yaug, x1, shape_est, scale_est, K)) {
    LL = ifelse(x1 < x0, x1, LL);
    RR = ifelse(x1 >= x0, x1, RR);
    for (int i1 = 0; i1 < k; i1++) {
      x1[i1] = R::runif(LL[i1], RR[i1]);
    }
  }
  return x1;
}

// Define the log-density function for the Laplace distribution
double log_dlaplace(double x, double mu, double b) {
  return -std::abs(x - mu) / b - std::log(2 * b);
}

// [[Rcpp::export]]
double sum_log_dlaplace(NumericVector x, NumericVector mu, double b) {
  double aux = 0.0;
  for (int i = 0; i < x.size(); i++) {
    aux += log_dlaplace(x[i], mu[i], b);
  }
  return aux;
}

// [[Rcpp::export]]
List DPstep(NumericMatrix Yaug, NumericVector alpha, NumericVector S0,
            NumericVector nSS2, NumericVector shape_est,
            NumericVector scale_est, int n2, int K, double a, double b) {
  int k = Yaug.ncol();
  NumericVector ystar(k);
  NumericVector S0star(k);
  
  double test0;
  double paccept;
  
  for (int i = 0; i < n2 / K; i++) {
    
    for (int j = 0; j < k; j++) {
      ystar[j] = R::rgamma(alpha[j], 1.0);
    }
    ystar = ystar / sum(ystar);
    NumericVector test = - K * Rcpp::log(Rcpp::ifelse(Yaug(i, _) < a, a, Yaug(i, _))) / n2 +
      K * Rcpp::log(Rcpp::ifelse(ystar < a, a, ystar)) / n2;
    
    for (int j = 0; j < k; j++) {
      S0star[j] = S0[j] + test[j];
    }
    // Rcpp::Rcout << "a " << a << std::endl; // Print paccept value
    //Rcpp::Rcout << "K " << K << std::endl; // Print paccept value
    //Rcpp::Rcout << "n2 " << n2 << std::endl; // Print paccept value

    paccept = std::min(exp(sum_log_dlaplace(nSS2, S0star, b) - sum_log_dlaplace(nSS2, S0, b)), 1.0);
    // Rcpp::Rcout << "S0: " << S0 << std::endl; // Print paccept value
    // Rcpp::Rcout << "S0star: " << S0star << std::endl; // Print paccept value

    if (R::rbinom(1, paccept) == 1) {
      Yaug(i, _) = ystar;
      S0 = Rcpp::clone(S0star);
    }
  }
  
  // Create a list to return Yaug and S0
  List result;
  result["Yaug"] = Yaug;
  result["S0"] = S0;
  
  return result;
}

// Calculate the mean of each column in a NumericMatrix
NumericVector colMeansCpp(NumericMatrix mat) {
  int ncol = mat.ncol();
  int nrow = mat.nrow();
  NumericVector means(ncol);
  
  for (int j = 0; j < ncol; j++) {
    double colSum = 0.0;
    for (int i = 0; i < nrow; i++) {
      colSum += mat(i, j);
    }
    means[j] = colSum / nrow;
  }
  
  return means;
}

// Calculate the mean of each column in a NumericMatrix
NumericVector colSumsCpp(NumericMatrix mat) {
  int ncol = mat.ncol();
  int nrow = mat.nrow();
  NumericVector sums(ncol);
  
  for (int j = 0; j < ncol; j++) {
    double colSum = 0.0;
    for (int i = 0; i < nrow; i++) {
      colSum += mat(i, j);
    }
    sums[j] = colSum;
  }
  
  return sums;
}


//Gibbs sampler
// [[Rcpp::export]]
List runGibbsSampler(int nsave, NumericMatrix Yaug, NumericVector alpha,
                     NumericVector S0, NumericVector nSS2,
                     NumericVector shape_est, NumericVector scale_est,
                     int n2, int K, double a, double b, int ch) {
  
  int ncol = Yaug.ncol();
  NumericMatrix ALPHA3(nsave, ncol);
  NumericMatrix SS0(nsave, ncol);
  
  NumericMatrix log_Yaug(n2/K, ncol);
  for (int i = 0; i < n2 / K; i++) {
    for (int j = 0; j < ncol; j++) {
      log_Yaug(i, j) = Yaug(i, j) < a ? std::log(a) : std::log(Yaug(i, j));;
    }
  }
  S0 = colMeansCpp(log_Yaug);
  // Rcpp::Rcout << "S0 v2: " << S0 << std::endl; // Print paccept value
  
  
  for (int J = 0; J < nsave; J++) {
    alpha = slice_sampler(alpha, Yaug, shape_est, scale_est, K);
    ALPHA3(J, _) = alpha;
    
    List result = DPstep(Yaug, alpha, S0, nSS2, shape_est, scale_est, n2, K, a, b);
    Yaug = as<NumericMatrix>(result["Yaug"]);
    
    S0 = as<NumericVector>(result["S0"]);

    SS0(J, _) = S0;
    
    if ((J+1) % 1000 == 0)
      Rcpp::Rcout << J+1 << "  chain " << ch << std::endl;
  }
  
  return List::create(Named("ALPHA3") = ALPHA3,
                      Named("SS0") = SS0);
}

//Gibbs sampler with no privacy
// [[Rcpp::export]]
NumericMatrix runGibbsSampler_noprivate(int nsave, NumericMatrix Yaug, NumericVector alpha,
                                        NumericVector shape_est, NumericVector scale_est,
                                        int K, int ch) {
  
  int ncol = Yaug.ncol();
  NumericMatrix ALPHA3(nsave, ncol);
  
  for (int J = 0; J < nsave; J++) {
    alpha = slice_sampler(alpha, Yaug, shape_est, scale_est, K);
    ALPHA3(J, _) = alpha;
    
    if ((J+1) % 1000 == 0)
      Rcpp::Rcout << J+1 << "  chain " << ch << std::endl;
  }
  
  return ALPHA3;
}


// Log density of multivariate normal distribution with provided inverse of Sigma and determinant of Sigma
// [[Rcpp::export]]
double log_dmvn(const NumericVector& x, const NumericVector& mu, 
                const NumericMatrix& invSigma, double detSigma) {
  int n = x.size();
  
  // Calculate the centered data vector (x - mu)
  NumericVector centered_x = x - mu;
  
  // Calculate the quadratic term -0.5 * (x - mu)^T * invSigma * (x - mu)
  double quad_term = 0.0;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      quad_term += -0.5 * centered_x[i] * invSigma(i, j) * centered_x[j];
    }
  }
  
  // Calculate the constant term -0.5 * n * log(2*pi) - 0.5 * log(det(Sigma))
  double const_term = -0.5 * n * log(2.0 * M_PI) - 0.5 * log(detSigma);
  
  // Calculate the log density value
  double log_density = quad_term + const_term;
  
  return log_density;
}

// Define the log-density function for the Gaussian copula
// [[Rcpp::export]]
double log_gaussian_copula_with_gamma_marginals(NumericVector alpha,
                                                NumericVector shape, NumericVector scale,
                                                const NumericMatrix& invSigma, double detSigma) {
  int k = alpha.size(); // Dimension of the copula
  NumericVector u = NumericVector(k, 0.0);
  NumericVector z = NumericVector(k, 0.0);
  NumericVector mu0(k, 0.0);
  NumericMatrix Sigma0(k, k);
  for (int i = 0; i < k; i++) {
    Sigma0(i, i) = 1.0;
  }
  
  // double loglik_copula = 0.0;
  double loglik_marginals = 0.0;
  
  // Calculate the log likelihood of the marginal gamma distributions
  for (int j = 0; j < k; j++) {
    u[j] = R::pgamma(alpha[j], shape[j], scale[j], true, false);
    z[j] = R::qnorm(u[j], 0.0, 1.0, true, false);
    loglik_marginals += R::dgamma(alpha[j], shape[j], scale[j], true);//log_gamma_density(alpha[j], shape[j], scale[j]);
  }
  // Rcpp::Rcout << "loglik_marginals " << loglik_marginals << std::endl; // Print paccept value
  // Rcpp::Rcout << "u " << u << std::endl; // Print paccept value
  // Rcpp::Rcout << "z " << z << std::endl; // Print paccept value
  
  double loglik_copula = log_dmvn(z, mu0, invSigma, detSigma) 
    -  log_dmvn(z, mu0, Sigma0, 1.0);
  // Rcpp::Rcout << "loglik_copula " << loglik_copula << std::endl; // Print paccept value
  
  double loglik_total = loglik_copula + loglik_marginals;

  
  // Check if loglik_total is NaN
  if (std::isnan(loglik_total)) {
    loglik_total = logprior(alpha, shape, scale);
  }
  
  return loglik_total;
}

// Posterior copula
double logpost_copula(NumericMatrix Yaug, NumericVector alpha, 
                      NumericVector shape_est, NumericVector scale_est, 
                      const NumericMatrix& invSigma, double detSigma, int K) {
  double logpost = K * loglikelihood(Yaug, alpha) + 
    log_gaussian_copula_with_gamma_marginals(alpha, shape_est, scale_est, invSigma, detSigma);
  
  return logpost;
}

// Slice sampler function WITH COPULA
// [[Rcpp::export]]
NumericVector slice_sampler_copula(
    NumericVector x0, NumericMatrix Yaug, 
    NumericVector shape_est, NumericVector scale_est, 
    const NumericMatrix& invSigma, double detSigma, int K) {
  
  int k = x0.size();
  NumericVector x1 = NumericVector(k, 0.0);
  NumericVector w = NumericVector(k, 1.0);
  
  double fx0 = log(R::runif(0, 1)) + 
    logpost_copula(Yaug, x0, shape_est, scale_est, invSigma, detSigma, K);
  
  NumericVector u = Rcpp::runif(k);
  NumericVector LL = x0 - w * u;
  NumericVector RR = LL + w;
  LL = ifelse(LL < 0, 0, LL);
  RR = ifelse(RR < 0, 0, RR);
  
  for (int i1 = 0; i1 < k; i1++) {
    x1[i1] = R::runif(LL[i1], RR[i1]);
  }
  int z = 0;
  
  while (fx0 >= logpost_copula(Yaug, x1, shape_est, scale_est, invSigma, detSigma, K)) {
    LL = ifelse(x1 < x0, x1, LL);
    RR = ifelse(x1 >= x0, x1, RR);
    for (int i1 = 0; i1 < k; i1++) {
      x1[i1] = R::runif(LL[i1], RR[i1]);
    }
  }
  return x1;
}



//Gibbs sampler with copula
// [[Rcpp::export]]
List runGibbsSampler_copula(int nsave, NumericMatrix Yaug, NumericVector alpha,
                            NumericVector S0, NumericVector nSS2,
                            NumericVector shape_est, NumericVector scale_est, 
                            const NumericMatrix& invSigma, double detSigma,
                            int n2, int K, double a, double b, int ch) {
  
  int ncol = Yaug.ncol();
  NumericMatrix ALPHA3(nsave, ncol);
  NumericMatrix SS0(nsave, ncol);
  
  NumericMatrix log_Yaug(n2/K, ncol);
  for (int i = 0; i < n2 / K; i++) {
    for (int j = 0; j < ncol; j++) {
      log_Yaug(i, j) = Yaug(i, j) < a ? std::log(a) : std::log(Yaug(i, j));;
    }
  }
  S0 = colMeansCpp(log_Yaug);
  // Rcpp::Rcout << "S0 v2: " << S0 << std::endl; // Print paccept value
  
  for (int J = 0; J < nsave; J++) {
    alpha = slice_sampler_copula(alpha, Yaug, shape_est, scale_est, invSigma, detSigma, K);
    ALPHA3(J, _) = alpha;
    
    List result = DPstep(Yaug, alpha, S0, nSS2, shape_est, scale_est, n2, K, a, b);
    Yaug = as<NumericMatrix>(result["Yaug"]);
    S0 = as<NumericVector>(result["S0"]);
    
    SS0(J, _) = S0;
    
    if ((J+1) % 1000 == 0)
      Rcpp::Rcout << J+1 << "  chain " << ch << std::endl;
  }
  
  return List::create(Named("ALPHA3") = ALPHA3,
                      Named("SS0") = SS0);
}

//Unbounded DPstep
// [[Rcpp::export]]
List UnboundedDPstep(NumericMatrix Yaug, NumericVector alpha, NumericVector S0,
            NumericVector nSS2, int n_dp, NumericVector shape_est,
            NumericVector scale_est, int n, int K, double a, double b, double b_n) {
  int k = Yaug.ncol();
  NumericVector ystar(k);
  NumericVector S0star(k);
  
  double test0;
  double paccept;
  
  for (int i = 0; i < n / K; i++) {
    
    for (int j = 0; j < k; j++) {
      ystar[j] = R::rgamma(alpha[j], 1.0);
    }
    ystar = ystar / sum(ystar);
    NumericVector test = -Rcpp::log(Rcpp::ifelse(Yaug(i, _) < a, a, Yaug(i, _))) +
      Rcpp::log(Rcpp::ifelse(ystar < a, a, ystar));
    
    for (int j = 0; j < k; j++) {
      S0star[j] = S0[j] + test[j];
    }
    // Rcpp::Rcout << "a " << a << std::endl; // Print paccept value
    //Rcpp::Rcout << "K " << K << std::endl; // Print paccept value
    //Rcpp::Rcout << "n " << n << std::endl; // Print paccept value
    
    paccept = std::min(exp(sum_log_dlaplace(nSS2, S0star, b) - sum_log_dlaplace(nSS2, S0, b)), 1.0);
    // Rcpp::Rcout << "S0: " << S0 << std::endl; // Print paccept value
    // Rcpp::Rcout << "S0star: " << S0star << std::endl; // Print paccept value
    
    if (R::rbinom(1, paccept) == 1) {
      Yaug(i, _) = ystar;
      S0 = Rcpp::clone(S0star);
    }
  }
  
  // if (true) {
    // RJMCMC
    // Rcpp::Rcout << "n_before: " << n << std::endl; // Print paccept value
    // Sample a random value that is either +1 or -1
    int random_value = (rand() % 2 == 0) ? 1 : -1;
    int n_prop = n + random_value;
  if (true) {
      // n+1
    if (random_value == 1){
      
      for (int j = 0; j < k; j++) {
        ystar[j] = R::rgamma(alpha[j], 1.0);
      }
      
      ystar = ystar / sum(ystar);
      NumericVector test = Rcpp::log(Rcpp::ifelse(ystar < a, a, ystar));
      
      for (int j = 0; j < k; j++) {
        S0star[j] = S0[j] + test[j];
      }
      paccept = std::min(exp(sum_log_dlaplace(nSS2, S0star, b) +  log_dlaplace(n_dp, n_prop, b_n)
                               - sum_log_dlaplace(nSS2, S0, b) - log_dlaplace(n_dp, n, b_n)), 1.0);
      
      // Rcpp::Rcout << "sum_log_dlaplace(nSS2, S0star, b): " << sum_log_dlaplace(nSS2, S0star, b) << std::endl; // Print paccept value
      // Rcpp::Rcout << "log_dlaplace(n_dp, n_prop, 1): " << log_dlaplace(n_dp, n_prop, 1) << std::endl; // Print paccept value
      // Rcpp::Rcout << "sum_log_dlaplace(nSS2, S0, b): " << sum_log_dlaplace(nSS2, S0, b) << std::endl; // Print paccept value
      // Rcpp::Rcout << "log_dlaplace(n_dp, n, 1)): " << log_dlaplace(n_dp, n, 1) << std::endl; // Print paccept value
      
      // paccept = std::min(exp(log_dlaplace(n_dp, n_prop, 1) - log_dlaplace(n_dp, n, 1)), 1.0);
      
      // Rcpp::Rcout << "++++++++ 1" << std::endl; // Print paccept value
      // Rcpp::Rcout << "n before " << n << std::endl; // Print paccept value
      // Rcpp::Rcout << "Yaug(n-1, _) before " << Yaug(n-1, 0) << std::endl; 
      // Rcpp::Rcout << "Yaug(n, _) before " << Yaug(n, 0) << std::endl; 
      // Rcpp::Rcout << "Yaug(n+1, _) before " << Yaug(n+1, 0) << std::endl; 

      if (R::rbinom(1, paccept) == 1) {
        n = n_prop;
        Yaug(n-1, _) = ystar;
        S0 = Rcpp::clone(S0star);
        // Rcpp::Rcout << "n_prop " << n_prop << std::endl; // Print paccept value
      }
      // Rcpp::Rcout << "n after" << n << std::endl; // Print paccept value
      // Rcpp::Rcout << "Yaug(n-1, _) " << Yaug(n-1, 0) << std::endl; 
      // Rcpp::Rcout << "Yaug(n, _) " << Yaug(n, 0) << std::endl; 
      // Rcpp::Rcout << "Yaug(n+1, _) " << Yaug(n+1, 0) << std::endl; 
    }
    // n-1
    if (random_value == -1){
      
      ystar = Yaug(n-1, _);
      
      //ystar = ystar / sum(ystar);
      NumericVector test = Rcpp::log(Rcpp::ifelse(ystar < a, a, ystar));
      
      for (int j = 0; j < k; j++) {
        S0star[j] = S0[j] - test[j];
      }
      paccept = std::min(exp(sum_log_dlaplace(nSS2, S0star, b) +  log_dlaplace(n_dp, n_prop, b_n)
                               - sum_log_dlaplace(nSS2, S0, b) - log_dlaplace(n_dp, n, b_n)), 1.0);
      
      // Rcpp::Rcout << "sum_log_dlaplace(nSS2, S0star, b): " << sum_log_dlaplace(nSS2, S0star, b) << std::endl; // Print paccept value
      // Rcpp::Rcout << "log_dlaplace(n_dp, n_prop, 1): " << log_dlaplace(n_dp, n_prop, 1) << std::endl; // Print paccept value
      // Rcpp::Rcout << "sum_log_dlaplace(nSS2, S0, b): " << sum_log_dlaplace(nSS2, S0, b) << std::endl; // Print paccept value
      // Rcpp::Rcout << "log_dlaplace(n_dp, n, 1)): " << log_dlaplace(n_dp, n, 1) << std::endl; // Print paccept value
      
      // paccept = std::min(exp(log_dlaplace(n_dp, n_prop, 1) - log_dlaplace(n_dp, n, 1)), 1.0);
      
      // Rcpp::Rcout << "-------- 1" << std::endl; // Print paccept value
      // Rcpp::Rcout << "n before " << n << std::endl; // Print paccept value
      // Rcpp::Rcout << "Yaug(n-1, _) before " << Yaug(n-1, 0) << std::endl; 
      // Rcpp::Rcout << "Yaug(n, _) before " << Yaug(n, 0) << std::endl; 
      // Rcpp::Rcout << "Yaug(n+1, _) before " << Yaug(n+1, 0) << std::endl; 
      
      
      if (R::rbinom(1, paccept) == 1) {
        n = n_prop;
        for (int j = 0; j < k; j++) {
          Yaug(n, j) = 1.0;
        }
        S0 = Rcpp::clone(S0star);
        // Rcpp::Rcout << "n_prop " << n_prop << std::endl; // Print paccept value
      }
      // Rcpp::Rcout << "n after " << n << std::endl; // Print paccept value
      // Rcpp::Rcout << "Yaug(n-1, _) " << Yaug(n-1, 0) << std::endl; 
      // Rcpp::Rcout << "Yaug(n, _) " << Yaug(n, 0) << std::endl; 
      // Rcpp::Rcout << "Yaug(n+1, _) " << Yaug(n+1, 0) << std::endl; 
      
    }
  }
  // Rcpp::Rcout << "random_value: " << random_value << std::endl; // Print paccept value
  // Rcpp::Rcout << "n_prop: " << n_prop << std::endl; // Print paccept value
  // Rcpp::Rcout << "paccept: " << paccept << std::endl; // Print paccept value
  // Rcpp::Rcout << "n_after: " << n << std::endl; // Print paccept value
  
  // Create a list to return Yaug and S0
  List result;
  result["Yaug"] = Yaug;
  result["S0"] = S0;
  result["n"] = n;
  
  return result;
}

//Gibbs sampler unbounded
// [[Rcpp::export]]
List runUnboundedGibbsSampler(int nsave, NumericMatrix Yaug, NumericVector alpha,
                     NumericVector S0, NumericVector nSS2, int n_dp,
                     NumericVector shape_est, NumericVector scale_est,
                     int n, int K, double a, double b, double b_n, int ch) {
  
  int ncol = Yaug.ncol();
  NumericMatrix ALPHA3(nsave, ncol);
  NumericMatrix SS0(nsave, ncol);
  NumericVector nn(nsave);
  
  NumericMatrix log_Yaug(n/K, ncol);
  for (int i = 0; i < n/ K; i++) {
    for (int j = 0; j < ncol; j++) {
      log_Yaug(i, j) = Yaug(i, j) < a ? std::log(a) : std::log(Yaug(i, j));;
    }
  }
  S0 = colSumsCpp(log_Yaug);
  // Rcpp::Rcout << "S0 v2: " << S0 << std::endl; // Print paccept value
  
  
  for (int J = 0; J < nsave; J++) {
    
    // Create a new matrix with the first n rows
    NumericMatrix YaugSubset(n, ncol);
    
    // Copy the first n rows from Yaug to YaugSubset
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < ncol; j++) {
        YaugSubset(i, j) = Yaug(i, j);
      }
    }
    //Rcpp::Rcout << "YaugSubset " << YaugSubset(n-1,0) << std::endl; // Print paccept value
    alpha = slice_sampler(alpha, YaugSubset, shape_est, scale_est, K);
    ALPHA3(J, _) = alpha;
    
    List result = UnboundedDPstep(Yaug, alpha, S0, nSS2, n_dp,
                                  shape_est, scale_est, n, K, a, b, b_n);
    Yaug = as<NumericMatrix>(result["Yaug"]);
    S0 = as<NumericVector>(result["S0"]);
    n = Rcpp::as<int>(result["n"]);
    //Rcpp::Rcout << "n " << n << std::endl; // Print paccept value
    
    SS0(J, _) = S0;
    nn(J) = n;
    
    
    if ((J+1) % 1000 == 0)
      Rcpp::Rcout << J+1 << "  chain " << ch << std::endl;
  }
  
  return List::create(Named("ALPHA3") = ALPHA3,
                      Named("SS0") = SS0,
                      Named("n") = nn);
}
